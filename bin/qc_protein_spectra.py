import argparse
from pathlib import Path
import re
import csv
import logging
from pyteomics import mzml

import numpy as np
import pandas as pd
from scipy import signal
from pybaselines import Baseline
from scipy import signal, interpolate


def find_integer_at_end(string):
    return int(re.search(r'\d+$', string).group()) if re.search(r'\d+$', string) else 'N/A'


def baseline_als(y, lam=1e5, p=0.01):
    """Asymmetric Least Squares Smoothing for baseline correction (MicrobeMS uses AsLS)."""
    baseline_fitter = Baseline(y)
    baseline, _ = baseline_fitter.asls(y, lam=lam, p=p)
    return baseline

def apply_high_mz_offset_correction(mz, intensity, threshold=10000):
    """
    Finds the mean intensity in the high m/z region (e.g., > 10,000) 
    and subtracts it from the entire spectrum.
    """
    # Identify the indices for the high m/z region
    high_mz_mask = mz > threshold
    
    if np.any(high_mz_mask):
        # Calculate the average intensity in that "quiet" region
        offset = np.mean(intensity[high_mz_mask])
        # Subtract the offset from the whole spectrum
        corrected_intensity = intensity - offset
        return corrected_intensity
    
    return intensity # Return as is if threshold is never reached

def calculate_resolving_power(mz, intensity, peak_indices):
    """
    Calculates m/z divided by FWHM for the top 10 most intense peaks
    using m/z difference rather than point-spacing.
    """
    if len(peak_indices) == 0: 
        return 0
        
    # Sort by intensity and take top 10
    top_peaks = sorted(peak_indices, key=lambda x: intensity[x], reverse=True)[:10]
    powers = []
    
    for idx in top_peaks:
        peak_mz = mz[idx]
        half_max = intensity[idx] / 2
        
        # Find the left bound of the FWHM
        left_idx = idx
        while left_idx > 0 and intensity[left_idx] > half_max:
            left_idx -= 1
            
        # Find the right bound of the FWHM
        right_idx = idx
        while right_idx < len(intensity) - 1 and intensity[right_idx] > half_max:
            right_idx += 1
            
        # Linear interpolation for more precision at the crossing point
        # Left m/z at half-max
        mz_l = np.interp(half_max, [intensity[left_idx], intensity[left_idx+1]], 
                                    [mz[left_idx], mz[left_idx+1]])
        # Right m/z at half-max
        mz_r = np.interp(half_max, [intensity[right_idx], intensity[right_idx-1]], 
                                    [mz[right_idx], mz[right_idx-1]])
        
        fwhm = mz_r - mz_l
        
        if fwhm > 0:
            powers.append(peak_mz / fwhm)
            
    return np.mean(powers) if powers else 0

def calculate_microbe_ms_style_noise_score(norm_intensity, window=99, poly=3):
    n = len(norm_intensity)
    if n < window:
        return 0.0, 6.0

    # 1. Smooth the vector to get the "baseline/trend"
    # We smooth first to define what 'noise' is deviating from
    smoothed = signal.savgol_filter(norm_intensity, window, poly)
    diff = norm_intensity - smoothed
    
    # 2. STRICT REMOVAL: Remove top 30% and bottom 2% of the residuals
    # This specifically removes the peak signals and spikes
    sorted_diff = np.sort(diff)
    lower_idx = int(n * 0.02)
    upper_idx = int(n * 0.70)
    
    trimmed_diff = sorted_diff[lower_idx:upper_idx]
    
    # 3. Determine Absolute Noise on remaining values
    if len(trimmed_diff) == 0:
        return 0.0, 6.0
        
    noise_std = np.std(trimmed_diff)
    
    # 4. Mapping [0.2, 6.0] -> [100, 0]
    noise_score = 100 * (6.0 - noise_std) / (6.0 - 0.2)
    return np.clip(noise_score, 0, 100), noise_std

def calculate_microbe_ms_style_peak_score(num_peaks):
    """
    Calculates the peak score using logarithmic binning limits as described.
    Limits: [14, 16, 19, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70]
    """
    # 14 bin limits define 15 possible bins
    limits = [14, 16, 19, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70]
    
    # If fewer than the first limit, score is 0
    if num_peaks < limits[0]:
        return 0.0
    
    # If equal or more than the highest limit (70), score is 100
    if num_peaks >= limits[-1]:
        return 100.0

    # Determine which bin the number falls into
    # We find how many limits the peak count has surpassed
    for i, limit in enumerate(limits):
        if num_peaks < limits[i+1]:
            # Each step is 100 / 14 (approx 7.1428)
            # This matches the example: 28 peaks falls at index 5 -> 5 * 7.14... = 42.85
            return (i + 1) * (100 / 14)
            
    return 100.0

def calculate_logistic_threshold(mz, noise_std, sensitivity_factor=3):
    """
    Creates a dynamic threshold curve using a generalized logistic function.
    Models higher sensitivity/noise at low m/z and lower at high m/z.
    """
    # Parameters for the logistic curve (typical for MALDI-ToF protein spectra)
    L = noise_std * 5   # Starting threshold (High at low m/z)
    K = noise_std * 2   # Ending threshold (Low at high m/z)
    x0 = 5000           # The inflection point (where the drop is steepest)
    k = 0.001           # Steepness of the curve
    
    # Generalized Logistic Formula: f(x) = K + (L - K) / (1 + exp(k * (x - x0)))
    threshold_curve = K + (L - K) / (1 + np.exp(k * (mz - x0)))
    
    return threshold_curve * sensitivity_factor

def microbe_ms_style_qc(mz, intensity, weights={'peaks': 0.55, 'noise': 0.30, 'baseline': 0.00, 'res': 0.15}):
    """Implements a QC scoring system inspired by MicrobeMS metrics for protein spectra.
    
    Metrics:
    1. Noise Quality (30% weight): Based on the standard deviation of the noise after smoothing.
    2. Baseline Quality (15% weight): Area under the baseline curve, with lower being better.
    3. Number of Peaks (40% weight): More peaks generally indicate better quality, but with diminishing returns.
    4. Resolving Power (15% weight): Average m/z divided by FWHM for the top 10 most intense peaks.
    
    Scoring:
    Each metric is scored from 0 to 100 based on thresholds derived from typical Microbe
    MS spectra, and then combined into a weighted total score. The final status is assigned as:
    - GREEN: Total Score > 45
    - YELLOW: 30 <= Total Score <= 45
    - RED: Total Score < 30
    """
    # MANDATORY CHECK: Weights must sum to 1.0
    if not np.isclose(sum(weights.values()), 1.0):
        raise ValueError("Error: The sum of weightings must equal 100% (1.0).")

    # STEP 1: Preprocessing in exact order
    raw_baseline = baseline_als(intensity)
    
    # "Cutting the spectra between two m/z values, usually between m/z 2000 and 13000" -- we will do 3k
    mask = (mz >= 3000) & (mz <= 13000)
    mz_c, int_c, base_c = mz[mask], intensity[mask], raw_baseline[mask]
    
    # Baseline subtraction and Normalization
    corrected = int_c - base_c
    norm_factor = np.sum(np.abs(corrected))
    norm_intensity = (corrected / norm_factor) * 100_000  # Scale to a more typical range for noise calculation
    
    # Offset correction (using high m/z region)
    high_mz_mask = mz_c > (mz_c.max() - 500)
    offset = np.mean(norm_intensity[high_mz_mask])
    final_intensity = norm_intensity - offset

    # STEP 2: Noise Score (Trimmed SD logic)
    noise_score, noise_std = calculate_microbe_ms_style_noise_score(final_intensity)

    ### DEBUG:
    print("noise_std:", noise_std, flush=True)
    
    # STEP 3: Baseline Score (Integral of normalized baseline)
    norm_baseline_curve = base_c / norm_factor
    baseline_area = np.trapz(norm_baseline_curve, mz_c)
    baseline_score = np.clip(100 * (1 - (baseline_area - 0.15) / (40 - 0.15)), 0, 100)

    # STEP 4: Peaks Score (Log Binning)
    # Note: Description mentions a 'logistic threshold function' here
    # We will use the simplified noise-based detection for now.
    threshold_curve = calculate_logistic_threshold(mz_c, noise_std)
    
    # Find peaks that exceed the dynamic curve
    peaks, _ = signal.find_peaks(final_intensity, height=threshold_curve)
    num_peaks = len(peaks)
    peak_score = calculate_microbe_ms_style_peak_score(num_peaks)
    
    # 5. Resolving Power Test
    # The description says resolving power uses peaks found ABOVE the intensity threshold
    res_power = calculate_resolving_power(mz_c, final_intensity, peaks)
    res_score = np.interp(res_power, [200, 1200], [0, 100])

    # Final Weighted Calculation
    total_score = (peak_score * weights['peaks']) + \
                  (noise_score * weights['noise']) + \
                  (baseline_score * weights['baseline']) + \
                  (res_score * weights['res'])

    if np.isnan(total_score):
        total_score = 0

    # Traffic Light Status (Thresholds 30 and 45)
    status = "GREEN" if total_score > 45 else ("RED" if total_score < 30 else "YELLOW")
    

    return {
        "Total QC Score": round(total_score, 1),
        "Status": status,
        "Sub-Scores": {
            "Peaks": round(peak_score, 1),
            "Noise": round(noise_score, 1),
            "Baseline": round(baseline_score, 1),
            "Resolving Power": round(res_score, 1)
        }
    }

def main():
    parser = argparse.ArgumentParser(description="QC for protein spectra using MicrobeMS-style metrics")
    parser.add_argument('--input_spectra', help='Path to input spectra file (e.g., mzML)')
    parser.add_argument('--output_path', help='Path to save QC .tsv report')
    args = parser.parse_args()

    input_file = Path(args.input_spectra)
    output_file = Path(args.output_path)

    if not input_file.exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist.")
    if not output_file.parent.exists():
        output_file.parent.mkdir(parents=True, exist_ok=True)

    with mzml.MzML(str(input_file)) as reader:
        with open(output_file, 'w', encoding='utf-8') as output_csv:
            headers = ['original_filename', 'scan', 'Total QC Score', 'Status', 'Peaks Score', 'Noise Score', 'Baseline Score', 'Resolving Power Score']
            output_writer = csv.DictWriter(output_csv, fieldnames=headers)
            output_writer.writeheader()
            for scan in reader:
                mz = scan['m/z array']
                intensity = scan['intensity array']
                try:
                    qc_results = microbe_ms_style_qc(mz, intensity)
                except Exception as e:
                    logging.error(f"Error processing scan {scan['id']} in file {input_file.name}: {e}")
                    qc_results = {
                        'Total QC Score': 'Error',
                        'Status': 'Error',
                        'Sub-Scores': {
                            'Peaks': 'Error',
                            'Noise': 'Error',
                            'Baseline': 'Error',
                            'Resolving Power': 'Error'
                        }
                    }
                output_writer.writerow({
                    'original_filename': input_file.name,
                    'scan': find_integer_at_end(scan['id']),
                    'Total QC Score': qc_results['Total QC Score'],
                    'Status': qc_results['Status'],
                    'Peaks Score': qc_results['Sub-Scores']['Peaks'],
                    'Noise Score': qc_results['Sub-Scores']['Noise'],
                    'Baseline Score': qc_results['Sub-Scores']['Baseline'],
                    'Resolving Power Score': qc_results['Sub-Scores']['Resolving Power']
                })

if __name__ == "__main__":
    main()