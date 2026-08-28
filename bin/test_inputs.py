import argparse
import os
import sys
import numpy as np
from pyteomics import mzml
import lxml
import csv
import re
import logging


USER_ERROR_PREFIX = "IDBAC_USER_ERROR:"


def fail_user(message):
    """Emit a machine-detectable user error and stop the validation process."""
    print(f"{USER_ERROR_PREFIX} {message}", file=sys.stderr, flush=True)
    sys.exit(1)


def find_integer_at_end(string):
    return int(re.search(r'\d+$', string).group()) if re.search(r'\d+$', string) else 'N/A'

def breaking_errors(input_file:str):
    """For checking all errors so critical, NextFlow should terminate immediately."""
    # Ensure file is .mzML
    if not input_file.lower().endswith('.mzml'):
        fail_user(f"Input file '{os.path.basename(input_file)}' is not an mzML file. Please provide a valid mzML file.")
    
    # Check that this is not an ESI file (we won't check that it's MALDI because I don't trust MSConvert Enough)
    # Check for minimal criteria for <cvParam accession="MS:1000073"...>
    root = lxml.etree.parse(input_file).getroot()
    esi_params = root.findall('.//{*}cvParam[@accession="MS:1000073"]')
    if len(esi_params) > 0:
        fail_user(f"Input file '{os.path.basename(input_file)}' appears to be an ESI file. Please provide a valid MALDI mzML file.")


def validate_protein_mass_range(input_file, mass_range_lower, mass_range_upper, spectra=None):
    """Require at least one usable protein peak inside the configured search range.

    This catches uncalibrated or incorrectly converted axes (for example, a
    time-of-flight axis labelled as m/z) before downstream processing turns the
    file into an empty binned spectrum.
    """
    if mass_range_lower >= mass_range_upper:
        fail_user(
            f"The protein mass range is invalid: lower bound {mass_range_lower:g} "
            f"must be less than upper bound {mass_range_upper:g}."
        )

    observed_min = np.inf
    observed_max = -np.inf
    usable_peak_count = 0

    reader = spectra if spectra is not None else mzml.MzML(input_file)
    try:
        for scan in reader:
            mz_array = np.asarray(scan.get('m/z array', []), dtype=float)
            intensity_array = np.asarray(scan.get('intensity array', []), dtype=float)

            finite_mz = mz_array[np.isfinite(mz_array)]
            if finite_mz.size:
                observed_min = min(observed_min, float(finite_mz.min()))
                observed_max = max(observed_max, float(finite_mz.max()))

            if mz_array.size != intensity_array.size:
                continue

            usable = (
                np.isfinite(mz_array)
                & np.isfinite(intensity_array)
                & (intensity_array > 0)
                & (mz_array > mass_range_lower)
                & (mz_array < mass_range_upper)
            )
            usable_peak_count += int(np.count_nonzero(usable))
    finally:
        if spectra is None:
            reader.close()

    if usable_peak_count == 0:
        if np.isfinite(observed_min) and np.isfinite(observed_max):
            observed = f"{observed_min:g} to {observed_max:g}"
        else:
            observed = "no finite m/z values"
        fail_user(
            f"Protein spectrum '{os.path.basename(input_file)}' has no finite, positive-intensity peaks "
            f"inside the configured mass range {mass_range_lower:g} to {mass_range_upper:g} m/z. "
            f"Observed m/z range: {observed}. Regenerate the mzML with a calibrated m/z axis or "
            "correct the configured protein mass range."
        )


def validate_file(input_file:str, output_file:str)->int:
    """This function will validate the input file for common errors and output them to a csv file if errors are found."""
    
    # Check that the file is not empty
    if os.path.getsize(input_file) == 0:
        with open(output_file, 'w', encoding='utf-8') as output_csv:
            headers = ['error_level', 'scan', 'original_filename', 'error']
            output_writer = csv.DictWriter(output_csv, fieldnames=headers)
            output_writer.writerow({'error_level': 'critical', 
                                    'scan': 'N/A',
                                    'original_filename': input_file, 
                                    'error': 'File is empty. If the file was uploaded, please check the uploaded file to ensure it was successful.'})
        return 1
    
    output_status = 0
    with mzml.MzML(input_file) as reader:
        with open(output_file, 'w', encoding='utf-8') as output_csv:
            headers = ['error_level', 'scan', 'original_filename', 'error']
            output_writer = csv.DictWriter(output_csv, fieldnames=headers)           
            
            # Check to make sure scans are readable
            if len(reader) == 0:
                output_writer.writerow({'original_filename': input_file, 'error': 'No scans found'})
                output_status = 1
            # Check to make sure scans have nonzero intensity arrays
            scan_status = np.zeros(len(reader), dtype=bool)
            for scan_idx, scan in enumerate(reader):
                if 'intensity array' not in scan:
                    output_writer.writerow({'error_level': 'warning',
                                            'scan': find_integer_at_end(scan['id']),
                                            'original_filename': input_file,
                                            'error': f"Scan {scan['id']} is missing intensity array"})
                    output_status = 1
                    scan_status[scan_idx] = True
                    logging.info(f"The following scan is missing intensity array: {scan['id']}")
                    logging.info(str(scan))

                if len(scan['intensity array']) == 0:
                    output_writer.writerow({'error_level': 'warning',
                                            'scan': find_integer_at_end(scan['id']),
                                            'original_filename': input_file,
                                            'error': f"Scan {scan['id']} has empty intensity array"})
                    output_status = 1
                    scan_status[scan_idx] = True
                    logging.info(f"The following scan reported empty intensity array: {scan['id']}")
                    logging.info(str(scan))

                if max(scan['intensity array']) == 0:
                    output_writer.writerow({'error_level': 'warning',
                                            'scan': find_integer_at_end(scan['id']),
                                            'original_filename': input_file,
                                            'error': f"Scan {scan['id']} contains no peaks"})
                    output_status = 1
                    scan_status[scan_idx] = True
                    logging.info(f"The following scan reported zero intensity: {scan['id']}")
                    logging.info(str(scan))

            if np.all(scan_status):
                output_writer.writerow({'error_level': 'critical',
                                        'scan': 'N/A',
                                        'original_filename': input_file,
                                        'error': 'All scans in this file are unreadable or empty. Please check the file and try again.'})
                output_status = 1
            
    return output_status

def main():
    parser = argparse.ArgumentParser(description='Check mzML files for common errors. If errors are found, output them to a csv and exit with 1.')
    parser.add_argument('--input_file', help='Path to the input file')
    parser.add_argument('--output_file', help='Path to the output file')
    parser.add_argument('--spectrum_type', choices=['protein', 'other'], default='other')
    parser.add_argument('--mass_range_lower', type=float)
    parser.add_argument('--mass_range_upper', type=float)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    for var, value in vars(args).items():
        logging.info(f'Argument {var}: {value}')
    
    breaking_errors(args.input_file)

    if args.spectrum_type == 'protein':
        if args.mass_range_lower is None or args.mass_range_upper is None:
            fail_user("Protein spectrum validation requires both mass-range bounds.")
        validate_protein_mass_range(args.input_file, args.mass_range_lower, args.mass_range_upper)

    status = validate_file(args.input_file, args.output_file)
    # sys.exit(status) # Nextflow doesn't have a provision to output files if the process fails. If this comes in the future, it will be a good way to warn users

if __name__ == "__main__":
    main()
