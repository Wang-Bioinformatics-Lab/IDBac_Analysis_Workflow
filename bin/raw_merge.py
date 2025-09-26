import argparse
import os

import numpy as np
from pyteomics import mzml
from psims.mzml.writer import MzMLWriter
from scipy.interpolate import interp1d

def merge_and_interpolate_mzml(input_file, output_file, interpolation_interval=0.1):
    # Read all spectra and collect m/z ranges
    mz_min, mz_max = float('inf'), float('-inf')
    spectra = []
    with mzml.MzML(input_file) as reader:
        for spec in reader:
            mzs = spec.get('m/z array')
            ints = spec.get('intensity array')
            if mzs is not None and ints is not None and len(mzs) > 0:
                mz_min = min(mz_min, mzs[0])
                mz_max = max(mz_max, mzs[-1])
                spectra.append((mzs, ints))

    # Define common m/z grid
    mz_grid = np.arange(mz_min, mz_max, interpolation_interval)
    merged_intensity = np.zeros_like(mz_grid)

    # Interpolate and sum intensities
    for mzs, ints in spectra:
        interp = interp1d(mzs, ints, kind='linear', bounds_error=False, fill_value=0)
        merged_intensity += interp(mz_grid)

    # Write to mzML
    with MzMLWriter(output_file) as writer:
        writer.controlled_vocabularies()
        with writer.run():
            with writer.spectrum_list(count=1):
                writer.write_spectrum(
                    mz_grid,
                    merged_intensity,
                    id=f"scan=1", params=[
                        "MS1 Spectrum",
                        {"ms level": 1},
                        {"total ion current": sum(merged_intensity)}
                    ])

def main():
    parser = argparse.ArgumentParser(description='Merge and interpolate mzML spectra')
    parser.add_argument('--input_file', required=True, help='Input mzML file')
    parser.add_argument('--output_folder', required=True, help='Output folder for merged mzML')
    args = parser.parse_args()

    output_file = os.path.join(args.output_folder, os.path.basename(args.input_file))
    merge_and_interpolate_mzml(args.input_file, output_file)

if __name__ == "__main__":
    main()