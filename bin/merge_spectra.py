import sys
import os
import argparse
import pandas as pd
import uuid
import json
from massql import msql_fileloading
from pyteomics import mzxml, mzml
from psims.mzml.writer import MzMLWriter
from tqdm import tqdm
import glob


def load_data(input_filename):
    try:
        ms1_df, ms2_df = msql_fileloading.load_data(input_filename)

        return ms1_df, ms2_df
    except:
        print("Error loading data, falling back on default")

    MS_precisions = {
        1: 5e-6,
        2: 20e-6,
        3: 20e-6,
        4: 20e-6,
        5: 20e-6,
        6: 20e-6,
        7: 20e-6
    }
    
    ms1_df = pd.DataFrame()
    ms2_df = pd.DataFrame()

    all_mz = []
    all_i = []
    all_scan = []
    
    # TODO: read the mzML directly
    with mzml.read(input_filename) as reader:
        for spectrum in tqdm(reader):
            try:
                scan = spectrum["id"].replace("scanId=", "").split("scan=")[-1]
            except:
                scan = spectrum["id"]

            mz = spectrum["m/z array"]
            intensity = spectrum["intensity array"]

            all_mz += list(mz)
            all_i += list(intensity)
            all_scan += len(mz) * [scan]

            print(spectrum["id"])
            
    if len(all_mz) > 0:
        ms1_df['i'] = all_i
        ms1_df['mz'] = all_mz
        ms1_df['scan'] = all_scan

    return ms1_df, ms2_df


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input_folder')
    parser.add_argument('output_folder')
    parser.add_argument('--merge_replicates', default="No")
    parser.add_argument('--bin_size', default=10.0, type=float)

    args = parser.parse_args()

    bin_size = args.bin_size

    # Lets read all the spectra that are coming out of the input_folder
    all_input_files = glob.glob(os.path.join(args.input_folder, "*.mzML"))
    all_input_files.sort()

    all_spectra_df_list = []

    for input_filename in all_input_files:
        print("Loading data from {}".format(input_filename))
        ms1_df, ms2_df = load_data(input_filename)

        max_mz = 15000.0

        # Filtering m/z
        ms1_df = ms1_df[ms1_df['mz'] < max_mz]

        # Bin the MS1 Data by m/z within each spectrum
        ms1_df['bin'] = (ms1_df['mz'] / bin_size).astype(int)

        # Now we need to group by scan and bin
        ms1_df = ms1_df.groupby(['scan', 'bin']).agg({'i': 'sum'}).reset_index()
        ms1_df["mz"] = ms1_df["bin"] * bin_size
        ms1_df["bin_name"] = "BIN_" + ms1_df["bin"].astype(str)
        
        # Turning each scan into a 1d vector that is the intensity value for each bin
        spectra_binned_df = ms1_df.pivot(index='scan', columns='bin_name', values='i').reset_index()
        spectra_binned_df["filename"] = os.path.basename(input_filename)

        bins_to_remove = []
        # merging replicates
        if args.merge_replicates == "Yes":
            # Lets do the merge
            all_bins = [x for x in spectra_binned_df.columns if x.startswith("BIN_")]
            for bin in all_bins:
                all_values = spectra_binned_df[bin]

                # Count non-zero values
                non_zero_count = len(all_values[all_values > 0])

                # Calculate percent non-zero
                percent_non_zero = non_zero_count / len(all_values)

                if percent_non_zero < 0.5:
                    bins_to_remove.append(bin)

            # Removing the bins
            spectra_binned_df = spectra_binned_df.drop(bins_to_remove, axis=1)

            # Now lets get the mean for each bin
            spectra_binned_df = spectra_binned_df.groupby("filename").mean().reset_index()
            spectra_binned_df["scan"] = "merged"

        # Writing an mzML file with the merged spectra
        output_filename = os.path.join(args.output_folder, os.path.basename(input_filename))
        with MzMLWriter(open(output_filename, 'wb'), close=True) as out:
            # Add default controlled vocabularies
            out.controlled_vocabularies()
            # Open the run and spectrum list sections
            with out.run(id="my_analysis"):
                spectrum_count = len(spectra_binned_df)

                spectrum_list = spectra_binned_df.to_dict(orient="records")

                with out.spectrum_list(count=spectrum_count):
                    scan = 1
                    for spectrum_dict in spectrum_list:

                        all_keys = list(spectrum_dict.keys())

                        mz_array = [float(key.replace("BIN_", "")) * bin_size for key in all_keys if key.startswith("BIN_") if spectrum_dict[key] > 0]
                        intensity_array = [spectrum_dict[key] for key in all_keys if key.startswith("BIN_") if spectrum_dict[key] > 0]

                        # Write scan
                        out.write_spectrum(
                            mz_array, intensity_array,
                            id="scan={}".format(scan), params=[
                                "MS1 Spectrum",
                                {"ms level": 1},
                                {"total ion current": sum(intensity_array)}
                            ])
                        
                        scan += 1


if __name__ == '__main__':
    main()