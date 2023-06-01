import sys
import os
import argparse
import pandas as pd
import uuid
import json
from massql import msql_fileloading
from pyteomics import mzxml, mzml
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
    parser.add_argument('output_summary')

    args = parser.parse_args()

    # Lets read all the spectra that are coming out of the input_folder
    all_input_files = glob.glob(os.path.join(args.input_folder, "*.mzML"))

    output_scan_list = []
    for input_file in all_input_files:
        ms1_df, _ = load_data(input_file)

        all_scans = ms1_df["scan"].unique()

        for scan in all_scans:
            scan_dict = {}
            scan_dict["scan"] = scan
            scan_dict["filename"] = os.path.basename(input_file)
        
            output_scan_list.append(scan_dict)

    output_df = pd.DataFrame(output_scan_list)
    output_df.to_csv(args.output_summary, sep="\t", index=False)


if __name__ == '__main__':
    main()