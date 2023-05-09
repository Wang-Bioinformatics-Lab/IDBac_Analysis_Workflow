import sys
import os
import argparse
import pandas as pd
import requests
import json
from pyteomics import mzxml, mzml
from massql import msql_fileloading

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
    parser.add_argument('query_file')
    parser.add_argument('database_mzML')
    parser.add_argument('database_json')
    parser.add_argument('results_folder')
    parser.add_argument('--merge_replicates', default="No")

    args = parser.parse_args()

    # Reading Query
    ms1_df, _ = load_data(args.query_file)

    bin_size = 10.0
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
    spectra_binned_df["filename"] = os.path.basename(args.query_file)

    bins_to_remove = []
    # merging replicates
    if args.merge_replicates == "Yes":
        # Lets do the merge
        all_bins = [x for x in spectra_binned_df.columns if x.startswith("BIN_")]
        for bin in all_bins:
            all_values = spectra_binned_df[bin]
            #print(bin, all_values)

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

    # Reading Database
    ms1_df_db, _ = load_data(args.database_mzML)

    # Filtering m/z
    ms1_df_db = ms1_df_db[ms1_df_db['mz'] < max_mz]

    # Bin the MS1 Data by m/z within each spectrum
    ms1_df_db['bin'] = (ms1_df_db['mz'] / bin_size).astype(int)

    # Now we need to group by scan and bin
    ms1_df_db = ms1_df_db.groupby(['scan', 'bin']).agg({'i': 'sum'}).reset_index()
    ms1_df_db["mz"] = ms1_df_db["bin"] * bin_size
    ms1_df_db["bin_name"] = "BIN_" + ms1_df_db["bin"].astype(str)
    
    # Turning each scan into a 1d vector that is the intensity value for each bin
    spectra_binned_db_df = ms1_df_db.pivot(index='scan', columns='bin_name', values='i').reset_index()
    spectra_binned_db_df["filename"] = os.path.basename(args.database_mzML)

    query_numerical_columns = [x for x in spectra_binned_df.columns if x.startswith("BIN_")]
    db_numerical_columns = [x for x in spectra_binned_db_df.columns if x.startswith("BIN_")]

    merged_numerical_columns = list(set(query_numerical_columns + db_numerical_columns))

    # Fill in the missing values with 0
    numerical_columns = [x for x in spectra_binned_df.columns if x.startswith("BIN_")]
    spectra_binned_df[numerical_columns] = spectra_binned_df[numerical_columns].fillna(0)
    
    # Add missing columns
    for column in merged_numerical_columns:
        if column not in spectra_binned_df.columns:
            spectra_binned_df[column] = 0

    query_data_np = spectra_binned_df[merged_numerical_columns].to_numpy()

    # Fill in the missing values with 0
    numerical_columns = [x for x in spectra_binned_db_df.columns if x.startswith("BIN_")]
    spectra_binned_db_df[numerical_columns] = spectra_binned_db_df[numerical_columns].fillna(0)

    # Add missing columns
    for column in merged_numerical_columns:
        if column not in spectra_binned_db_df.columns:
            spectra_binned_db_df[column] = 0

    database_data_np = spectra_binned_db_df[merged_numerical_columns].to_numpy()

    # Now lets do pairwise cosine similarity
    from sklearn.metrics.pairwise import cosine_similarity
    similarity_matrix = cosine_similarity(query_data_np, database_data_np)

    print(similarity_matrix)
    for i, row in enumerate(similarity_matrix):
        for j, item in enumerate(row):
            query_index = i
            database_index = j
            similarity = item

    # TODO: We will need to map the index back to the query scan and filename
    # TODO: We will need to map the database index back to the original

            





    

if __name__ == '__main__':
    main()