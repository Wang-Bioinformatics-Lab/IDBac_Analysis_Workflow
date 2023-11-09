import sys
import os
import argparse
import pandas as pd
import requests
import json
from pyteomics import mzxml, mzml
from massql import msql_fileloading
from tqdm import tqdm
import glob
from functools import lru_cache
import numpy as np

bin_size = 10.0

# Create an LRU Cache from functools
@lru_cache(maxsize=1000)
def _retreive_kb_metadata(database_id):
    url = "https://idbac-kb.gnps2.org/api/spectrum"
    params = {}
    params["database_id"] = database_id

    r = requests.get(url, params=params)

    spectrum_dict = r.json()

    return spectrum_dict

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

def load_database(database_filtered_json):
    all_database_spectra = json.load(open(database_filtered_json))

    formatted_database_spectra = []
    for spectrum in all_database_spectra:
        formatted_spectrum = {}
        formatted_spectrum["database_id"] = spectrum["database_id"]

        for peak in spectrum["peaks"]:
            formatted_spectrum["BIN_" + str(int(peak["mz"] / bin_size))] = peak["i"]

        formatted_database_spectra.append(formatted_spectrum)

    return pd.DataFrame(formatted_database_spectra)
    

    # Now we can make this into a dataframe that we can use as we did before
    
    # Now we need to create consensus spectra in the database
    # max_mz = 15000.0

    # # Filtering m/z
    # db_spectra = db_spectra[db_spectra['mz'] < max_mz]

    # # Bin the MS1 Data by m/z within each spectrum
    # db_spectra['bin'] = (db_spectra['mz'] / bin_size).astype(int)

    # # Now we need to group by scan and bin
    # db_spectra = db_spectra.groupby(['scan', 'bin']).agg({'i': 'sum'}).reset_index()
    # db_spectra["mz"] = db_spectra["bin"] * bin_size
    # db_spectra["bin_name"] = "BIN_" + db_spectra["bin"].astype(str)

    # # Turning each scan into a 1d vector that is the intensity value for each bin
    # spectra_binned_df = db_spectra.pivot(index='scan', columns='bin_name', values='i').reset_index()

    #return merged_spectra_df


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input_folder')
    parser.add_argument('database_filtered_json')
    parser.add_argument('output_results_tsv')
    parser.add_argument('--merge_replicates', default="Yes")
    parser.add_argument('--score_threshold', default=0.7, type=float)
    
    args = parser.parse_args()

    # Loading the database, this will also merge the spectra
    database_df = load_database(args.database_filtered_json)

    # Create a row count column, starting at 0, counting all the way up, will be useful for keeping track of things when we do matrix multiplication
    database_df["row_count"] = np.arange(len(database_df))

    # Prepping database for matching
    db_numerical_columns = [x for x in database_df.columns if x.startswith("BIN_")]

    all_input_files = glob.glob(os.path.join(args.input_folder, "*.mzML"))

    output_results_list = []
    #for input_filename in all_input_files[:5]:
    for input_filename in all_input_files:
        # Reading Query
        ms1_df, _ = load_data(input_filename)

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
        
        # Formatting Database to fill NAs
        
        # Turning each scan into a 1d vector that is the intensity value for each bin
        spectra_binned_db_df = database_df

        query_numerical_columns = [x for x in spectra_binned_df.columns if x.startswith("BIN_")]

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

        #print(similarity_matrix)
        for i, row in enumerate(similarity_matrix):
            for j, item in enumerate(row):
                query_index = i
                database_index = j
                similarity = item

                if similarity < args.score_threshold:
                    continue

                result_dict = {}
                result_dict["query_index"] = query_index
                result_dict["database_index"] = database_index
                result_dict["similarity"] = similarity
                result_dict["query_filename"] = os.path.basename(input_filename)

                output_results_list.append(result_dict)

    small_database_df = database_df[["row_count", "database_id"]]
    small_database_df["database_scan"] = small_database_df["database_id"]

    output_results_df = pd.DataFrame(output_results_list)
    if len(output_results_df) == 0:
        print("No matches found")
        open(args.output_results_tsv, "w").write("\n")

        exit(0)
    
    # Here we want to take the results from the search
    output_results_df = output_results_df.merge(small_database_df, left_on="database_index", right_on="row_count", how="left")
    output_results_df["database_id"] = output_results_df["database_scan"]
    output_results_df["database_scan"] = output_results_df["row_count"] + 1
                
    # Enrich the library information by hitting the web api
    all_database_metadata_json = requests.get("https://idbac-kb.gnps2.org/api/spectra").json()
    all_database_metadata_df = pd.DataFrame(all_database_metadata_json)
    all_database_metadata_df = all_database_metadata_df[["database_id", "Strain name", "Culture Collection", "Sample name", "Genbank accession"]]
    
    # lets merge the results with the metadata
    output_results_df = output_results_df.merge(all_database_metadata_df, left_on="database_id", right_on="database_id", how="left")

    # rename columns
    output_results_df = output_results_df.rename(columns={"Strain name": "db_strain_name", "Culture Collection": "db_culture_collection", "Sample name": "db_sample_name", "Genbank accession": "db_genbank_accession"})

    # Output data
    output_results_df.to_csv(args.output_results_tsv, sep="\t", index=False)

    



if __name__ == '__main__':
    main()