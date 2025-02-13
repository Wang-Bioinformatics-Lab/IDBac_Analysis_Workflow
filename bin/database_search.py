import sys
import os
import argparse
import pandas as pd
import requests
import json
from tqdm import tqdm
import glob
from functools import lru_cache
import numpy as np
import logging
from sklearn.metrics.pairwise import cosine_distances, euclidean_distances

from utils import load_data, spectrum_binner, compute_distances_binned, write_spectra_df_to_mzML

# Create an LRU Cache from functools
@lru_cache(maxsize=1000)
def _retreive_kb_metadata(database_id):
    url = "https://idbac.org/api/spectrum"
    params = {}
    params["database_id"] = database_id

    r = requests.get(url, params=params)

    spectrum_dict = r.json()

    return spectrum_dict

def load_database(database_filtered_json, mz_min, mz_max, bin_size):
    all_database_spectra = json.load(open(database_filtered_json))

    formatted_database_spectra = []
    for spectrum in all_database_spectra:
        formatted_spectrum = {}
        formatted_spectrum["database_id"] = spectrum["database_id"]
        formatted_spectrum["genus"] = spectrum["genus"]
        formatted_spectrum["species"] = spectrum["species"]

        for peak in spectrum["peaks"]:
            if peak["mz"] < mz_min or peak["mz"] > mz_max:
                continue
            formatted_spectrum["BIN_" + str(int(peak["mz"] / bin_size))] = peak["i"]

        formatted_database_spectra.append(formatted_spectrum)

    return pd.DataFrame(formatted_database_spectra)

def compute_db_db_distance(database_df, db_numerical_columns, output_path, distance_metric="cosine"):

    """ This function computes the pairwise distance between the the database results in the first 
    argument, and themselves. This is used to fill out the remainder of the distance matrix.
    
    Args:
    database_df: pd.DataFrame
        The database results dataframe
    db_numerical_columns: list
        The list of numerical columns in the database dataframe
    output_path: str
        The path to save the results
    distance_metric: str (default: "cosine")
        The distance metric to use for the search
        
    Returns:
    output_results_df: pd.DataFrame
        The results of the database distance search with the following columns:
        ["left_index", "right_index", "database_id_left", "database_id_right", "database_scan_left", "database_scan_right", "distance"]
    """
    database_data_np = database_df[db_numerical_columns].to_numpy()
    db_db_distance = compute_distances_binned(database_data_np, distance_metric=distance_metric)
    
    output_results_list = []
    
    for i, row in enumerate(db_db_distance):
        for j in range (i+1, len(row)):         # Skip the diagonal, to save some time
            distance = db_db_distance[i][j]
            left_index = i
            right_index = j
            
            # Upper diagonal
            result_dict = {}
            result_dict["left_index"] = left_index
            result_dict["right_index"] = right_index
            result_dict["database_id_left"] = database_df.iloc[left_index]["database_id"]
            result_dict["database_id_right"] = database_df.iloc[right_index]["database_id"]
            result_dict["database_scan_left"] = database_df.iloc[left_index]["database_scan"]
            result_dict["database_scan_right"] = database_df.iloc[right_index]["database_scan"]
            result_dict["distance"] = distance

            output_results_list.append(result_dict)
            
            # Lower diagonal
            result_dict = {}
            result_dict["left_index"] = right_index
            result_dict["right_index"] = left_index
            result_dict["database_id_left"] = database_df.iloc[right_index]["database_id"]
            result_dict["database_id_right"] = database_df.iloc[left_index]["database_id"]
            result_dict["database_scan_left"] = database_df.iloc[right_index]["database_scan"]
            result_dict["database_scan_right"] = database_df.iloc[left_index]["database_scan"]
            result_dict["distance"] = distance
            
            output_results_list.append(result_dict)
            
    if len(output_results_list) == 0:
        # Save an empty dataframe with a header
        output_results_df = pd.DataFrame(columns=["left_index", "right_index", "database_id_left", "database_id_right", "database_scan_left", "database_scan_right", "distance"])
        output_results_df.to_csv(output_path, sep="\t", index=False)
        return output_results_df
                    
    output_results_df = pd.DataFrame(output_results_list)
    output_results_df.to_csv(output_path, sep="\t", index=False)
    
    return output_results_df

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input_folder')
    parser.add_argument('database_filtered_json')
    parser.add_argument('output_results_tsv')
    parser.add_argument('complete_output_results_tsv')
    parser.add_argument('output_db_db_distance_tsv')
    parser.add_argument('--merge_replicates', default="Yes")
    parser.add_argument('--score_threshold', default=0.7, type=float)
    parser.add_argument('--distance', default="cosine", help="The distance metric to use for the search", choices=["cosine", "euclidean", "presence"])
    parser.add_argument('--bin_size', type=float, help="Size of the spectra bins for distance calculations.", required=True)
    parser.add_argument('--mass_range_lower', default=2000.0, type=float, help="Minimum m/z value to consider for binning.")
    parser.add_argument('--mass_range_upper', default=20000.0, type=float, help="Maximum m/z value to consider for binning.")
    # I actually don't know how the checkbox form input works, so this is a placeholder for now TODO, Potential BUG
    parser.add_argument('--seed_genera', default='', help="Comma separated list of genera to seed the search with.")
    parser.add_argument('--seed_species', default='', help="Comma separated list of species to seed the search with.")
    parser.add_argument('--debug', action='store_true')
    
    args = parser.parse_args()
    
    bin_size = float(args.bin_size)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
        # Log all params
        for arg in vars(args):
            logging.debug("%s: %s", arg, getattr(args, arg))
    else:
        logging.basicConfig(level=logging.INFO)

    # Loading the database, this will also merge the spectra
    database_df = load_database(args.database_filtered_json, args.mass_range_lower, args.mass_range_upper, bin_size)
    logging.debug("Database shape: {}".format(database_df.shape))
    logging.debug("database_df {}".format(database_df))

    # Create a row count column, starting at 0, counting all the way up, will be useful for keeping track of things when we do matrix multiplication
    database_df["row_count"] = np.arange(len(database_df))

    # Prepping database for matching
    db_numerical_columns = [x for x in database_df.columns if x.startswith("BIN_")]
    db_metadata = [x for x in database_df.columns if x not in db_numerical_columns]
    for required_col in ['genus', 'species', 'database_id']:
        if required_col not in database_df.columns:
            raise ValueError("Database metadata does not contain required column: {}".format(required_col))
    db_metadata = database_df.loc[:, db_metadata]

    # Get the index of selected genera and species seeds
    db_metadata.reset_index(inplace=True)   # Reset index to ensure contiguous numbering
    if args.seed_genera:
        seed_genera = args.seed_genera.split(";")
        seed_genera_indices = db_metadata[db_metadata["genus"].isin(seed_genera)].index
    else:
        seed_genera_indices = []
    if args.seed_species:
        seed_species = args.seed_species.split(";")
        seed_species_indices = db_metadata[db_metadata["species"].isin(seed_species)].index
    else:
        seed_species_indices = []

    all_input_files = glob.glob(os.path.join(args.input_folder, "*.mzML"))

    output_results_list = []

    for input_filename in all_input_files:
        # Reading Query
        ms1_df, _ = load_data(input_filename)
        
        spectra_binned_df = spectrum_binner(ms1_df,
                                            input_filename,
                                            bin_size=bin_size,
                                            min_mz=args.mass_range_lower,
                                            max_mz=args.mass_range_upper,
                                            merge_replicates="Yes")

        output_filename = os.path.join('query_spectra',  os.path.basename(input_filename))
        write_spectra_df_to_mzML(output_filename, spectra_binned_df, bin_size)
        
        # Formatting Database to fill NAs
        
        # Turning each scan into a 1d vector that is the intensity value for each bin
        spectra_binned_db_df = database_df

        query_numerical_columns = [x for x in spectra_binned_df.columns if x.startswith("BIN_")]
        merged_numerical_columns = list(set(query_numerical_columns + db_numerical_columns))

        # Fill in the missing values with 0
        numerical_columns = [x for x in spectra_binned_df.columns if x.startswith("BIN_")]
        spectra_binned_df[numerical_columns] = spectra_binned_df[numerical_columns].fillna(0)
        
        # Add missing columns
        missing_columns = [column for column in merged_numerical_columns if column not in spectra_binned_df.columns]
        missing_df = pd.DataFrame(0, index=spectra_binned_df.index, columns=missing_columns)
        spectra_binned_df = pd.concat([spectra_binned_df, missing_df], axis=1)


        query_data_np = spectra_binned_df[merged_numerical_columns].to_numpy()
        logging.debug("query_data_np shape: {}".format(query_data_np.shape))
        logging.debug("query_data_np: {}".format(query_data_np))
        logging.debug("query_data_np max: {}".format(np.max(query_data_np)))
        logging.debug("query_data_np min: {}".format(np.min(query_data_np)))
        logging.debug("query_data_np mean: {}".format(np.mean(query_data_np)))

        # Fill in the missing values with 0
        numerical_columns = [x for x in spectra_binned_db_df.columns if x.startswith("BIN_")]
        spectra_binned_db_df[numerical_columns] = spectra_binned_db_df[numerical_columns].fillna(0)

        # Add missing columns
        for column in merged_numerical_columns:
            if column not in spectra_binned_db_df.columns:
                spectra_binned_db_df[column] = 0

        database_data_np = spectra_binned_db_df[merged_numerical_columns].to_numpy()

        # Now lets do pairwise cosine distance
        distance_matrix = compute_distances_binned(query_data_np, database_data_np, distance_metric=args.distance)

        for i, row in enumerate(distance_matrix):
            for j, item in enumerate(row):
                query_index = i
                database_index = j
                distance = item
               
                if len(seed_genera_indices) == 0 and len(seed_species_indices) == 0:
                    if distance < args.score_threshold:
                        result_dict = {}
                        result_dict["query_index"] = query_index
                        result_dict["database_index"] = database_index
                        result_dict["distance"] = distance
                        result_dict["query_filename"] = os.path.basename(input_filename)

                        output_results_list.append(result_dict)
                else:
                    is_seed = False
                    if database_index in seed_genera_indices or \
                        database_index in seed_species_indices:
                        is_seed = True
                    if is_seed:
                        result_dict = {}
                        result_dict["query_index"] = query_index
                        result_dict["database_index"] = database_index
                        result_dict["distance"] = distance
                        result_dict["query_filename"] = os.path.basename(input_filename)

                        output_results_list.append(result_dict)
                    
        # We will also need a complete distance_matrix for the query-db pairs (in contrast to the above on that's trimmed by the score_threshold)
        complete_output_results_list = []
        # Get columns where there is at least one match less than threshold
        thresholded_indices = np.where(np.any(distance_matrix < args.score_threshold, axis=0))[0]
        for i, row in enumerate(distance_matrix):
            for j in thresholded_indices:
                distance = row[j]
                query_index = i
                database_index = j
                
                result_dict = {}
                result_dict["query_index"] = query_index
                result_dict["database_index"] = database_index
                result_dict["distance"] = distance
                result_dict["query_filename"] = os.path.basename(input_filename)
                
                complete_output_results_list.append(result_dict)
                

    small_database_df = database_df[["row_count", "database_id"]]
    small_database_df["database_scan"] = small_database_df["database_id"]

    output_results_df = pd.DataFrame(output_results_list, columns=["query_index", "database_index", "distance", "query_filename"])
    complete_output_results_df = pd.DataFrame(complete_output_results_list, columns=["query_index", "database_index", "distance", "query_filename"])
    if len(output_results_df) == 0:
        print("No matches found")
        # Write headers only
        with open(args.output_results_tsv, "w") as f:
            f.write("query_index\tdatabase_index\tdistance\tquery_filename\n")
        with open(args.complete_output_results_tsv, "w") as f:
            f.write("query_index\tdatabase_index\tdistance\tquery_filename\n")
        with open(args.output_db_db_distance_tsv, "w") as f:
            f.write("left_index\tright_index\tdatabase_id_left\tdatabase_id_right\tdatabase_scan_left\tdatabase_scan_right\tdistance\n")
        exit(0)
    
    # Here we want to take the results from the search
    output_results_df = output_results_df.merge(small_database_df, left_on="database_index", right_on="row_count", how="left")
    output_results_df["database_id"] = output_results_df["database_scan"]
    output_results_df["database_scan"] = output_results_df["row_count"] + 1
    complete_output_results_df = complete_output_results_df.merge(small_database_df, left_on="database_index", right_on="row_count", how="left")
    
    # Sanity Check: We assume these to be unique below
    if len(database_df.database_id.unique()) != len(database_df.database_id):
        raise ValueError("Database ID is not unique")
    
    # Perform pairwise distance of database result hits to fill out the remainder of the distance matrix.
    database_results_df = output_results_df[["database_id", "database_scan"]].drop_duplicates()

    # Add numerical columns to the database results, for the distance calculation
    database_results_df = database_results_df.merge(database_df, left_on="database_id", right_on="database_id", how="left")
    
    compute_db_db_distance(database_results_df, db_numerical_columns, args.output_db_db_distance_tsv, distance_metric=args.distance)
                
    # Enrich the library information by hitting the web api
    all_database_metadata_json = requests.get("https://idbac.org/api/spectra").json()
    all_database_metadata_df = pd.DataFrame(all_database_metadata_json)
    all_database_metadata_df = all_database_metadata_df[["database_id", "Strain name", "Culture Collection", "Sample name", "Genbank accession"]]
    
    # lets merge the results with the metadata
    output_results_df = output_results_df.merge(all_database_metadata_df, left_on="database_id", right_on="database_id", how="left")
    complete_output_results_df = complete_output_results_df.merge(all_database_metadata_df, left_on="database_id", right_on="database_id", how="left")

    # rename columns
    output_results_df = output_results_df.rename(columns={"Strain name": "db_strain_name", "Culture Collection": "db_culture_collection", "Sample name": "db_sample_name", "Genbank accession": "db_genbank_accession"})

    # Output data
    output_results_df.to_csv(args.output_results_tsv, sep="\t", index=False)
    complete_output_results_df.to_csv(args.complete_output_results_tsv, sep="\t", index=False)


if __name__ == '__main__':
    main()