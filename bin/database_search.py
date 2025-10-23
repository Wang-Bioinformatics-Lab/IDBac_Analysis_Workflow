import os
import argparse
import pandas as pd
import requests
import json
import glob
from functools import lru_cache
import numpy as np
import logging
from sklearn.metrics.pairwise import cosine_distances, euclidean_distances
from collections import defaultdict
from utils import load_data, spectrum_binner, compute_distances_binned, write_spectra_df_to_mzML, peak_filtering
from pathlib import Path
import yaml
from time import time

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

        if "peaks" in spectrum and "ml_embedding" in spectrum:
            raise ValueError("Spectrum should not contain both 'peaks' and 'ml_embedding'. Please check the database file.")

        # If ml_embedding is present, use it instead of peaks (convert list values into BIN_ keys)
        if "ml_embedding" in spectrum:
            for i, value in enumerate(spectrum["ml_embedding"]):
                formatted_spectrum["BIN_" + str(i)] = value
        else:
            for peak in spectrum["peaks"]:
                if peak["mz"] < mz_min or peak["mz"] > mz_max:
                    continue
                formatted_spectrum["BIN_" + str(int(peak["mz"] / bin_size))] = peak["i"]



        formatted_database_spectra.append(formatted_spectrum)

    return pd.DataFrame(formatted_database_spectra)

def compute_db_db_distance(database_df, db_numerical_columns, output_path, distance_metric="cosine", unmatched_peak_penalty=1.0):

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
    start_time = time()
    logging.info("Computing database to database distance matrix.")
    db_db_distance = compute_distances_binned(database_data_np, distance_metric=distance_metric, penalty=unmatched_peak_penalty)
    logging.info("Computed database to database distance matrix in {:.2f} seconds".format(time() - start_time))

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

def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument(
        '--input_folder'
        )
    parser.add_argument(
        '--database_filtered_json'
        )
    parser.add_argument(
        '--output_search_results_tsv',
        help="Output file for database search results",
        )
    parser.add_argument(
        '--complete_output_results_tsv',
        help="Output file for database search results containing all query-db distances for anything returned in the database search"
        )
    parser.add_argument(
        '--output_db_db_distance_tsv',
        help="Output file for database to database distances",
        )
    parser.add_argument(
        '--output_query_query_distances_tsv',
        help="Output file for query to query distances",
        )
    parser.add_argument(
        '--merge_replicates',
        default="Yes"
        )
    parser.add_argument(
        '--score_threshold',
        default=0.7,
        type=float
        )
    parser.add_argument(
        '--distance',
        default="cosine",
        help="The distance metric to use for the search",
        choices=["cosine", "euclidean", "presence", "reverse_cosine", "reverse_presence"]
        )
    parser.add_argument(
        '--unmatched_peak_penalty',
        default=1.0,
        type=float,
        help="The penalty to apply for unmatched peaks when using reverse_cosine or reverse_presence distance metrics. Ranges from 0.0 to 1.0."
    )
    parser.add_argument(
        '--bin_size',
        type=float,
        help="Size of the spectra bins for distance calculations.",
        required=True
        )
    parser.add_argument(
        '--mass_range_lower',
        default=2000.0,
        type=float,
        help="Minimum m/z value to consider for binning."
        )
    parser.add_argument(
        '--mass_range_upper',
        default=20000.0,
        type=float,
        help="Maximum m/z value to consider for binning."
        )
    parser.add_argument(
        '--seed_genera',
        default='',
        help="Comma separated list of genera to seed the search with."
        )
    parser.add_argument(
        '--seed_species',
        default='',
        help="Comma separated list of species to seed the search with."
        )
    parser.add_argument(
        '--debug',
        action='store_true'
        )
    parser.add_argument(
        '--ml',
        action='store_true',
        help="Use the ML database for the search. If not set, the standard database will be used."
        )
    parser.add_argument(
        '--MALDI_instrument',
        default='general',
        help="The type of MALDI instrument used. This affects the peak filtering step."
        )
    parser.add_argument(
        '--config',
        default=None,
        help="Path to the config file containing the relative intensity thresholds for different MALDI instruments."
        )
    
    # Dump all args
    for arg in vars(parser.parse_args()):
        logging.info("%s: %s", arg, getattr(parser.parse_args(), arg))

    # Check unmatched_peak_penalty in [0.0, 1.0]
    if parser.parse_args().unmatched_peak_penalty < 0.0 or parser.parse_args().unmatched_peak_penalty > 1.0:
        raise ValueError("unmatched_peak_penalty must be between 0.0 and 1.0")

    return parser.parse_args()


def load_and_prepare_database(database_filtered_json, mass_range_lower, mass_range_upper, bin_size):
    """
    Load the database from a JSON file and prepare it for matching.
    
    Args:
    database_filtered_json: str
        Path to the JSON file containing the database spectra.
    mass_range_lower: float
        Minimum m/z value to consider for binning.
    mass_range_upper: float
        Maximum m/z value to consider for binning.
    bin_size: float
        Size of the spectra bins for distance calculations.
    
    Returns:
    database_df: pd.DataFrame
        DataFrame containing the prepared database spectra.
    """
    # Loading the database, this will also merge the spectra
    database_df = load_database(database_filtered_json, mass_range_lower, mass_range_upper, bin_size)
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
    db_metadata.reset_index(inplace=True)   # Reset index to ensure contiguous numbering

    return database_df, db_numerical_columns, db_metadata

def get_seed_indices(db_metadata, seed_genera=None, seed_species=None):
    """
    Get the indices of the database entries that match the seed genera and species.
    
    Args:
    db_metadata: pd.DataFrame
        DataFrame containing the database spectra.
    seed_genera: str
        List of genera to seed the search with.
    seed_species: str
        List of species to seed the search with.
    
    Returns:
    seed_genera_indices: pd.Index
        Indices of the database entries that match the seed genera.
    seed_species_indices: pd.Index
        Indices of the database entries that match the seed species.
    """
    seed_genera_indices = pd.Index([])
    seed_species_indices = pd.Index([])

    if seed_genera:
        seed_genera = seed_genera.split(";")
        seed_genera = [x.strip() for x in seed_genera]
        seed_genera_indices = db_metadata[db_metadata["genus"].isin(seed_genera)].index
        
    if seed_species:
        seed_species = seed_species.split(";")
        seed_species = [x.strip() for x in seed_species]
        seed_species_indices = db_metadata[db_metadata["species"].isin(seed_species)].index

    return seed_genera_indices, seed_species_indices

def spectrum_iterator(input_data, bin_size, args):
    is_mzmls = all([x.endswith(".mzML") for x in input_data])
    is_feather = len(input_data) == 1 and input_data[0].endswith(".feather")
    print("is_mzmls", is_mzmls, "is_feather", is_feather, 'input_data', input_data)
    if is_mzmls:
        for input_filename in input_data:
            ms1_df, _ = load_data(input_filename)
            spectra_binned_df = spectrum_binner(ms1_df,
                                                input_filename,
                                                bin_size=bin_size,
                                                min_mz=args.mass_range_lower,
                                                max_mz=args.mass_range_upper,
                                                merge_replicates="Yes")

            # Load the config
            relative_intensity=0.0
            if args.config is not None:
                logging.info("Loading config file from {}".format(args.config))
                with open(args.config, "r", encoding='utf-8') as f:
                    config = yaml.safe_load(f)
                    relative_intensity = config[args.MALDI_instrument].get("relative_intensity", 0.0)
                    print("Using relative_intensity", relative_intensity, "for MALDI instrument", args.MALDI_instrument)

                spectra_binned_df = peak_filtering(spectra_binned_df, relative_intensity)
            else:
                logging.info("No config file provided, skipping peak filtering step.")

            yield input_filename, spectra_binned_df
    elif is_feather:
        df = pd.read_feather(input_data[0])
        print("feather df shape", df.shape)
        df['filename'] = df["database_id"]
        df = df.set_index("database_id")
        # Prepend BIN_ to the columns
        new_names = {col : "BIN_" + str(col) for col in df.columns if col not in ["filename", "genus", "species"]}
        df.rename(columns=new_names, inplace=True)
        
        # TODO probably need genus and species columns
        for database_id in df.index:
            print("Yielding database_id", database_id, df.loc[[database_id]])
            yield database_id, df.loc[[database_id]]
    else:
        raise ValueError("Unsupported input_data type")

def database_search(input_paths, bin_size, database_df,
                    db_numerical_columns, seed_genera_indices, seed_species_indices, args):
    print("Running database search with the following parameters:")
    output_results_list = []
    complete_output_results = defaultdict(list)

    iterator = spectrum_iterator(input_paths, bin_size, args)

    for input_filename, spectra_binned_df in iterator:
        # Reading Query
        # ms1_df, _ = load_data(input_filename)
        
        # spectra_binned_df = spectrum_binner(ms1_df,
        #                                     input_filename,
        #                                     bin_size=bin_size,
        #                                     min_mz=args.mass_range_lower,
        #                                     max_mz=args.mass_range_upper,
        #                                     merge_replicates="Yes",)
        # print("standard spectra_binned_df cols", spectra_binned_df.columns)

        output_filename = os.path.join('query_spectra',  Path(input_filename).stem + ".mzML")
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
        logging.info("Computing distance matrix for query {} against database.".format(input_filename))
        start_time = time()
        distance_matrix = compute_distances_binned(query_data_np, database_data_np, distance_metric=args.distance, penalty=args.unmatched_peak_penalty)
        logging.info("Computed distance matrix for query {} against database in {:.2f} seconds".format(input_filename, time() - start_time))

        assert distance_matrix.shape[0] == len(spectra_binned_df), "Distance matrix rows should match the number of queries, but got {} and {}".format(distance_matrix.shape[0], len(spectra_binned_df))
        assert distance_matrix.shape[1] == len(spectra_binned_db_df), "Distance matrix columns should match the number of database entries, but got {} and {}".format(distance_matrix.shape[1], len(spectra_binned_db_df))

        for i, row in enumerate(distance_matrix):
            for j, item in enumerate(row):
                query_index = i
                database_index = j
                distance = item

                result_dict = {}
                result_dict["query_index"] = query_index
                result_dict["database_index"] = database_index
                result_dict["distance"] = distance
                result_dict["query_filename"] = os.path.basename(input_filename)
        
                complete_output_results[database_index].append(result_dict)
               
                if len(seed_genera_indices) == 0 and len(seed_species_indices) == 0:
                    if distance < args.score_threshold:
                        output_results_list.append(result_dict)
                else:
                    is_seed = False
                    if database_index in seed_genera_indices or \
                        database_index in seed_species_indices:
                        is_seed = True
                    if is_seed:
                        output_results_list.append(result_dict)
                    
    # We will also need a complete distance_matrix for the query-db pairs (in contrast to the above one that's trimmed by the score_threshold)
    matched_db_indices = np.unique([x["database_index"] for x in output_results_list])

    complete_output_results_list = []
    for idx in matched_db_indices:
        relevant_matches = complete_output_results[idx]
        complete_output_results_list.extend(relevant_matches)

    return output_results_list, complete_output_results_list

def within_query_distance(input_paths, bin_size, args):
    """ This function computes the pairwise distance between all queries in the input folder.
    
    Args:
    input_paths: list
        List of paths to the input files.
    bin_size: float
        Size of the spectra bins for distance calculations.
    args: argparse.Namespace
        Parsed command line arguments.
    Returns:
    output_results_df: pd.DataFrame
        DataFrame containing the pairwise distances between all queries.
    """
    logging.info("Computing pairwise distances between all {} queries in the input folder.".format(len(input_paths)))
    all_queries = []

    iterator = spectrum_iterator(input_paths, bin_size, args)

    for input_filename, spectra_binned_df in iterator:
        if 'scan' in spectra_binned_df.columns:
            spectra_binned_df.drop(['scan'], axis=1, inplace=True)
        all_queries.append(spectra_binned_df)

    logging.debug("Sample of queries: {}".format(all_queries[0]))

    # Concatenate all queries into a single DataFrame
    all_queries_df = pd.concat(all_queries, ignore_index=True)
    all_queries_df.set_index("filename", inplace=True)

    # Now lets do pairwise cosine distance
    query_data_np = all_queries_df.to_numpy()
    query_data_np = np.nan_to_num(query_data_np, 0.0)  # Fill NaNs with 0 for distance calculation  (means no peak)
    logging.debug("query_data_np shape: {}".format(query_data_np.shape))
    logging.debug("query_data_np: {}".format(query_data_np))
    logging.debug("query_data_np max: {}".format(np.max(query_data_np)))
    logging.debug("query_data_np min: {}".format(np.min(query_data_np)))
    logging.debug("query_data_np contains NaNs: {}".format(np.isnan(query_data_np).any()))

    # Query-Query distance
    start_time = time()
    distance_matrix = compute_distances_binned(query_data_np, distance_metric=args.distance, penalty=args.unmatched_peak_penalty)
    logging.info("Computed query to query distance matrix in {:.2f} seconds".format(time() - start_time))

    # Create a DataFrame from the distance matrix
    distance_df = pd.DataFrame(distance_matrix, index=all_queries_df.index, columns=all_queries_df.index)
    # Convert to adjacency list
    distance_df_melted = distance_df.melt(
        ignore_index=False, var_name="query_filename_right", value_name="distance"
    ).reset_index()
    distance_df_melted.columns = ["query_filename_left", "query_filename_right", "distance"]

    # For now, let's do a sanity check on random indices
    if len(distance_df_melted) > 0:
        vals_a = distance_df_melted["query_filename_left"].sample(5, replace=True).values
        vals_b = distance_df_melted["query_filename_right"].sample(5, replace=True).values

        for a, b in zip(vals_a, vals_b):
            # Check that the distance is symmetric
            if 'reverse' not in args.distance:
                assert distance_df_melted.loc[(distance_df_melted["query_filename_left"] == a) & (distance_df_melted["query_filename_right"] == b), "distance"].values[0] == \
                    distance_df_melted.loc[(distance_df_melted["query_filename_left"] == b) & (distance_df_melted["query_filename_right"] == a), "distance"].values[0], \
                    f"Distance between {a} and {b} is not the same as between {b} and {a} (not symmetric)."
                
            # Check that it's the same in the melted DataFrame
            assert distance_df_melted.loc[(distance_df_melted["query_filename_left"] == a) & (distance_df_melted["query_filename_right"] == b), "distance"].values[0] == \
                distance_df.loc[a, b], \
                f"Distance between {a} and {b} in melted DataFrame does not match the original DataFrame: {distance_df.loc[a, b]}"
                
            # Also compare against a manual distance calculation
            manual_distance = compute_distances_binned(np.nan_to_num(all_queries_df.loc[a].to_numpy()).reshape(1, -1),
                                                       np.nan_to_num(all_queries_df.loc[b].to_numpy()).reshape(1, -1),
                                                       distance_metric=args.distance,
                                                       penalty=args.unmatched_peak_penalty)
            assert np.isclose(distance_df_melted.loc[(distance_df_melted["query_filename_left"] == a) & (distance_df_melted["query_filename_right"] == b), "distance"].values[0],
                               manual_distance[0][0]), \
                f"Distance between {a} and {b} is not the same as the manual distance calculation: {manual_distance[0][0]}"
    else:
        logging.warning("No queries found in the input folder. Returning empty distance DataFrame.")

    logging.info("Done computing pairwise distances between all queries.")
    return distance_df_melted
    

def mock_output_results(args):
    # Write headers only
    with open(args.output_search_results_tsv, "w") as f:
        f.write("query_index\tdatabase_index\tdistance\tquery_filename\n")
    with open(args.complete_output_results_tsv, "w") as f:
        f.write("query_index\tdatabase_index\tdistance\tquery_filename\n")
    with open(args.output_db_db_distance_tsv, "w") as f:
        f.write("left_index\tright_index\tdatabase_id_left\tdatabase_id_right\tdatabase_scan_left\tdatabase_scan_right\tdistance\n")
    exit(0)

def main():
    
    args = parse_args()
    
    bin_size = float(args.bin_size)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
        logging.getLogger().setLevel(logging.DEBUG)  # Force set the level
        # Log all params
        for arg in vars(args):
            logging.debug("%s: %s", arg, getattr(args, arg))
    else:
        logging.basicConfig(level=logging.INFO)
        logging.getLogger().setLevel(logging.INFO)  # Force set the level

    # Print log level to sys.stderr
    print("Log level:", logging.getLevelName(logging.getLogger().level), file=os.sys.stderr)

    database_df, db_numerical_columns, db_metadata = load_and_prepare_database(args.database_filtered_json,
                                                                               args.mass_range_lower,
                                                                               args.mass_range_upper,
                                                                               bin_size)

    # Get the index of selected genera and species seeds.
    seed_genera_indices, seed_species_indices = get_seed_indices(db_metadata,
                                                                 seed_genera=args.seed_genera,
                                                                 seed_species=args.seed_species)

    if args.ml:
        all_input_files = glob.glob(os.path.join(args.input_folder, "*.feather"))
    else:
        all_input_files = glob.glob(os.path.join(args.input_folder, "*.mzML"))


    # Database search    
    logging.info("Starting database search for {} input files.".format(len(all_input_files)))
    output_results_list, complete_output_results_list = database_search(all_input_files, bin_size, database_df,
                                                                 db_numerical_columns, seed_genera_indices, seed_species_indices, args)

    assert len(output_results_list) <= len(complete_output_results_list), f"Output results list should be less than or equal to complete output results list but got {len(output_results_list)} and {len(complete_output_results_list)} respectively."

    small_database_df = database_df[["row_count", "database_id"]]
    small_database_df["database_scan"] = small_database_df["database_id"]

    # Query-Query Distance
    logging.info("Starting query-query distance calculation for {} input files.".format(len(all_input_files)))
    within_query_distance_df = within_query_distance(all_input_files, bin_size, args)
    within_query_distance_df.to_csv(args.output_query_query_distances_tsv, sep="\t", index=False)

    output_results_df = pd.DataFrame(output_results_list, columns=["query_index", "database_index", "distance", "query_filename"])
    complete_output_results_df = pd.DataFrame(complete_output_results_list, columns=["query_index", "database_index", "distance", "query_filename"])
    if len(output_results_df) == 0:
        print("No matches found")
        mock_output_results(args)    
        
    # Here we want to take the results from the search
    output_results_df = output_results_df.merge(small_database_df, left_on="database_index", right_on="row_count", how="left")
    output_results_df["database_id"] = output_results_df["database_scan"]
    output_results_df["database_scan"] = output_results_df["row_count"] + 1
    complete_output_results_df = complete_output_results_df.merge(small_database_df, left_on="database_index", right_on="row_count", how="left")
    complete_output_results_df["database_id"] = complete_output_results_df["database_scan"]
    complete_output_results_df["database_scan"] = complete_output_results_df["row_count"] + 1

    # Sanity Check: We assume these to be unique below
    if len(database_df.database_id.unique()) != len(database_df.database_id):
        raise ValueError("Database ID is not unique")
    
    # Perform pairwise distance of database result hits to fill out the remainder of the distance matrix.
    database_results_df = output_results_df[["database_id", "database_scan"]].drop_duplicates()

    # Add numerical columns to the database results, for the distance calculation
    logging.info("Merging database results with database for DB-DB distance calculation.")
    database_results_df = database_results_df.merge(database_df, left_on="database_id", right_on="database_id", how="left")

    # Compute DB-DB distances and save to file immediately
    if 'reverse' not in args.distance:
        logging.info("Starting database-database distance calculation for {} database results.".format(len(database_results_df)))
        compute_db_db_distance(database_results_df, db_numerical_columns, args.output_db_db_distance_tsv, distance_metric=args.distance, unmatched_peak_penalty=args.unmatched_peak_penalty)

    # Enrich the library information by hitting the web api
    all_database_metadata_json = requests.get("https://idbac.org/api/spectra").json()
    all_database_metadata_df = pd.DataFrame(all_database_metadata_json)
    all_database_metadata_df = all_database_metadata_df[["database_id", "Strain name", "Culture Collection", "Sample name", "Genbank accession"]]
    
    # lets merge the results with the metadata
    output_results_df = output_results_df.merge(all_database_metadata_df, left_on="database_id", right_on="database_id", how="left")
    complete_output_results_df = complete_output_results_df.merge(all_database_metadata_df, left_on="database_id", right_on="database_id", how="left")

    # rename columns
    output_results_df = output_results_df.rename(columns={"Strain name": "db_strain_name", "Culture Collection": "db_culture_collection", "Sample name": "db_sample_name", "Genbank accession": "db_genbank_accession"})
    complete_output_results_df = complete_output_results_df.rename(columns={"Strain name": "db_strain_name", "Culture Collection": "db_culture_collection", "Sample name": "db_sample_name", "Genbank accession": "db_genbank_accession"})

    # Output data
    output_results_df.to_csv(args.output_search_results_tsv, sep="\t", index=False)
    complete_output_results_df.to_csv(args.complete_output_results_tsv, sep="\t", index=False)


if __name__ == '__main__':
    main()