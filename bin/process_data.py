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
    parser.add_argument('output_filename')

    args = parser.parse_args()

    # Lets read all the spectra that are coming out of the input_folder
    all_input_files = glob.glob(os.path.join(args.input_folder, "*.mzML"))

    all_spectra_df_list = []

    for input_filename in all_input_files:
        print("Loading data from {}".format(input_filename))
        ms1_df, ms2_df = load_data(input_filename)

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
        spectra_binned_df["filename"] = input_filename

        all_spectra_df_list.append(spectra_binned_df)

        # DEBUG SKIP
        # break

    all_spectra_df = pd.concat(all_spectra_df_list)

    numerical_columns = [x for x in all_spectra_df.columns if x.startswith("BIN_")]

    # Fill in the missing values with 0
    all_spectra_df[numerical_columns] = all_spectra_df[numerical_columns].fillna(0)

    data_np = all_spectra_df[numerical_columns].to_numpy()

    # lets convert each row to a numpy array
    print(all_spectra_df)
    print(data_np)

    # Now lets do pairwise cosine similarity
    from sklearn.metrics.pairwise import cosine_similarity
    similarity_matrix = cosine_similarity(data_np)

    print(similarity_matrix)

    # Lets make this into a dendrogram
    import plotly.figure_factory as ff

    # Creating labels
    all_spectra_df["label"] = all_spectra_df["filename"].apply(lambda x: os.path.basename(x)) + ":" + all_spectra_df["scan"].astype(str)
    all_labels_list = all_spectra_df["label"].to_list()

    dendro = ff.create_dendrogram(similarity_matrix, orientation='left', labels=all_labels_list)
    dendro.update_layout(width=800, height=15*len(all_labels_list))
    dendro.write_html(args.output_filename)


    




if __name__ == '__main__':
    main()