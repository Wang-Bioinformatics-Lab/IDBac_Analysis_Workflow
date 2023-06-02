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

def parse_metadata(metadata_filename):
    # check extension
    if metadata_filename.endswith(".xlsx"):
        metadata_df = pd.read_excel(metadata_filename)
    elif metadata_filename.endswith(".csv"):
        metadata_df = pd.read_csv(metadata_filename)
    elif metadata_filename.endswith(".tsv"):
        metadata_df = pd.read_csv(metadata_filename, sep="\t")
    else:
        metadata_df = pd.read_csv(metadata_filename, sep=None)

    return metadata_df

    

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input_folder')
    parser.add_argument('metadata_filename') # We will use this metadata to paint the plots better
    parser.add_argument('output_basic_html_plot')
    parser.add_argument('output_metadata_html_plot')
    parser.add_argument('--merge_replicates', default="No")
    parser.add_argument('--similarity', default="cosine")
    parser.add_argument('--metadata_column', default="None")

    args = parser.parse_args()

    metadata_df = parse_metadata(args.metadata_filename)

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

        all_spectra_df_list.append(spectra_binned_df)


    all_spectra_df = pd.concat(all_spectra_df_list)

    numerical_columns = [x for x in all_spectra_df.columns if x.startswith("BIN_")]

    # Fill in the missing values with 0
    all_spectra_df[numerical_columns] = all_spectra_df[numerical_columns].fillna(0)

    data_np = all_spectra_df[numerical_columns].to_numpy()

    # Now lets do pairwise cosine similarity
    from sklearn.metrics.pairwise import cosine_distances, euclidean_distances

    selected_distance_fun = None
    if args.similarity == "cosine":
        selected_distance_fun = cosine_distances
    elif args.similarity == "euclidean":
        selected_distance_fun = euclidean_distances
    elif args.similarity == "presence":
        selected_distance_fun = cosine_distances
        
        # Update the data to be 1 or 0
        data_np[data_np > 0] = 1

    # Lets make this into a dendrogram
    import plotly.figure_factory as ff

    # Creating labels
    all_spectra_df["label"] = all_spectra_df["filename"].apply(lambda x: os.path.basename(x)) + ":" + all_spectra_df["scan"].astype(str)
    all_labels_list = all_spectra_df["label"].to_list()

    dendro = ff.create_dendrogram(data_np, orientation='left', labels=all_labels_list, distfun=selected_distance_fun)
    dendro.update_layout(width=800, height=max(15*len(all_labels_list), 150))
    dendro.write_html(args.output_basic_html_plot)

    # Making using metadata
    try:
        # Selecting the column to visualize
        # coloring_column = "Strain / Isolate Source"
        coloring_column = args.metadata_column

        all_spectra_df = all_spectra_df.merge(metadata_df, how="left", left_on="filename", right_on="Filename")
        all_spectra_df["label"] = all_spectra_df[coloring_column]
        all_labels_list = all_spectra_df["label"].to_list()

        # Define color palette
        color_palette = ['red', 'green', 'blue', 'yellow']

        all_categories = list(all_spectra_df[coloring_column])
        unique_categories = list(set(all_categories))
        leaf_colors = [color_palette[unique_categories.index(category)] if unique_categories.index(category) < len(color_palette) else 'gray' for category in all_categories]

        # Creating Dendrogram
        dendro = ff.create_dendrogram(data_np, orientation='left', labels=all_labels_list, distfun=selected_distance_fun, colorscale=leaf_colors)
        dendro.update_layout(width=800, height=max(15*len(all_labels_list), 150))

        # Setting the leaf colors
        dendro['data'][0]['textfont']['color'] = leaf_colors
        
        # Create legend annotations
        import plotly.graph_objects as go

        legend_annotations = []
        for i, category in enumerate(unique_categories):
            annotation = go.layout.Annotation(
                x=1,
                y=i,
                xref='paper',
                yref='paper',
                text=category,
                showarrow=False,
                font=dict(color=color_palette[unique_categories.index(category)] if unique_categories.index(category) < len(color_palette) else 'gray')
            )
            legend_annotations.append(annotation)

        # Create legend shapes
        legend_shapes = []
        for i in range(len(unique_categories)):
            shape = go.layout.Shape(
                type='rect',
                x0=1.05,
                y0=1 - (i + 1) * 0.05,
                x1=1.1,
                y1=1 - i * 0.05,
                fillcolor=color_palette[unique_categories.index(category)] if unique_categories.index(category) < len(color_palette) else 'gray',
                line=dict(width=0)
            )
            legend_shapes.append(shape)


        # Add legend annotations and shapes to the layout
        dendro.update_layout(
            annotations=legend_annotations,
            #shapes=legend_shapes,
            showlegend=False  # Disable default legend
        )

        dendro.write_html(args.output_metadata_html_plot)

    except:
        pass

if __name__ == '__main__':
    main()