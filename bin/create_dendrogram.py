import sys
import os
import argparse
import pandas as pd
import uuid
import json
import numpy as np

from utils import load_data, spectrum_binner, compute_distances_binned, load_metadata_file
import logging

from tqdm import tqdm
import glob

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input_folder')
    parser.add_argument('metadata_filename') # We will use this metadata to paint the plots better
    parser.add_argument('input_database_hits')
    
    parser.add_argument('output_basic_html_plot')
    parser.add_argument('output_metadata_html_plot')
    parser.add_argument('output_db_html_plot')
    parser.add_argument('output_distance_table')

    # These are outputs for recreating the histogram
    parser.add_argument('output_histogram_data_directory')

    parser.add_argument('--merge_replicates', default="No")
    parser.add_argument('--distance', default="cosine")
    parser.add_argument('--metadata_column', default="None")
    parser.add_argument('--bin_size', default=10.0, type=float, help="Size of the spectra bins for distance calculations.")
    parser.add_argument('--mass_range_lower', default=2000.0, type=float, help="Minimum m/z value to consider for binning.")
    parser.add_argument('--mass_range_upper', default=20000.0, type=float, help="Maximum m/z value to consider for binning.")
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
        # Log all params
        for arg in vars(args):
            logging.debug("%s: %s", arg, getattr(args, arg))
    else:
        logging.basicConfig(level=logging.INFO)

    try:
        metadata_df = load_metadata_file(args.metadata_filename)
    except:
        metadata_df = None

    # Lets read all the spectra that are coming out of the input_folder
    all_input_files = glob.glob(os.path.join(args.input_folder, "*.mzML"))
    all_input_files.sort()

    all_spectra_df_list = []

    for input_filename in all_input_files:
        print("Loading data from {}".format(input_filename))
        ms1_df, _ = load_data(input_filename)
        logging.debug("Loaded data from {}".format(input_filename))
        logging.debug("MS1 shape: {}".format(ms1_df.shape))

        spectra_binned_df = spectrum_binner(ms1_df,
                                            input_filename,
                                            bin_size=args.bin_size,
                                            min_mz=args.mass_range_lower,
                                            max_mz=args.mass_range_upper,
                                            merge_replicates=args.merge_replicates)
        logging.debug("Binned data shape: {}".format(spectra_binned_df.shape))

        all_spectra_df_list.append(spectra_binned_df)


    all_spectra_df = pd.concat(all_spectra_df_list)
    logging.debug("All spectra shape: {}".format(all_spectra_df.shape))

    numerical_columns = [x for x in all_spectra_df.columns if x.startswith("BIN_")]
    logging.debug("Got {} numerical columns".format(len(numerical_columns)))

    # Fill in the missing values with 0
    all_spectra_df[numerical_columns] = all_spectra_df[numerical_columns].fillna(0)
    logging.debug("all_spectra_df {}".format(all_spectra_df))

    data_np = all_spectra_df[numerical_columns].to_numpy()
    logging.debug("data_np shape: {}".format(data_np.shape))
    logging.debug("data_np: {}".format(data_np))

    # Creating labels
    all_spectra_df["label"] = all_spectra_df["filename"].apply(lambda x: os.path.basename(x)) + ":" + all_spectra_df["scan"].astype(str)
    all_spectra_df["label"] = all_spectra_df["label"].apply(lambda x: x.replace(":merged", ""))
    all_labels_list = all_spectra_df["label"].to_list()

    output_scores_list = []
    if len(data_np) > 1:
        # Calculating the distances between all the spectra
        distance_matrix = compute_distances_binned(data_np, distance_metric=args.distance)
        logging.debug("Distance matrix shape: {}".format(distance_matrix.shape))
        logging.debug("Distance matrix: {}".format(distance_matrix))

        # Merge in the labels
        for index_i, label_i in enumerate(all_labels_list):
            for index_j, label_j in enumerate(all_labels_list):
                if index_i == index_j:
                    continue

                if index_i > index_j:
                    continue
                
                output_dict = {}
                output_dict["label_i"] = label_i
                output_dict["label_j"] = label_j
                output_dict["distance"] = distance_matrix[index_i][index_j]

                output_scores_list.append(output_dict)

    output_scores_df = pd.DataFrame(output_scores_list, columns=['label_i', 'label_j', 'distance'])
    output_scores_df.to_csv(args.output_distance_table, sep="\t", index=False)

    # Lets make this into a dendrogram
    import plotly.figure_factory as ff

    from sklearn.metrics.pairwise import cosine_distances, euclidean_distances
    if args.distance == "cosine":
        selected_distance_fun = cosine_distances
    elif args.distance == "euclidean":
        selected_distance_fun = euclidean_distances
    elif args.distance == "presence":
        selected_distance_fun = cosine_distances
        data_np[data_np > 0] = 1

    if data_np.shape[0] > 1:
        # Otherwise, there's no distances to compute
        dendro = ff.create_dendrogram(data_np, orientation='left', labels=all_labels_list, distfun=selected_distance_fun)
        dendro.update_layout(width=800, height=max(15*len(all_labels_list), 350))
        dendro.write_html(args.output_basic_html_plot)

    # Saving the output data to create the dendrograms
    output_numerical_spectra_filename = os.path.join(args.output_histogram_data_directory, "numerical_spectra.npy")
    output_labels_spectra_filename = os.path.join(args.output_histogram_data_directory, "labels_spectra.tsv")

    # Saving the data
    with open(output_numerical_spectra_filename, "wb") as output_file:
        np.save(output_file, data_np)

    # Saving the labels
    all_spectra_df.to_csv(output_labels_spectra_filename, sep="\t", index=False)


    # Making using metadata
    try:
        # Selecting the column to visualize
        # coloring_column = "Strain / Isolate Source"
        coloring_column = args.metadata_column

        all_spectra_df = all_spectra_df.merge(metadata_df, how="left", left_on="filename", right_on="Filename")
        all_spectra_df["label"] = all_spectra_df[coloring_column].fillna("No Metadata")
        all_spectra_df["label"] = all_spectra_df["label"] + " - " + all_spectra_df["filename"]
        all_labels_list = all_spectra_df["label"].to_list()

        # Define color palette
        color_palette = ['red', 'green', 'blue', 'yellow']

        all_categories = list(all_spectra_df[coloring_column])
        unique_categories = list(set(all_categories))
        leaf_colors = [color_palette[unique_categories.index(category)] if unique_categories.index(category) < len(color_palette) else 'gray' for category in all_categories]

        # Creating Dendrogram
        dendro = ff.create_dendrogram(data_np, orientation='left', labels=all_labels_list, distfun=selected_distance_fun, colorscale=leaf_colors)
        dendro.update_layout(width=800, height=max(15*len(all_labels_list), 350))

        # Setting the leaf colors
        # dendro['data'][0]['textfont']['color'] = leaf_colors
        
        # Create legend annotations
        # import plotly.graph_objects as go

        # legend_annotations = []
        # for i, category in enumerate(unique_categories):
        #     annotation = go.layout.Annotation(
        #         x=1,
        #         y=i,
        #         xref='paper',
        #         yref='paper',
        #         text=category,
        #         showarrow=False,
        #         font=dict(color=color_palette[unique_categories.index(category)] if unique_categories.index(category) < len(color_palette) else 'gray')
        #     )
        #     legend_annotations.append(annotation)

        # # Create legend shapes
        # legend_shapes = []
        # for i in range(len(unique_categories)):
        #     shape = go.layout.Shape(
        #         type='rect',
        #         x0=1.05,
        #         y0=1 - (i + 1) * 0.05,
        #         x1=1.1,
        #         y1=1 - i * 0.05,
        #         fillcolor=color_palette[unique_categories.index(category)] if unique_categories.index(category) < len(color_palette) else 'gray',
        #         line=dict(width=0)
        #     )
        #     legend_shapes.append(shape)


        # # Add legend annotations and shapes to the layout
        # dendro.update_layout(
        #     annotations=legend_annotations,
        #     #shapes=legend_shapes,
        #     showlegend=False  # Disable default legend
        # )

        dendro.write_html(args.output_metadata_html_plot)

    except:
        pass

    

    # Making using database hits
    try:
        input_database_hits_df = pd.read_csv(args.input_database_hits, sep="\t")
        
        # sort by distance
        input_database_hits_df = input_database_hits_df.sort_values(by="distance", ascending=False)

        # Picking the best one
        input_database_hits_df = input_database_hits_df.groupby("query_filename").first().reset_index()

        all_spectra_df = all_spectra_df.merge(input_database_hits_df, how="left", left_on="filename", right_on="query_filename")
        all_spectra_df["label"] = all_spectra_df["filename"]
        all_labels_list = all_spectra_df["label"].to_list()

        refined_labels_list = []
        for label in all_labels_list:
            try:
                refined_labels_list.append(label.split("; ")[-2] + "; " + label.split("; ")[-1])
            except:
                refined_labels_list.append(label)

        # Creating Dendrogram
        dendro = ff.create_dendrogram(data_np, orientation='left', labels=refined_labels_list, distfun=selected_distance_fun)
        dendro.update_layout(width=800, height=max(15*len(all_labels_list), 350))

        dendro.write_html(args.output_db_html_plot)

    except:
        pass

if __name__ == '__main__':
    main()