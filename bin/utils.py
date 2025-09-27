import os
import numpy as np
import pandas as pd
from tqdm import tqdm

from massql import msql_fileloading
from pyteomics import mzxml, mzml
from psims.mzml.writer import MzMLWriter

from sklearn.metrics.pairwise import cosine_distances, euclidean_distances

import logging
import yaml

from pathlib import Path

def load_data(input_filename):
    """
    Reads an mzML file and converts it to a pandas dataframe. 
    
    Args:
    input_filename: str, path to the mzML file
    
    Returns:
    ms1_df: pd.DataFrame, MS1 data
    ms2_df: pd.DataFrame, MS2 data
    """
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
                scan = spectrum["id"].replace("scanId=", "").split("scan=")[-1]+f"_{spectrum['index']}"
            except:
                scan = spectrum["id"] + str(spectrum['index'])

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

def peak_filtering(database_df, relative_intensity):
    """ Perform any instrument-specific postprocessing on merged spectra. 
    Current operations include: relative intensity peak filtering based on instrument type.
    
    Args:
    database_df: pd.DataFrame, binned MS1 data
    relative_intensity: float, relative intensity threshold (between 0 and 1)

    Returns:
    database_df: pd.DataFrame, filtered binned MS1 data
    """

    bin_cols = [x for x in database_df.columns if x.startswith("BIN_")]

    if relative_intensity <= 0.0:
        return database_df

    def _helper(row):
        """Set all bin_cols to zero if they are below the relative intensity threshold
        for the given instrument

        Inputs:
            row: a row of the dataframe
        Returns:
            row: the modified row
        """
        max_intensity = row[bin_cols].max()
        if pd.isna(max_intensity) or max_intensity == 0:
            # If all intensities are zero, do nothing
            return row

        relative_thresh = max_intensity * relative_intensity
        row[bin_cols] = row[bin_cols].apply(lambda x: x if x >= relative_thresh else 0)

        return row        

    database_df = database_df.apply(_helper, axis=1)

    return database_df

def spectrum_binner(ms1_df:pd.DataFrame, input_filename:str, bin_size=1.0, min_mz=2000.0, max_mz=20000.0, merge_replicates="No"):
    """
    Bins MS1 dataframe into a 1d vector that is the intensity value for each bin
    
    Args:
    ms1_df: pd.DataFrame, MS1 data
    input_filename: str, path to the mzML file
    bin_size: float, size of the bin
    min_mz: float, minimum m/z value to consider
    max_mz: float, maximum m/z value to consider
    merge_replicates: str, whether to merge replicates or not
    
    Returns:
    spectra_binned_df: pd.DataFrame, binned MS1 data
    """
    # Filtering m/z
    ms1_df = ms1_df.loc[ms1_df['mz'] > min_mz]
    ms1_df = ms1_df.loc[ms1_df['mz'] < max_mz]

    # Bin the MS1 Data by m/z within each spectrum
    ms1_df['bin'] = (ms1_df['mz'] / bin_size).astype(int)   # Since we have a min_mz, the bins won't start at zero, but this isn't an issue

    # Now we need to group by scan and bin
    ms1_df = ms1_df.groupby(['scan', 'bin']).agg({'i': 'sum'}).reset_index()
    ms1_df["mz"] = ms1_df["bin"] * bin_size
    ms1_df["bin_name"] = "BIN_" + ms1_df["bin"].astype(str)
    
    # Turning each scan into a 1d vector that is the intensity value for each bin
    spectra_binned_df = ms1_df.pivot(index='scan', columns='bin_name', values='i').reset_index()
    spectra_binned_df["filename"] = os.path.basename(input_filename)

    bins_to_remove = []
    # merging replicates
    if merge_replicates == "Yes":
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
        spectra_binned_df = spectra_binned_df.drop("scan", axis=1).groupby("filename").mean().reset_index()
        spectra_binned_df["scan"] = "merged"
    
    return spectra_binned_df

def compute_distances_binned(np_data_X:np.ndarray, np_data_Y:np.ndarray=None, distance_metric:str='cosine'):
    
    if distance_metric not in ['cosine', 'euclidean', 'presence']:
        raise ValueError(f'Invalid distance metric. Expected "cosine", "euclidean", or "presence", but got {distance_metric}')
    
    if distance_metric == "cosine":
        selected_distance_fun = cosine_distances
    elif distance_metric == "euclidean":
        selected_distance_fun = euclidean_distances
    elif distance_metric == "presence":
        selected_distance_fun = cosine_distances
        np_data_X[np_data_X > 0] = 1
        if np_data_Y is not None:
            np_data_Y[np_data_Y > 0] = 1
            
        
    return selected_distance_fun(np_data_X, np_data_Y)


def load_metadata_file(metadata_path:str):
    """
    Reads a metadata file and converts it to a pandas dataframe. 
    
    Args:
    metadata_path: str, path to the metadata file

    Returns:
    metadata_df: pd.DataFrame, metadata

    Raises:
    ValueError: if the metadata file is not a CSV, XLSX, XLS, or TSV file
    """
    metadata_path = str(metadata_path)
    if metadata_path.endswith('.csv'):
        metadata_df = pd.read_csv(metadata_path)
    elif metadata_path.endswith('.xlsx'):
        metadata_df = pd.read_excel(metadata_path, sheet_name=None)
        # If it contains multiple tables, get the one named "Metadata sheet"
        if isinstance(metadata_df, dict):
            metadata_df = pd.read_excel(metadata_path, sheet_name=None)
        if 'Metadata sheet' in metadata_df:
            metadata_df = metadata_df['Metadata sheet']
        elif 'Metadata template' in metadata_df:
            metadata_df = metadata_df['Metadata template']
        else:
            # Pop "instructions" or Instructions or anything like that
            all_keys = list(metadata_df.keys())
            for key in all_keys:
                if key.lower().strip() in ['instructions', 'instruction', 'metadata instructions']:
                    metadata_df.pop(key)
            # If there is only one sheet, use that
            if len(metadata_df) == 1:
                metadata_df = list(metadata_df.values())[0]
            else:
                raise ValueError(f"Excel file should contain only one sheet, or one named 'Metadata sheet' or 'Metadata template'. Instead found {metadata_df.keys()}")
    elif metadata_path.endswith('.xls'):
        metadata_df = pd.read_excel(metadata_path, sheet_name=None)
        # If it contains multiple tables, get the one named "Metadata sheet"
        if isinstance(metadata_df, dict):
            metadata_df = pd.read_excel(metadata_path, sheet_name=None)
        if 'Metadata sheet' in metadata_df:
            metadata_df = metadata_df['Metadata sheet']
        elif 'Metadata template' in metadata_df:
            metadata_df = metadata_df['Metadata template']
        else:
            all_keys = list(metadata_df.keys())
            for key in all_keys:
                if key.lower().strip() in ['instructions', 'instruction', 'metadata instructions']:
                    metadata_df.pop(key)
            # If there is only one sheet, use that
            if len(metadata_df) == 1:
                metadata_df = list(metadata_df.values())[0]
            else:
                raise ValueError(f"Excel file should contain only one sheet, or one named 'Metadata sheet' or 'Metadata template'. Instead found {metadata_df.keys()}")
    elif metadata_path.endswith('.tsv'):
        metadata_df = pd.read_csv(metadata_path, sep='\t')
    else:
        raise ValueError(f'Metadata file must be a CSV, XLSX, XLS, or TSV file, but got {metadata_path} instead.')
    
    return metadata_df

def write_spectra_df_to_mzML(output_filename:str, spectra_binned_df:pd.DataFrame, bin_size):
    """The function writes a dataframe of binned spectra to an mzML file.

    Args:
    output_filename: str, path to the output mzML file
    spectra_binned_df: pd.DataFrame, binned MS1 data
    bin_size: float, size of the bin

    Returns:
    None
    """
    output_filename = Path(output_filename)
    output_filename.parent.mkdir(parents=True, exist_ok=True)
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