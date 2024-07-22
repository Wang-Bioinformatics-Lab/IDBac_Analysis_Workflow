import argparse
import os
import pandas as pd
from pyteomics import mzxml, mzml
from psims.mzml.writer import MzMLWriter
# from merge_spectra import load_data

def load_data(input_filepath):
    
    ms1_df_list = []
    ms2_df_list = []
    with mzml.MzML(input_filepath) as reader:
        for spectrum in reader:
            if spectrum['ms level'] == 1:
                ms1_df_list.append(pd.DataFrame({'mz': spectrum['m/z array'], 'i': spectrum['intensity array'], 'scan': spectrum['id']}))
            elif spectrum['ms level'] == 2:
                ms2_df_list.append(pd.DataFrame({'mz': spectrum['m/z array'], 'i': spectrum['intensity array'], 'scan': spectrum['id']}))
            else:
                raise ValueError(f"Unsupported MS level {spectrum['ms level']}")
            
    ms1_df = pd.concat(ms1_df_list, ignore_index=True)
    if len(ms2_df_list) > 0:
        ms2_df = pd.concat(ms2_df_list, ignore_index=True)
    else:
        ms2_df = None
    return ms1_df, ms2_df
    

def media_control(mzml_file, media_control_file, output_file):
    """This function will iterate over the media control file, collecting all peaks that exist in all scans 
    (discretized by 0.001 m/z) and remove them from the mzML file."""
    
    media_control_df, _ = load_data(media_control_file)
 
    ## DEBUG
    print("Initial media control df")
    print(media_control_df)
    
    # Discretize mz values to 3 decimal places
    media_control_df['mz'] = media_control_df['mz'].round(3).astype(str)    # Convert to string to groupby
    # Sum intensities for descritized m/z values in the same scan
    media_control_df = media_control_df.groupby(['scan', 'mz']).sum().reset_index()
    # Now that we've summed intensities we can convert mz back to float
    media_control_df['mz'] = media_control_df['mz'].astype(float)
    # Normalize intensities per scan (i.e., divide by max per scan)
    media_control_df['i'] = media_control_df.groupby('scan')['i'].transform(lambda x: x / x.max())
    # Remove peaks less than 0.5 intensity
    org_len = len(media_control_df)
    media_control_df = media_control_df.loc[media_control_df.i > 0.05]
    print(f"Removed {org_len - len(media_control_df)} peaks with intensity less than 0.05 intensity")
    
    # Get m/z values that occur in all scans
    media_control_df = media_control_df.groupby('mz').filter(lambda x: len(x) == media_control_df['scan'].nunique())

    ## DEBUG
    print("Final media control df")
    print(media_control_df)

    with mzml.MzML(mzml_file) as mzml_reader:
        with MzMLWriter(open(output_file, 'wb'), close=True) as out:
                # Add default controlled vocabularies
                out.controlled_vocabularies()
                # Open the run and spectrum list sections
                with out.run(id="my_analysis"):
                    spectrum_count = len(mzml_reader)
                    with out.spectrum_list(count=spectrum_count):
                        scan = 1
                        for spectrum in mzml_reader:

                            assert spectrum['ms level'] == 1, "Only MS1 spectra are supported"
                            mz_array = spectrum['m/z array']
                            intensity_array = spectrum['intensity array']
                            
                            # Remove peaks that are in the media control file
                            org_len = len(mz_array)
                            indices = [i for i, mz in enumerate(mz_array) if not any(abs(mz - x) < 0.001 for x in media_control_df["mz"])]
                            mz_array = [mz_array[i] for i in indices]
                            intensity_array = [intensity_array[i] for i in indices]
                            print(f"Removed {org_len - len(mz_array)} peaks from scan {scan} that were in the media control file")

                            # Write scan
                            out.write_spectrum(
                                mz_array, intensity_array,
                                id=f"scan={scan}", params=[
                                    "MS1 Spectrum",
                                    {"ms level": 1},
                                    {"total ion current": sum(intensity_array)}
                                ])
                            
                            scan += 1

def main():
    parser = argparse.ArgumentParser(description='Remove Media and Control Samples from mzML files')
    parser.add_argument('--small_molecule_file', help='Path to the small molecule file')
    parser.add_argument('--metadata_file', help='Path to the metadata file')
    parser.add_argument('--media_control_dir', help='Path to the media control dir')
    parser.add_argument('--output_file', help='Path to the output file')
    args = parser.parse_args()
    
    # Create output folder
    output_folder = os.path.dirname(args.output_file)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)
        
    if args.metadata_file.endswith('.csv'):
        metadata_df = pd.read_csv(args.metadata_file)
    elif args.metadata_file.endswith('.xlsx'):
        metadata_df = pd.read_excel(args.metadata_file)
        # If it contains multiple tables, get the one named "Metadata sheet"
        if isinstance(metadata_df, dict):
            metadata_df = metadata_df['Metadata sheet']
    elif args.metadata_file.endswith('.xls'):
        metadata_df = pd.read_excel(args.metadata_file)
        # If it contains multiple tables, get the one named "Metadata sheet"
        if isinstance(metadata_df, dict):
            metadata_df = metadata_df['Metadata sheet']
    elif args.metadata_file.endswith('.tsv'):
        metadata_df = pd.read_csv(args.metadata_file, sep='\t')
    else:
        raise ValueError(f'Metadata file must be a CSV, XLSX, XLS, or TSV file, but got {args.metadata_file} instead.')
    
    # Since sometimes the metadata has asterisks in the 'Scan/Coordinate' column that extend beyond the remainder of the values, we need to remove these rows
    metadata_df = metadata_df.dropna(subset=[x for x in metadata_df.columns if 'Scan/Coordinate' != x], how='all')
    
    # Get relevant media control file for the mzML file
    try:
        media_control_file = metadata_df[metadata_df['Small molecule file name'] == os.path.basename(args.small_molecule_file)]['Blank filename'].values[0] # Column names subject to change
    except IndexError:
        print(metadata_df)
        print(metadata_df[metadata_df['Small molecule file name'] == os.path.basename(args.small_molecule_file)])
        raise ValueError(f"Could not find media control file for '{os.path.basename(args.small_molecule_file)}' Is it in the metadata file?")
       
    media_control(mzml_file=args.small_molecule_file, media_control_file=os.path.join(args.media_control_dir,media_control_file), output_file=args.output_file)

if __name__ == "__main__":
    main()