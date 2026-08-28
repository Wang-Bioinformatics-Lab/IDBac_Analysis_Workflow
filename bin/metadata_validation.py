import argparse
import pandas as pd

def metadata_validation(metadata_table:pd.DataFrame):
    """Validates the metadata table against the spectrum file. This function checks the following:
    1. There are no duplicate values in the "Filename" column. If there are, the fucntion will raise an error and halt execution.

    Any changes to this function should also be reflected in the IDBac Interactive Interface code.
    """

    if 'Filename' not in metadata_table.columns:
        raise ValueError("IDBAC_USER_ERROR: The metadata table is missing the required 'Filename' column.")

    # Check for duplicates in the metadata table
    duplicated_rows = metadata_table[metadata_table['Filename'].duplicated(keep=False)]
    if not duplicated_rows.empty:
        duplicate_names = duplicated_rows['Filename'].astype(str).tolist()
        raise ValueError(
            "IDBAC_USER_ERROR: The metadata table contains duplicate values in the 'Filename' column. "
            f"Please remove duplicates and try again. Duplicate filenames: {duplicate_names}"
        )



def main():
    parser = argparse.ArgumentParser(description='Check metadata table for common errors. If errors are found, output them to a csv and exit with 1.')
    parser.add_argument('--metadata_table', help='Path to the metadata table')
    args = parser.parse_args()
    
    metadata_table = pd.read_csv(args.metadata_table)
    
    metadata_validation(metadata_table)
