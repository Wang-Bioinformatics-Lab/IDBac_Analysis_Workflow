import argparse
import os
from pathlib import Path

from utils import load_metadata_file

def main():
    parser = argparse.ArgumentParser(description="Preprocess metadata file")
    parser.add_argument("--output_histogram_data_directory", type=str, help="Path to the output histogram data directory")
    parser.add_argument("--metadata_path", type=str, help="Path to the metadata file")
    args = parser.parse_args()
    
    output_histogram_data_directory = Path(args.output_histogram_data_directory)

    if not output_histogram_data_directory.exists():
        output_histogram_data_directory.mkdir(parents=True, exist_ok=True)

    metadata_path = Path(args.metadata_path)
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found at {metadata_path}")

    metadata_df = load_metadata_file(metadata_path)
    if metadata_df is not None:
        metadata_df.to_csv((output_histogram_data_directory / "metadata.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()