import argparse
import pathlib
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Merge a csv by columns.')
    parser.add_argument('--input_folder', type=str, required=True, help='The folder containing the csv files to merge.')
    parser.add_argument('--output_folder', type=str, required=True, help='The folder to save the merged csv file.') 
    args = parser.parse_args()

    input_folder = pathlib.Path(args.input_folder)
    output_folder = pathlib.Path(args.output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    all_input_files = list(input_folder.glob("*.csv"))

    df = pd.read_csv(all_input_files[0], index_col=0)

    for input_file in all_input_files[1:]:
        _df = pd.read_csv(input_file, index_col=0)
        df = df.merge(_df, left_index=True, right_index=True, how='outer')
        

    df.to_csv(output_folder / "merged.csv")

if __name__ == "__main__":
    main()