import argparse
import sys
from pyteomics import mzml
import csv

def validate_file(input_file:str, output_file:str)->int:
    """This function will validate the input file for common errors and output them to a csv file if errors are found."""
    output_status = 0
    with mzml.MzML(input_file) as reader:
        with open(output_file, 'w', encoding='utf-8') as output_csv:
            headers = ['original_filename', 'error']
            output_writer = csv.DictWriter(output_csv, fieldnames=headers)           
            
            # Check to make sure scans are readable
            if len(reader) == 0:
                output_writer.writerow({'original_filename': input_file, 'error': 'No scans found'})
                output_status = 1
            # Check to make sure scans have nonzero intensity arrays
            for scan in reader:
                if 'intensity array' not in scan:
                    output_writer.writerow({'original_filename': input_file, 'error': f"Scan {scan['id']} is missing intensity array"})
                    output_status = 1
                if len(scan['intensity array']) == 0:
                    output_writer.writerow({'original_filename': input_file, 'error': f"Scan {scan['id']} has empty intensity array"})
                    output_status = 1
                if max(scan['intensity array']) == 0:
                    output_writer.writerow({'original_filename': input_file, 'error': f"Scan {scan['id']} contains no peaks"})
                    output_status = 1
            
    return output_status

def main():
    parser = argparse.ArgumentParser(description='Check mzML files for common errors. If errors are found, output them to a csv and exit with 1.')
    parser.add_argument('--input_file', help='Path to the input file')
    parser.add_argument('--output_file', help='Path to the output file')
    args = parser.parse_args()
    
    status = validate_file(args.input_file, args.output_file)
    # sys.exit(status) # Nextflow doesn't have a provision to output files if the process fails. If this comes in the future, it will be a good way to warn users

if __name__ == "__main__":
    main()