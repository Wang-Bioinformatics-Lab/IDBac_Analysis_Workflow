import argparse
import os
import sys
import numpy as np
from pyteomics import mzml
import csv
import re
import logging


def find_integer_at_end(string):
    return int(re.search(r'\d+$', string).group()) if re.search(r'\d+$', string) else 'N/A'


def validate_file(input_file:str, output_file:str)->int:
    """This function will validate the input file for common errors and output them to a csv file if errors are found."""
    
    # Check that the file is not empty
    if os.path.getsize(input_file) == 0:
        with open(output_file, 'w', encoding='utf-8') as output_csv:
            headers = ['error_level', 'scan', 'original_filename', 'error']
            output_writer = csv.DictWriter(output_csv, fieldnames=headers)
            output_writer.writerow({'error_level': 'critical', 
                                    'scan': 'N/A',
                                    'original_filename': input_file, 
                                    'error': 'File is empty. If the file was uploaded, please check the uploaded file to ensure it was successful.'})
        return 1
    
    output_status = 0
    with mzml.MzML(input_file) as reader:
        with open(output_file, 'w', encoding='utf-8') as output_csv:
            headers = ['error_level', 'scan', 'original_filename', 'error']
            output_writer = csv.DictWriter(output_csv, fieldnames=headers)           
            
            # Check to make sure scans are readable
            if len(reader) == 0:
                output_writer.writerow({'original_filename': input_file, 'error': 'No scans found'})
                output_status = 1
            # Check to make sure scans have nonzero intensity arrays
            scan_status = np.zeros(len(reader), dtype=bool)
            for scan_idx, scan in enumerate(reader):
                if 'intensity array' not in scan:
                    output_writer.writerow({'error_level': 'warning',
                                            'scan': find_integer_at_end(scan['id']),
                                            'original_filename': input_file,
                                            'error': f"Scan {scan['id']} is missing intensity array"})
                    output_status = 1
                    scan_status[scan_idx] = True
                    logging.info(f"The following scan is missing intensity array: {scan['id']}")
                    logging.info(str(scan))

                if len(scan['intensity array']) == 0:
                    output_writer.writerow({'error_level': 'warning',
                                            'scan': find_integer_at_end(scan['id']),
                                            'original_filename': input_file,
                                            'error': f"Scan {scan['id']} has empty intensity array"})
                    output_status = 1
                    scan_status[scan_idx] = True
                    logging.info(f"The following scan reported empty intensity array: {scan['id']}")
                    logging.info(str(scan))

                if max(scan['intensity array']) == 0:
                    output_writer.writerow({'error_level': 'warning',
                                            'scan': find_integer_at_end(scan['id']),
                                            'original_filename': input_file,
                                            'error': f"Scan {scan['id']} contains no peaks"})
                    output_status = 1
                    scan_status[scan_idx] = True
                    logging.info(f"The following scan reported zero intensity: {scan['id']}")
                    logging.info(str(scan))

            if np.all(scan_status):
                output_writer.writerow({'error_level': 'critical',
                                        'scan': 'N/A',
                                        'original_filename': input_file,
                                        'error': 'All scans in this file are unreadable or empty. Please check the file and try again.'})
                output_status = 1
            
    return output_status

def main():
    parser = argparse.ArgumentParser(description='Check mzML files for common errors. If errors are found, output them to a csv and exit with 1.')
    parser.add_argument('--input_file', help='Path to the input file')
    parser.add_argument('--output_file', help='Path to the output file')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    for var, value in vars(args).items():
        logging.info(f'Argument {var}: {value}')
    
    status = validate_file(args.input_file, args.output_file)
    # sys.exit(status) # Nextflow doesn't have a provision to output files if the process fails. If this comes in the future, it will be a good way to warn users

if __name__ == "__main__":
    main()