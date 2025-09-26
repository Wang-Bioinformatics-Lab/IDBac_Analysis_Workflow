import argparse
from pathlib import Path
import os
import logging
from psims.mzml.writer import MzMLWriter
from pyteomics import mzml
import sys
import subprocess
import ijson
import logging
from decimal import Decimal
from tqdm import tqdm
from joblib import Parallel, delayed
import numpy as np
from collections import defaultdict

import math
import json5
import json
from typing import List
import re
from multiprocessing import Pool
import pandas as pd

TEMP_MZML_DIR = Path('./temp_mzml')
INTERMEDIATE_MZML_DIR =Path('./temp_processed_mzml')
GROUP_SIZE = 100

#################### From ML side ####################
def convert_to_serializable(obj):
    """Convert Decimal types in the object to float."""
    if isinstance(obj, dict):
        return {key: convert_to_serializable(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_serializable(item) for item in obj]
    elif isinstance(obj, Decimal):
        return float(obj)  # Convert Decimal to float
    return obj

def process_line(line_output):
    line, output_mzML_dir = line_output
    if not line.strip():
        return
    obj = json.loads(line)
    obj = convert_to_serializable(obj)

    db_id = obj['database_id']
    all_scans = obj["spectrum"]

    with MzMLWriter(open(output_mzML_dir / f'{db_id}.mzML', 'wb'), close=True) as writer:
        writer.controlled_vocabularies()
        with writer.run(id='my_analysis'):
            with writer.spectrum_list(count=len(all_scans)):
                for scan_idx, scan in enumerate(all_scans):
                    mz_array, intensity_array = zip(*scan)
                    writer.write_spectrum(
                        mz_array,
                        intensity_array,
                        id=f'scan={scan_idx}',
                        params=[
                            "MS1 Spectrum",
                            {"ms level": 1},
                            {"total ion current": sum(intensity_array)}
                        ]
                    )

def write_mzML_files_from_json(json_input: Path, output_mzML_dir: Path, workers: int = 4) -> None:
    with open(json_input, 'r', encoding='utf-8') as input_file:
        lines = [(line, output_mzML_dir) for line in input_file if line.strip()]

    with Pool(processes=workers) as pool:
        list(tqdm(pool.imap_unordered(process_line, lines), total=len(lines), desc="Writing mzML files"))

def process_with_maldi_quant(script_path, input_path: Path, output_path: Path, n_jobs:int=-1):
    debug = logging.getLogger().getEffectiveLevel() == logging.DEBUG
    subprocess_output_path = subprocess.DEVNULL
    subprocess_check = False

    if debug:
        n_jobs = 1
        subprocess_output_path = sys.stdout
        subprocess_check = True

    
    def _run_rscript(spectrum_paths: List[Path],):
        for spectrum_path in spectrum_paths:
            _output_path = output_path / f"{spectrum_path.stem}.mzML"
            print(f"Processing {spectrum_path} to {_output_path}")
            subprocess.run(['Rscript', script_path, str(spectrum_path), str(_output_path)], 
                           stdout = subprocess_output_path, 
                           stderr = subprocess_output_path,
                           check  = subprocess_check)

    # Performs peak picking, baseline correction, binning, and merging
    all_spectrum_paths = list(input_path.glob("*.mzML"))

    # Split into groups to reduce overhead
    all_spectrum_paths = np.array_split(all_spectrum_paths, len(all_spectrum_paths) // GROUP_SIZE + 1)

    # Run the R script in parallel
    logging.info("Running MaldiQuant in Parallel")
    logging.info("Outputting files to %s", output_path)
    Parallel(n_jobs=n_jobs)(delayed(_run_rscript)(spectrum_paths) for spectrum_paths in tqdm(all_spectrum_paths, desc="Running MALDIquant"))

#################### Novel Stuff ####################
def convert_to_tensors(path: Path, output_dir: Path):
    all_spectrum_paths = list(path.glob("*.mzML"))
    logging.info(f"Found {len(all_spectrum_paths)} mzML files in {path}")
    

    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
    
    for spectrum_path in tqdm(all_spectrum_paths, desc="Converting to NPY"):
        output_path = output_dir/ f"{spectrum_path.stem}.npy"

        mzml_file = mzml.read(str(spectrum_path))
        # Average the scans
        num_scans = 0
        spectrum_dict = defaultdict(float)
        for scan in mzml_file:
            num_scans += 1
            for mz, intensity in zip(scan["m/z array"], scan["intensity array"]):
                spectrum_dict[mz] += intensity

        output_mz_array = []
        output_intensity_array = []
        for mz in spectrum_dict:
            output_mz_array.append(mz)
            output_intensity_array.append(spectrum_dict[mz] / num_scans)

        logging.info(f"Writing {spectrum_path.stem} to {output_path}")
        np.save(output_path, np.stack((np.array(output_mz_array), np.array(output_intensity_array))), allow_pickle=False)
def main():
    parser = argparse.ArgumentParser(description="Preprocess data for ML inference.")
    parser.add_argument("--input_json", help="Path to the input DB JSON file", required=False)
    parser.add_argument("--input_mzml", help="Path to the input mzML files", required=False)
    parser.add_argument("--output_dir", help="Directory to save the preprocessed data", required=True)
    parser.add_argument("--debug", action="store_true", help="Enable debug mode for additional logging")
    parser.add_argument("--rscript", help="Path to the R script for preprocessing", required=True)
    args = parser.parse_args()

    if not args.input_json and not args.input_mzml:
        raise ValueError("You must specify either --input_json or --input_mzml")
    if args.input_json and args.input_mzml:
        raise ValueError("You cannot specify both --input_json and --input_mzml. Please choose one.")

    logging.basicConfig(level=logging.INFO)

    for arg in vars(args):
        logging.info("%s: %s", arg, getattr(args, arg))

    if args.input_json:
        input_json = Path(args.input_json)
    else:
        input_json = None
    if args.input_mzml:
        input_mzml = Path(args.input_mzml)
    else:
        input_mzml = None
    output_dir = Path(args.output_dir)

    global TEMP_MZML_DIR

    TEMP_MZML_DIR.mkdir(parents=True, exist_ok=True)
    INTERMEDIATE_MZML_DIR.mkdir(parents=True, exist_ok=True)
    

    if input_json and not input_json.exists():
        raise ValueError(f"Input JSON file does not exist: {input_json}")
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    if args.input_json:
        write_mzML_files_from_json(input_json, TEMP_MZML_DIR)
    else:
        TEMP_MZML_DIR = input_mzml

    if args.debug:
        # Set logging level to DEBUG
        logging.getLogger().setLevel(logging.DEBUG)

    process_with_maldi_quant(args.rscript, TEMP_MZML_DIR, INTERMEDIATE_MZML_DIR, n_jobs=-1 if not args.debug else 1)

    # Write to npy
    convert_to_tensors(INTERMEDIATE_MZML_DIR, output_dir)



if __name__ == "__main__":
    main()