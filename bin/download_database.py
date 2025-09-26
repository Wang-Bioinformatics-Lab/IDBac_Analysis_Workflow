import sys
import os
import argparse
import requests
import json
from tqdm import tqdm
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

@retry(
    stop=stop_after_attempt(3),  # Retry up to 3 times
    wait=wait_exponential(multiplier=1, min=2, max=10),  # Exponential backoff (2s, 4s, 8s)
    retry=retry_if_exception_type(requests.exceptions.RequestException),  # Retry on request errors
)
def get_with_retries(url, params=None):
    """Fetch data from the given URL with retries."""
    r = requests.get(url, params=params)
    if r.status_code == 404:
        print(f"404 error for {r.url}")
        return None

    r.raise_for_status()  # Raise exception for HTTP errors (4xx, 5xx)
    return r

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--output_library_json', required=True)
    parser.add_argument('--download_bin_size', type=int, required=False)
    parser.add_argument('--ml', required=False, action='store_true')
    args = parser.parse_args()

    if not args.ml and not args.download_bin_size:
        print("You must specify either --ml or --download_bin_size")
        sys.exit(1)


    # Let's read the spectra from idbac kb


    # Now we'll iterate and get each individual spectrum
    output_spectra_list = []
    if args.ml:
        
        url = "https://idbac.org/api/spectra"
        all_spectra_list = get_with_retries(url).json()
        db_id_to_taxa = {spectrum["database_id"]: {'genus': spectrum["genus"], 'species': spectrum["species"]} for spectrum in all_spectra_list}


        url = "https://idbac.org/api/spectrum/ml_db"

        ml_db = get_with_retries(url)
        # Iterate over lines, load with json.loads
        lines = ml_db.text.splitlines()
        for line in tqdm(lines, desc="Processing spectra", unit="line"):
            if line:
                ml_spectrum_dict = json.loads(line)

                params = {}
                params["database_id"] = ml_spectrum_dict["database_id"]
                params["bin_width"] = str(10)

              

                spectrum_dict = {}

                # Add the database_id, genus, and species
                spectrum_dict["database_id"] = ml_spectrum_dict["database_id"]
                # Remove 'peaks'
                spectrum_dict.pop("peaks", None)
                spectrum_dict["ml_embedding"] = ml_spectrum_dict["spectrum"]
                spectrum_dict["genus"] = db_id_to_taxa.get(ml_spectrum_dict["database_id"], {}).get("genus", "Unknown")
                spectrum_dict["species"] = db_id_to_taxa.get(ml_spectrum_dict["database_id"], {}).get("species", "Unknown")
                
                # Add to json list
                output_spectra_list.append(spectrum_dict)

    
    else:
        url = "https://idbac.org/api/spectra"

        all_spectra_list = get_with_retries(url).json()

        for spectrum in all_spectra_list:
            database_id = spectrum["database_id"]

            url = "https://idbac.org/api/spectrum/filtered"
            params = {}
            params["database_id"] = database_id
            params["bin_width"] = str(args.download_bin_size)

            r = get_with_retries(url, params=params)
            if r is None:
                print(f"Skipping {database_id}")
                continue

            spectrum_dict = r.json()
            spectrum_dict["database_id"] = database_id
            spectrum_dict["genus"] = spectrum["genus"]
            spectrum_dict["species"] = spectrum["species"]

            # add to json list
            output_spectra_list.append(spectrum_dict)

    # Saving JSON
    with open(args.output_library_json, "w") as output_file:
        json.dump(output_spectra_list, output_file, indent=4)

    

if __name__ == '__main__':
    main()