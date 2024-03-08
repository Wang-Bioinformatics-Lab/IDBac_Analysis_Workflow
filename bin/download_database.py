import sys
import os
import argparse
import requests
import json


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('output_library_json')

    args = parser.parse_args()

    # Let's read the spectra from idbac kb
    url = "https://idbac.org/api/spectra"

    r = requests.get(url)
    all_spectra_list = r.json()


    # Now we'll iterate and get each individual spectrum
    output_spectra_list = []

    for spectrum in all_spectra_list:
        database_id = spectrum["database_id"]

        # TODO: Download everything all at once

        url = "https://idbac.org/api/spectrum/filtered"
        params = {}
        params["database_id"] = database_id

        try:
            r = requests.get(url, params=params)

            # raise for status
            r.raise_for_status()

            spectrum_dict = r.json()
            spectrum_dict["database_id"] = database_id
        except:
            continue

        # add to json list
        output_spectra_list.append(spectrum_dict)

    # Saving JSON
    with open(args.output_library_json, "w") as output_file:
        json.dump(output_spectra_list, output_file, indent=4)

    

if __name__ == '__main__':
    main()