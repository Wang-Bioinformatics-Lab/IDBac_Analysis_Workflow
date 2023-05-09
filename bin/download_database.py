import sys
import os
import argparse
import pandas as pd
import requests
import json
from psims.mzml.writer import MzMLWriter


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('output_library_json')
    parser.add_argument('output_library_mzml')

    args = parser.parse_args()

    # Let's read the spectra from idbac kb
    url = "https://idbac-kb.gnps2.org/api/spectra"

    r = requests.get(url)
    all_spectra_list = r.json()

    # Now we'll iterate and get each individual spectrum
    output_spectra_list = []

    for spectrum in all_spectra_list:
        database_id = spectrum["database_id"]

        url = "https://idbac-kb.gnps2.org/api/spectrum"
        params = {}
        params["database_id"] = database_id

        r = requests.get(url, params=params)

        spectrum_dict = r.json()

        # add to json list
        output_spectra_list.append(spectrum_dict)

    # Saving mzML
    with MzMLWriter(open(args.output_library_mzml, 'wb'), close=True) as out:
        # Add default controlled vocabularies
        out.controlled_vocabularies()
        with out.run(id="my_analysis"):
            with out.spectrum_list(count=len(output_spectra_list)):
                scan = 1
                for spectrum in output_spectra_list:
                    # getting the mz peaks
                    peaks = spectrum["peaks"]

                    # unzip
                    mz_array = [x[0] for x in peaks]
                    intensity_array = [x[1] for x in peaks]

                    out.write_spectrum(
                        mz_array, intensity_array,
                        id="scan={}".format(scan), params=[
                            "MS1 Spectrum",
                            {"ms level": 1},
                            {"total ion current": sum(intensity_array)}
                        ])
                    
                    spectrum["scan"] = scan
                    
                    scan += 1

    # Saving JSON
    with open(args.output_library_json, "w") as output_file:
        json.dump(output_spectra_list, output_file)

    

if __name__ == '__main__':
    main()