import sys
import os
import argparse
import pandas as pd
import requests
import requests_cache
import json
import glob
import xmltodict

requests_cache.install_cache('idbac_cache')

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input_db_results')
    parser.add_argument('enriched_db_results')    

    args = parser.parse_args()

    # We might end up with a file without any entries
    try:
        results_df = pd.read_csv(args.input_db_results, sep="\t")
    except:
        print("Error reading input file, output empty file")
        open(args.enriched_db_results, "w").close()

        exit(0)

    results_list = results_df.to_dict(orient="records")

    for result_obj in results_list:
        print(result_obj)

        try:
            genbank_accession = result_obj["db_genbank_accession"]

            # Updating the URL
            mapping_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nucleotide&db=nucleotide&id={}&rettype=gb&retmode=xml".format(genbank_accession)

            r = requests.get(mapping_url)
            result_dictionary = xmltodict.parse(r.text)
            print(result_dictionary)
            try:
                nuccore_id = result_dictionary["eLinkResult"]["LinkSet"][0]["IdList"]["Id"]
            except:
                nuccore_id = result_dictionary["eLinkResult"]["LinkSet"]["IdList"]["Id"]

            # here we will use an API to get the information
            xml_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&retmode=xml".format(nuccore_id)

            r = requests.get(xml_url)
            result_dictionary = xmltodict.parse(r.text)

            # Getting taxonomy
            taxonomy = result_dictionary["GBSet"]["GBSeq"]["GBSeq_taxonomy"]
            organism = result_dictionary["GBSet"]["GBSeq"]["GBSeq_organism"]

            result_obj["db_taxonomy"] = taxonomy + "; " + organism
        except KeyboardInterrupt:
            raise
        except:
            print("Error Parsing")

    # outputting
    if len(results_list) == 0:  # Spoof a result so we show an empty table
        columns = ['query_filename', 'query_index', 'database_id', 'database_scan', 
                   'db_strain_name', 'db_culture_collection', 'db_taxonomy', 'distance']
        output_df = pd.DataFrame(columns=columns)
    else:
        output_df = pd.DataFrame(results_list)
    output_df.to_csv(args.enriched_db_results, sep="\t", index=False)



if __name__ == '__main__':
    main()