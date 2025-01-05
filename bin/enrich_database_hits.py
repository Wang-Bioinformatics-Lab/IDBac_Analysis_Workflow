
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input_db_results')
    parser.add_argument('enriched_db_results')
    parser.add_argument('--database_json', help='Path to the database json file')

    args = parser.parse_args()

    # We might end up with a file without any entries
    try:
        results_df = pd.read_csv(args.input_db_results, sep="\t")
    except:
        print("Error reading input file, output empty file")
        open(args.enriched_db_results, "w").close()

        exit(0)

    results_list = results_df.to_dict(orient="records")

    database_df = pd.read_json(args.database_json)

    for result_obj in results_list:

        result_db_id = result_obj["database_id"]

        result_obj["db_taxonomy"] = ""
        try:
            this_db_entry = database_df[database_df["database_id"] == result_db_id]
            if not len(this_db_entry) == 1:
                raise Exception(f"Expected one entry for database_id {result_db_id}, got {len(this_db_entry)}")
            this_db_entry = this_db_entry.iloc[0]

            genus = this_db_entry["genus"]
            if genus == None:
                genus = this_db_entry["16S Taxonomy"]

            delimited_taxonomy = this_db_entry["family"] + ";" + genus + ";" + this_db_entry["species"]

            result_obj["db_taxonomy"] = delimited_taxonomy

        except Exception as e:
            print(e)

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