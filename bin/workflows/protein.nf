#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/bin"

include { MLInferenceRawVectorsWorkflow as MLInferenceRawVectorsWorkflow } from "$baseDir/bin/ml_inference/ml_inference.nf"  addParams(output_dir: "./nf_output/ml_inference")

// Note: This is used for heatmaps, the outputs are not used for database search
process mergeInputSpectra {
    publishDir "./nf_output", mode: 'copy', pattern: "merged/*.mzML"
    publishDir "./nf_output", mode: 'copy', pattern: "merge_parameters.txt"
    publishDir "./nf_output", mode: 'copy', pattern: "bin_counts/*"

    cpus 2
    memory '16 GB'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file "input_spectra/*"
    val search_args

    output:
    file 'merged/*.mzML'
    file 'merge_parameters.txt'
    file 'bin_counts/bin_counts.csv' optional true   // Only outputs if it's merged
    file 'bin_counts/replicates.csv' optional true
    file 'bin_counts/binned_spectra.csv' optional true

    """
    mkdir merged
    python $TOOL_FOLDER/merge_spectra.py \
    input_spectra \
    merged \
    --mass_range_lower ${search_args.database_search_mass_range_lower} \
    --mass_range_upper ${search_args.database_search_mass_range_upper} \
    --bin_size ${search_args.heatmap_bin_size} \
    --debug

    # Dump merge parameters to a file
    echo "mass_range_lower: ${search_args.database_search_mass_range_lower}" >> merge_parameters.txt
    echo "mass_range_upper: ${search_args.database_search_mass_range_upper}" >> merge_parameters.txt
    echo "bin_size: ${search_args.heatmap_bin_size}" >> merge_parameters.txt
    """
}

process summarizeSpectra{
    publishDir "./nf_output/query", mode: 'copy'

    cpus 2
    memory '20 GB'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file "input_spectra/*"

    output:
    file 'summary.tsv'

    """
    python $TOOL_FOLDER/summarize_spectra.py \
    input_spectra \
    summary.tsv
    """
}

// Download spectra
process downloadDatabase {
    publishDir "./nf_output/library", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cpus 2
    memory '8 GB'

    input:
    val search_args

    output:
    file 'idbac_database.json'

    """
    ml_flag=""
    if [ "${search_args.ml_search_flag}" == "Yes" ]; then
        ml_flag="--ml"
    fi

    python $TOOL_FOLDER/download_database.py --output_library_json idbac_database.json \
                                             --download_bin_size ${search_args.search_bin_size} \${ml_flag}
    """
}

// Download database summary with retry up to 5 times
process downloadDatabaseSummary {
    
    cpus 2
    memory '8 GB'
    
    output:
    file 'idbac_database_summary.json'

    cache false

    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 5 https://idbac.org/api/spectra -O idbac_database_summary.json
    """
}

process databaseSearch {
    publishDir "./nf_output/search", mode: 'copy'

    cpus 2
    memory '20 GB'

    cache false

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file idbac_database_filtered_json
    file "input_spectra/*"
    val search_args

    output:
    file 'db_results.tsv'                       // Distance between query spectra and database hits
    path 'db_db_distance.tsv', optional: true   // Distance between database hits
    file 'complete_output_results.tsv'          // Output results with metadata
    file 'query_query_distances.tsv'            // Distance between query spectra
    file 'query_spectra/*.mzML'

    """
    ml_flag=""
    if [ "${search_args.ml_search_flag}" == "Yes" ]; then
        ml_flag="--ml"
    fi
    config_arg=""
    if [ "${search_args.config}" != "" ]; then   # TODO: Check if this is a URL or a file path
        wget --retry-connrefused --waitretry=1 --read-timeout=30 --timeout=30 -t 2 -O instrument_config.yaml ${search_args.config}
        config_arg="--config instrument_config.yaml"
    fi

    python $TOOL_FOLDER/database_search.py \
    --input_folder input_spectra \
    --database_filtered_json $idbac_database_filtered_json \
    --output_search_results_tsv db_results.tsv \
    --complete_output_results_tsv complete_output_results.tsv \
    --output_db_db_distance_tsv db_db_distance.tsv \
    --output_query_query_distances_tsv query_query_distances.tsv \
    --merge_replicates ${search_args.merge_replicates} \
    --score_threshold ${search_args.database_search_threshold} \
    --mass_range_lower ${search_args.database_search_mass_range_lower} \
    --mass_range_upper ${search_args.database_search_mass_range_upper} \
    --distance ${search_args.distance} \
    --unmatched_peak_penalty ${search_args.unmatched_peak_penalty} \
    --bin_size ${search_args.search_bin_size} \
    --seed_genera "${search_args.seed_genera}" \
    --seed_species "${search_args.seed_species}" \
    $search_args.debug_flag \
    \${ml_flag} \
    --MALDI_instrument "${search_args.MALDI_instrument}" \
    \${config_arg}

    # Manually check that db_db_distance.tsv exists, if distance is not reverse_*, since those do not produce db_db_distance.tsv
    if [ "${search_args.distance}" != "reverse_cosine" ] && [ "${search_args.distance}" != "reverse_presence" ] && [ ! -f db_db_distance.tsv ]; then
        echo "Error: db_db_distance.tsv was not created despite distance metric being set to '${search_args.distance}' which should produce this file." >&2
        exit 1
    fi
    """
}

process enrichCoreDatabaseSearch {
    publishDir "./nf_output/search", mode: 'copy'

    cpus 2
    memory '8 GB'

    conda "$TOOL_FOLDER/conda_env_enrichment.yml"

    cache false

    input:
    file input_results
    file idbac_database_summary
    
    output:
    file 'enriched_db_results.tsv'
    
    """
    python $TOOL_FOLDER/enrich_database_hits.py \
    $input_results \
    enriched_db_results.tsv \
    --database_json $idbac_database_summary
    """
}

process enrichCompleteDatabaseSearch {
    publishDir "./nf_output/search", mode: 'copy'

    cpus 2
    memory '8 GB'

    conda "$TOOL_FOLDER/conda_env_enrichment.yml"

    cache false

    input:
    file input_results
    file idbac_database_summary
    
    output:
    file 'complete_enriched_db_results.tsv'
    
    """
    python $TOOL_FOLDER/enrich_database_hits.py \
    $input_results \
    complete_enriched_db_results.tsv \
    --database_json $idbac_database_summary
    """
}

process createDendrogram {
    publishDir "./nf_output", mode: 'copy'

    cache true

    cpus 2
    memory '32 GB'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path input_file, stageAs: "input_spectra/*"
    file metadata_file
    file database_hits
    val  search_args

    output:
    file 'output.html' optional true
    file 'with_metadata.html' optional true
    file 'with_database.html' optional true
    file "output_distance_table.tsv" optional true
    file "output_histogram_data_directory/labels_spectra.tsv"  optional true
    file "output_histogram_data_directory/numerical_spectra.npy"  optional true

    script:
    
    if ("${input_file.fileName}" != "input_spectra/No_Query_Spectra")
        """
        mkdir output_histogram_data_directory

        python $TOOL_FOLDER/create_dendrogram.py \
        input_spectra \
        $metadata_file \
        $database_hits \
        output.html \
        with_metadata.html \
        with_database.html \
        output_distance_table.tsv \
        output_histogram_data_directory \
        --merge_replicates ${search_args.merge_replicates} \
        --distance ${search_args.distance} \
        --bin_size ${search_args.search_bin_size} \
        --metadata_column "${search_args.metadata_column}" \
        --mass_range_lower ${search_args.database_search_mass_range_lower} \
        --mass_range_upper ${search_args.database_search_mass_range_upper} \
        $search_args.debug_flag
        """
    else
        """
        mkdir -p output_histogram_data_directory
        # Single line with "filename" and a newline
        touch output_histogram_data_directory/labels_spectra.tsv
        echo "filename" > output_histogram_data_directory/labels_spectra.tsv
        # Empty numpy array
        python3 -c "import numpy as np; np.save('output_histogram_data_directory/numerical_spectra.npy', np.array([[]]))"
        """

}

workflow protein {

    take:
    baseline_query_spectra_ch
    input_mzml_files_ch
    input_metadata_file
    formatted_metadata_ch
    search_args

    main:
    // Merge already baseline and peak picked spectra
    (
        merged_spectra_ch,
        merge_params,
        count_tables,
        replicate_counts
    ) = mergeInputSpectra(
        baseline_query_spectra_ch.collect(),
        search_args
    )

    // Summarizing Input
    summarizeSpectra(
        merged_spectra_ch.collect()
    )

    // Downloading Database
    (
        output_idbac_database_ch
    ) = downloadDatabase(search_args)     // In ML mode, this will contain ML data automatically

    // Select the data we'll use for the DB search
    processed_query_data = baseline_query_spectra_ch.collect() // Query data will need to be preprocessed and embedded
    if (search_args.ml_search_flag == "Yes") {
        // If ML search is enabled, we run the ML inference workflow
        ml_inference_results_ch = MLInferenceRawVectorsWorkflow(
            input_mzml_files_ch.collect()
        )
        processed_query_data = ml_inference_results_ch.collect()
    }

    // Matching database to query spectra
    (
        core_search_results_ch,
        db_db_distances_ch,
        complete_search_results_ch,
        query_query_distance_ch,
        output_database_mzML
    ) = databaseSearch(
        output_idbac_database_ch,
        processed_query_data,
        search_args
    )

    // Enriching database search results
    db_summary = downloadDatabaseSummary()
    enriched_core_results_db_ch     = enrichCoreDatabaseSearch(
                                        core_search_results_ch,
                                        db_summary
                                    )
    enriched_complete_results_db_ch = enrichCompleteDatabaseSearch(
                                        complete_search_results_ch, db_summary)

    // Creating Spectra Dendrogram
    if (input_metadata_file != "") {
        metadata_file_ch = Channel.fromPath(input_metadata_file)
    } else {
        metadata_file_ch = Channel.fromPath("None")
    }
    
    // Logic to spoof outputs if no query spectra are present after merging
    collected                   = baseline_query_spectra_ch.collect()
    new_collected               = collected.ifEmpty(file("No_Query_Spectra")) 
    enriched_core_results_db_ch = enriched_core_results_db_ch.collect()
    enriched_results            = enriched_core_results_db_ch.ifEmpty(file("No_Enriched_Results"))

    createDendrogram(
        new_collected,
        metadata_file_ch,
        enriched_results,
        search_args
    )
}