#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra_folder = ""
params.input_small_molecule_folder = ""
params.input_media_control_folder = ""
params.input_metadata_file = "NO_FILE"

params.merge_replicates = "Yes"
params.distance = "cosine"
params.database_search_threshold = "0.7"
params.database_search_mass_range_lower = "2000"
params.database_search_mass_range_upper = "20000"
params.ml_search = "No" // If set to "Yes", it will use the ML database search, otherwise it will use the standard database search
params.metadata_column = "None"

// Extra peak processing parameters dependent on the MALDI instrument
params.MALDI_instrument = "general" // Default to general instrument settings
params.config = "https://idbac.org/analysis-utils/get_instrument_config" // URL or path to config

params.seed_genera = ""
params.seed_species = ""

params.heatmap_bin_size = 10  // Heatmaps use 1.0 Da bins by default
params.search_bin_size  = 10  // Database search uses 10.0 Da bins by default

params.debug_flag = "--debug" // Should be set to either "--debug" or ""

TOOL_FOLDER = "$baseDir/bin"

include { MLInferenceRawVectorsWorkflow as MLInferenceRawVectorsWorkflow } from "$baseDir/bin/ml_inference/ml_inference.nf"  addParams(output_dir: "./nf_output/ml_inference")

/* Simple sanity checks on mzML files that we can warn the users about
 * 1. Ensure each files has scans
 * 2. Ensure each scan has peaks with non-zero intensities (this tend to be a common issue with MALDI data)
 */
process preFlightCheck {
    conda "$TOOL_FOLDER/conda_env.yml"

    cpus 2
    memory '8 GB'

    input:
    each file(input_file )

    output:
    path('errors.csv'), optional: true

    """
    python3 $TOOL_FOLDER/test_inputs.py --input_file $input_file \
                                        --output_file errors.csv
    """
}

// Does some basic sanity checks on the metadata file
process metadataValidation {
    conda "$TOOL_FOLDER/conda_env.yml"

    cpus 2
    memory '8 GB'

    input:
    file metadata_file

    """
    python3 $TOOL_FOLDER/metadata_validation.py --metadata_table $metadata_file
    """
}

// This process outputs all of the errors from preFlightCheck to nf_output
process outputErrors {
    publishDir "./nf_output", mode: 'copy'

    cpus 2
    memory '8 GB'

    input:
    path input_files, stageAs: 'input_files/errors_*.csv'

    output:
    file 'errors.csv'

    """
    cat $input_files > errors.csv
    sed -i 1i"Error_Level,Scan,Filename,Error" errors.csv   # Add headers (since we are concatenating multiple files)
    """
}

process mergeForPlotting {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cpus 2
    memory '20 GB'
    cache false

    errorStrategy 'ignore'

    input:
    file input_file

    output:
    path 'raw_merged_for_plotting/*.mzML'

    """
    mkdir -p raw_merged_for_plotting
    python3 $TOOL_FOLDER/raw_merge.py \
        --input_file $input_file \
        --output_folder raw_merged_for_plotting \
    """
}

process baselineCorrection {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_maldiquant.yml"

    cpus 2
    memory '20 GB'

    errorStrategy 'ignore'

    input:
    file input_file 

    output:
    path 'baselinecorrected/*.mzML'

    """
    mkdir baselinecorrected
    Rscript $TOOL_FOLDER/baselineCorrection.R $input_file baselinecorrected/${input_file}
    """
}

// Optionally merges all spectra within the same file
// Note: This is only used for plotting, the outputs are not used for database search
process mergeInputSpectra {
    publishDir "./nf_output", mode: 'copy', pattern: "merged/*.mzML"
    publishDir "./nf_output", mode: 'copy', pattern: "merge_parameters.txt"
    publishDir "./nf_output", mode: 'copy', pattern: "bin_counts/*"

    cpus 2
    memory '16 GB'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file "input_spectra/*"

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
    --mass_range_lower ${params.database_search_mass_range_lower} \
    --mass_range_upper ${params.database_search_mass_range_upper} \
    --bin_size ${params.heatmap_bin_size} \
    --debug

    # Dump merge parameters to a file
    echo "mass_range_lower: ${params.database_search_mass_range_lower}" >> merge_parameters.txt
    echo "mass_range_upper: ${params.database_search_mass_range_upper}" >> merge_parameters.txt
    echo "bin_size: ${params.heatmap_bin_size}" >> merge_parameters.txt
    """
}

process baselineCorrectionSmallMolecule {
    publishDir "./nf_output/small_molecule/baseline_corrected", mode: 'copy'

    cpus 2
    memory '20 GB'

    conda "$TOOL_FOLDER/conda_maldiquant.yml"

    errorStrategy 'ignore'

    input:
    file input_file 

    output:
    file 'baselinecorrected/*.mzML'

    """
    mkdir baselinecorrected
    Rscript $TOOL_FOLDER/baselineCorrectionSM.R $input_file "baselinecorrected/${input_file.fileName}"
    """
}

process baselineCorrectionBlank {
    publishDir "./nf_output/small_molecule/baseline_corrected_blank", mode: 'copy'
    conda "$TOOL_FOLDER/conda_maldiquant.yml"

    cpus 2
    memory '20 GB'

    errorStrategy 'ignore'

    input:
    file input_file

    output:
    file 'baselinecorrected/*.mzML'
    """
    mkdir baselinecorrected
    Rscript $TOOL_FOLDER/baselineCorrectionSM.R $input_file baselinecorrected/$input_file
    """
}

process small_molecule_media_control {
    publishDir "./nf_output/small_molecule/media_control", mode: 'copy'
    cache true
    conda "$TOOL_FOLDER/conda_env.yml"

    cpus 2
    memory '8 GB'

    input:
    each small_molecule_file
    file metadata_file
    path media_control_dir, stageAs: "media_control_dir/*"

    output:
    file 'media_controled/*.mzML'

    """
    mkdir media_control
    python $TOOL_FOLDER/media_control.py \
        --small_molecule_file "${small_molecule_file}" \
        --metadata_file $metadata_file \
        --media_control_dir media_control_dir \
        --output_file "media_controled/${small_molecule_file.fileName}"
    """
}

process summarizeSmallMolecule {
    publishDir "./nf_output/small_molecule/", mode: 'copy'

    cache true

    cpus 2
    memory '8 GB'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file "input_spectra/*"

    output:
    file 'summary.json'

    """
    python $TOOL_FOLDER/summarize_small_molecule.py \
    --input_folder "input_spectra" \
    --output_file "summary.json"
    """

}

process formatMetadata {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cpus 2
    memory '8 GB'

    input:
    file metadata_file

    output:
    file 'output_histogram_data_directory/metadata.tsv'

    """
    python $TOOL_FOLDER/format_metadata.py \
    --output_histogram_data_directory "output_histogram_data_directory" \
    --metadata_path ${metadata_file}
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
        --merge_replicates ${params.merge_replicates} \
        --distance ${params.distance} \
        --bin_size ${params.search_bin_size} \
        --metadata_column "${params.metadata_column}" \
        --mass_range_lower ${params.database_search_mass_range_lower} \
        --mass_range_upper ${params.database_search_mass_range_upper} \
        $params.debug_flag
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


process downloadDatabase {
    publishDir "./nf_output/library", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cpus 2
    memory '8 GB'

    input:
    val x

    output:
    file 'idbac_database.json'

    """
    ml_flag=""
    if [ "${params.ml_search}" == "Yes" ]; then
        ml_flag="--ml"
    fi

    python $TOOL_FOLDER/download_database.py --output_library_json idbac_database.json \
                                             --download_bin_size ${params.search_bin_size} \${ml_flag} 
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

    output:
    file 'db_results.tsv'               // Distance between query spectra and database hits
    file 'db_db_distance.tsv'           // Distance between database hits
    file 'complete_output_results.tsv'  // Output results with metadata
    file 'query_query_distances.tsv'      // Distance between query spectra
    file 'query_spectra/*.mzML'

    """
    ml_flag=""
    if [ "${params.ml_search}" == "Yes" ]; then
        ml_flag="--ml"
    fi
    config_arg=""
    if [ "${params.config}" != "" ]; then   # TODO: Check if this is a URL or a file path
        wget --retry-connrefused --waitretry=1 --read-timeout=30 --timeout=30 -t 2 -O instrument_config.yaml ${params.config}
        config_arg="--config instrument_config.yaml"
    fi

    python $TOOL_FOLDER/database_search.py \
    --input_folder input_spectra \
    --database_filtered_json $idbac_database_filtered_json \
    --output_search_results_tsv db_results.tsv \
    --complete_output_results_tsv complete_output_results.tsv \
    --output_db_db_distance_tsv db_db_distance.tsv \
    --output_query_query_distances_tsv query_query_distances.tsv \
    --merge_replicates ${params.merge_replicates} \
    --score_threshold ${params.database_search_threshold} \
    --mass_range_lower ${params.database_search_mass_range_lower} \
    --mass_range_upper ${params.database_search_mass_range_upper} \
    --distance ${params.distance} \
    --bin_size ${params.search_bin_size} \
    --seed_genera "${params.seed_genera}" \
    --seed_species "${params.seed_species}" \
    $params.debug_flag \
    \${ml_flag} \
    --MALDI_instrument "${params.MALDI_instrument}" \
    \${config_arg}
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

workflow {
    input_mzml_files_ch = Channel.fromPath(params.input_spectra_folder + "/*.mzML")

    // Pre-flight check
    pre_flight_ch = input_mzml_files_ch
    if (params.input_small_molecule_folder != "") {
        small_mol_ch = Channel.fromPath(params.input_small_molecule_folder + "/*.mzML")
        pre_flight_ch = pre_flight_ch.concat(small_mol_ch)
    }
    if (params.input_media_control_folder != "") {
        blank_channel = Channel.fromPath(params.input_media_control_folder + "/*.mzML")
        pre_flight_ch = pre_flight_ch.concat(blank_channel)
    }
    if (params.input_metadata_file != "NO_FILE") {
        metadata_file_ch = Channel.fromPath(params.input_metadata_file)
        metadataValidation(metadata_file_ch)
    } else {
        metadata_file_ch = channel.empty()
    }

    preFlightFailures = preFlightCheck(pre_flight_ch)
    outputErrors(preFlightFailures.collect())

    // Summarizing small molecules
    if (params.input_small_molecule_folder != "") {
        baseline_corrected_small_molecule = baselineCorrectionSmallMolecule(small_mol_ch)
        if (params.input_media_control_folder != "" && params.input_metadata_file != "") {
            baseline_corrected_blank = baselineCorrectionBlank(blank_channel)
            baseline_corrected_small_molecule = small_molecule_media_control(baseline_corrected_small_molecule, Channel.fromPath(params.input_metadata_file), baseline_corrected_blank.collect())
        } else{
            if (params.input_media_control_folder != "") {
                error "An input metadata file is required for media control"
            }
        }
        summarizeSmallMolecule(baseline_corrected_small_molecule.collect())
    }

    // Do a direct merge of protein spectra for downstream plotting
    mergeForPlotting(input_mzml_files_ch)

    // Doing baseline correction
    // baseline_query_spectra_ch = Channel.empty()
    baseline_query_spectra_ch = baselineCorrection(input_mzml_files_ch)

    // Doing merging of spectra
    (merged_spectra_ch, merge_params, count_tables, replicate_counts) = mergeInputSpectra(baseline_query_spectra_ch.collect())

    // Summarizing Input
    summarizeSpectra(merged_spectra_ch.collect())

    // Downloading Database
    (output_idbac_database_ch) = downloadDatabase(1)            // In ML mode, this will contain ML data automatically

    // Select the data we'll use for the DB search
    processed_query_data = baseline_query_spectra_ch.collect() // Query data will need to be preprocessed and embedded
    if (params.ml_search == "Yes") {
        // If ML search is enabled, we run the ML inference workflow
        ml_inference_results_ch = MLInferenceRawVectorsWorkflow(input_mzml_files_ch.collect())
        processed_query_data = ml_inference_results_ch.collect()
    } 
    // Matching database to query spectra
    (core_search_results_ch, db_db_distances_ch, complete_search_results_ch, query_query_distance_ch, output_database_mzML) = databaseSearch(output_idbac_database_ch, processed_query_data)

    // Enriching database search results
    db_summary = downloadDatabaseSummary()
    enriched_core_results_db_ch     = enrichCoreDatabaseSearch(core_search_results_ch, db_summary)
    enriched_complete_results_db_ch = enrichCompleteDatabaseSearch(complete_search_results_ch, db_summary)

    // Format the metadata
    formatted_metadata_ch = formatMetadata(metadata_file_ch)

    // Creating Spectra Dendrogram
    if (params.input_metadata_file != "") {
        metadata_file_ch = Channel.fromPath(params.input_metadata_file)
    } else {
        metadata_file_ch = Channel.fromPath("None")
    }

    // If we have query spectra, we can create a dendrogram, otherwise we can't

    // collected_query_spectra = baseline_query_spectra_ch.collect()
    // baseline_query_spectra_ch.view()
    // baseline_query_spectra_ch = Channel.empty()
    
    // Logic to spoof outputs if no query spectra are present after merging
    collected = baseline_query_spectra_ch.collect()
    new_collected = collected.ifEmpty(file("No_Query_Spectra")) 
    enriched_core_results_db_ch = enriched_core_results_db_ch.collect()
    enriched_results = enriched_core_results_db_ch.ifEmpty(file("No_Enriched_Results"))

    // new_collected.view()
    // metadata_file_ch.view()
    // enriched_results.view()
    
    createDendrogram(new_collected, metadata_file_ch, enriched_results)


}