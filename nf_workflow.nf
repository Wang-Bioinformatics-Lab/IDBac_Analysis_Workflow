#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra_folder = ""
params.input_metadata_file = ""

params.merge_replicates = "No"
params.similarity = "presence"
params.database_search_threshold = "0.7"

params.metadata_column = "None"

TOOL_FOLDER = "$baseDir/bin"


process baselineCorrection {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_maldiquant.yml"

    errorStrategy 'ignore'

    input:
    file input_file 

    output:
    file 'baselinecorrected/*.mzML'

    """
    mkdir baselinecorrected
    Rscript $TOOL_FOLDER/baselineCorrection.R $input_file baselinecorrected/${input_file}
    """
}

process mergeInputSpectra {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file "input_spectra/*"

    output:
    file 'merged/*.mzML'

    """
    mkdir merged
    python $TOOL_FOLDER/merge_spectra.py \
    input_spectra \
    merged \
    --merge_replicates ${params.merge_replicates}
    """
}

process createDendrogram {
    publishDir "./nf_output", mode: 'copy'

    cache false

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file "input_spectra/*"
    file metadata_file
    file database_hits

    output:
    file 'output.html' optional true
    file 'with_metadata.html' optional true
    file 'with_database.html' optional true
    file "output_similarity_table.tsv" optional true
    file "output_histogram_data_directory"  optional true

    """
    mkdir output_histogram_data_directory

    python $TOOL_FOLDER/create_dendrogram.py \
    input_spectra \
    $metadata_file \
    $database_hits \
    output.html \
    with_metadata.html \
    with_database.html \
    output_similarity_table.tsv \
    output_histogram_data_directory \
    --merge_replicates ${params.merge_replicates} \
    --similarity ${params.similarity} \
    --metadata_column "${params.metadata_column}"
    """
}

process downloadDatabase {
    publishDir "./nf_output/library", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'idbac_database.json'

    """
    python $TOOL_FOLDER/download_database.py \
    idbac_database.json
    """
}

process databaseSearch {
    publishDir "./nf_output/search", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file idbac_database_filtered_json
    file "input_spectra/*"

    output:
    file 'db_results.tsv'
    file 'db_db_similarity.tsv'

    """
    python $TOOL_FOLDER/database_search.py \
    input_spectra \
    $idbac_database_filtered_json \
    db_results.tsv \
    db_db_similarity.tsv \
    --merge_replicates ${params.merge_replicates} \
    --score_threshold ${params.database_search_threshold}
    """
}

process enrichDatabaseSearch {
    publishDir "./nf_output/search", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env_enrichment.yml"

    input:
    file input_results
    
    output:
    file 'enriched_db_results.tsv'
    
    """
    python $TOOL_FOLDER/enrich_database_hits.py \
    $input_results \
    enriched_db_results.tsv
    """
}

process summarizeSpectra{
    publishDir "./nf_output/query", mode: 'copy'

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
    // Doing baseline correction
    input_mzml_files_ch = Channel.fromPath(params.input_spectra_folder + "/*.mzML")
    baseline_query_spectra_ch = baselineCorrection(input_mzml_files_ch)

    // Doing merging of spectra
    merged_spectra_ch = mergeInputSpectra(baseline_query_spectra_ch.collect())

    // Summarizing Input
    summarizeSpectra(merged_spectra_ch.collect())

    // Downloading Database
    (output_idbac_database_ch) = downloadDatabase(1)

    // Matching database to query spectra
    (search_results_ch, output_database_mzML) = databaseSearch(output_idbac_database_ch, baseline_query_spectra_ch.collect())

    // Enriching database search results
    enriched_results_db_ch = enrichDatabaseSearch(search_results_ch)

    // Creating Spectra Dendrogram
    metadata_file_ch = Channel.fromPath(params.input_metadata_file)
    createDendrogram(baseline_query_spectra_ch.collect(), metadata_file_ch, enriched_results_db_ch)

}