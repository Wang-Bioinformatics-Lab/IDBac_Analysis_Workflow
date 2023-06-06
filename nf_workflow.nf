#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra_folder = ""
params.input_metadata_file = ""

params.merge_replicates = "No"
params.similarity = "presence"

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

process baselineCorrection2 {
    publishDir "./nf_output/library", mode: 'copy'

    conda "$TOOL_FOLDER/conda_maldiquant.yml"

    input:
    file input_file 

    output:
    file 'baselinecorrected/*.mzML'

    """
    mkdir baselinecorrected
    Rscript $TOOL_FOLDER/baselineCorrection.R $input_file baselinecorrected/${input_file}
    """
}

process processData {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false

    input:
    file "input_spectra/*"
    file metadata_file

    output:
    file 'output.html'
    file 'with_metadata.html'

    """
    python $TOOL_FOLDER/process_data.py input_spectra $metadata_file \
    output.html \
    with_metadata.html \
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
    file 'idbac_database_scanmapping.tsv'
    file 'idbac_database.mzML'

    """
    python $TOOL_FOLDER/download_database.py \
    idbac_database.json \
    idbac_database_scanmapping.tsv \
    temp.mzML

    export LC_ALL=C && $TOOL_FOLDER/msconvert temp.mzML \
    --mzML --32 --outfile idbac_database.mzML
    """
}

process databaseSearch {
    publishDir "./nf_output/search", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file idbac_database_mzML
    file idbac_database_scan_mapping
    file "input_spectra/*"

    output:
    file 'db_results.tsv'

    """
    python $TOOL_FOLDER/database_search.py \
    input_spectra \
    $idbac_database_mzML \
    $idbac_database_scan_mapping \
    db_results.tsv \
    --merge_replicates ${params.merge_replicates}
    """
}

process summarizeBaselineResolvedSpectra{
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

    // Creating Spectra Dendrogram
    metadata_file_ch = Channel.fromPath(params.input_metadata_file)
    processData(baseline_query_spectra_ch.collect(), metadata_file_ch)

    // Summarizing Input
    summarizeBaselineResolvedSpectra(baseline_query_spectra_ch.collect())

    // Downloading Database
    (output_idbac_database_ch, output_scan_mapping_ch, output_idbac_mzML_ch) = downloadDatabase(1)

    // Baseline Correct the spectra
    baseline_corrected_mzML_ch = baselineCorrection2(output_idbac_mzML_ch)

    // Matching database to query spectra
    search_results_ch = databaseSearch(baseline_corrected_mzML_ch, output_scan_mapping_ch, baseline_query_spectra_ch.collect())

}