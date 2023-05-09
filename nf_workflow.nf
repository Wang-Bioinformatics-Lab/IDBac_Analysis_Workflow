#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra_folder = ""
params.input_metadata_file = ""

params.merge_replicates = "No"
params.similarity = "presence"

TOOL_FOLDER = "$baseDir/bin"


process baselineCorrection {
    publishDir "./nf_output", mode: 'copy'

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
    --similarity ${params.similarity}
    """
}

process downloadDatabase {
    publishDir "./nf_output/library", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    val x

    output:
    file 'idbac_database.json'
    file 'idbac_database.mzML'

    """
    python $TOOL_FOLDER/download_database.py \
    idbac_database.json \
    temp.mzML
    export LC_ALL=C && $TOOL_FOLDER/msconvert temp.mzML \
    --mzML --32 --outfile idbac_database.mzML
    """
}

process databaseSearch {
    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file idbac_database_mzML
    file idbac_database_json
    file query_mzML

    output:
    file 'search_results/*'

    """
    mkdir search_results
    python $TOOL_FOLDER/database_search.py \
    $query_mzML \
    $idbac_database_mzML \
    $idbac_database_json \
    search_results \
    --merge_replicates ${params.merge_replicates}
    """
}

workflow {
    // Doing baseline correction
    input_mzml_files_ch = Channel.fromPath(params.input_spectra_folder + "/*.mzML")
    normalized_spectra_ch = baselineCorrection(input_mzml_files_ch)

    // Creating Spectra Dendrogram
    metadata_file_ch = Channel.fromPath(params.input_metadata_file)
    processData(normalized_spectra_ch.collect(), metadata_file_ch)

    // Downloading Database
    (output_idbac_database_ch, output_idbac_mzML_ch) = downloadDatabase(1)

    // Baseline Correct the spectra
    baseline_corrected_mzML_ch = baselineCorrection2(output_idbac_mzML_ch)

    // Matching database to query spectra
    //search_results_ch = databaseSearch(baseline_corrected_mzML_ch, output_idbac_database_ch, normalized_spectra_ch)

}