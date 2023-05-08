#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra_folder = ""
params.input_metadata_file = ""

params.merge_replicates = "No"

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
}