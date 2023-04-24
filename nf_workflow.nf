#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra_folder = ""

TOOL_FOLDER = "$baseDir/bin"

process processData {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_spectra 

    output:
    file 'output.html'

    """
    python $TOOL_FOLDER/process_data.py $input_spectra output.html
    """
}

workflow {
    data = Channel.fromPath(params.input_spectra_folder)
    processData(data)
}