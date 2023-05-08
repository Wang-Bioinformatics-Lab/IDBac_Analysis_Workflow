#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra_folder = ""

TOOL_FOLDER = "$baseDir/bin"

process processData {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file "input_spectra/*" 

    output:
    file 'output.html'

    """
    python $TOOL_FOLDER/process_data.py input_spectra output.html
    """
}

process baselineCorrection {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_maldiquant.yml"

    input:
    file input_file 

    output:
    file '*baselineCorrection.mzML'

    """
    Rscript $TOOL_FOLDER/baselineCorrection.R $input_file ${input_file}.baselineCorrection.mzML
    """
}

workflow {
    // Doing baseline correction
    input_mzml_files_ch = Channel.fromPath(params.input_spectra_folder + "/*.mzML")
    normalized_spectra_ch = baselineCorrection(input_mzml_files_ch)

    // Creating Spectra
    processData(normalized_spectra_ch.collect())
}