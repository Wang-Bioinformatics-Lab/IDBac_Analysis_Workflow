#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/bin"

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

workflow small_mol {

    take:
    small_mol_ch                // Our strain-associated small molecule mzML files
    blank_channel               // Our strain-associated media control mzML files
    input_media_control_folder  // The directory containing those mzML files
    input_metadata_file         // Metadata file associated small molecule files with blank files

    main:
    baseline_corrected_small_molecule = baselineCorrectionSmallMolecule(small_mol_ch)

    // Conditional Media Control Correction
    if (input_media_control_folder != "" && input_metadata_file != "NO_FILE" && input_metadata_file != "") {
        
        // Check for required metadata when media control is present
        if (input_media_control_folder != "" && (input_metadata_file == "NO_FILE" || input_metadata_file == "")) {
             error "An input metadata file is required for media control"
        }

        baseline_corrected_blank = baselineCorrectionBlank(blank_channel)
        
        // Metadata is guaranteed to be present here based on the if condition
        metadata_file_ch_small_mol = Channel.fromPath(input_metadata_file) 
        
        baseline_corrected_small_molecule = small_molecule_media_control(
            baseline_corrected_small_molecule, 
            metadata_file_ch_small_mol, 
            baseline_corrected_blank.collect()
        )
    }

    summarizeSmallMolecule(baseline_corrected_small_molecule.collect())

    // Nothing to emit
}