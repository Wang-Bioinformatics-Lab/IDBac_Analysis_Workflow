#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/bin"

process baseline {
    publishDir "./nf_output/small_molecule/baseline_corrected", mode: 'copy'

    cpus 2
    memory '20 GB'

    conda "$TOOL_FOLDER/conda_maldiquant.yml"

    errorStrategy 'ignore'

    input:
    file input_file 
    val small_mol_contains_ms2

    output:
    file 'baselinecorrected/*.mzML'

    """
    mkdir baselinecorrected
    Rscript $TOOL_FOLDER/baselineCorrectionSM.R $input_file "baselinecorrected/${input_file.fileName}" $small_mol_contains_ms2
    """
}

process summarize {
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

process baseline_blank {
    publishDir "./nf_output/small_molecule/baseline_corrected_blank", mode: 'copy'
    conda "$TOOL_FOLDER/conda_maldiquant.yml"

    cpus 2
    memory '20 GB'

    errorStrategy 'ignore'

    input:
    file input_file
    val small_mol_contains_ms2

    output:
    file 'baselinecorrected/*.mzML'
    """
    mkdir baselinecorrected
    Rscript $TOOL_FOLDER/baselineCorrectionSM.R $input_file baselinecorrected/$input_file $small_mol_contains_ms2
    """
}

process media_control {
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

process convert_to_mgf {
    publishDir "./nf_output/small_molecule/gnps_mgf", mode: 'copy'
    cache true
    conda "$TOOL_FOLDER/conda_env.yml"

    cpus 2
    memory '8 GB'

    input:
    file small_molecule_file

    output:
    file 'mgf_output/*.mgf'

    """
    mkdir -p mgf_output
    python $TOOL_FOLDER/convert_mzml_to_mgf.py \
        --input_file "${small_molecule_file}" \
        --output_file "mgf_output/${small_molecule_file.simpleName}.mgf"
    """
}

workflow small_mol {

    take:
    small_mol_ch                // Our strain-associated small molecule mzML files
    blank_channel               // Our strain-associated media control mzML files
    input_media_control_folder  // The directory containing those mzML files
    input_metadata_file         // Metadata file associated small molecule files with blank files
    small_mol_contains_ms2      // Whether small molecule spectra are Centroided MS1+MS2 or standard profile MS1 data

    main:
    baseline_corrected_small_molecule = baseline(small_mol_ch, small_mol_contains_ms2)

    // Let's do some basic sanity checking to provide users friendly errors
    if ((input_media_control_folder != "" && input_media_control_folder != "NO_FILE") && (input_metadata_file == "NO_FILE" || input_metadata_file == "")) {
         error "An input metadata file is required for media control"
    }
    // Requires newer nextflow version than prod:
    // if (input_media_control_folder != ""&& blank_channel.size() == 0)  {
    //      error "No blank files found in the provided media control folder: ${input_media_control_folder}"
    // }

    // Conditional Media Control Correction
    if (input_media_control_folder != "" && input_metadata_file != "NO_FILE" && input_metadata_file != "") {
        
        // Check for required metadata when media control is present
        if ((input_media_control_folder != "" && input_media_control_folder != "NO_FILE") && (input_metadata_file == "NO_FILE" || input_metadata_file == "")) {
             error "An input metadata file is required for media control"
        }

        baseline_corrected_blank = baseline_blank(blank_channel, small_mol_contains_ms2)
        
        // Metadata is guaranteed to be present here based on the if condition
        metadata_file_ch_small_mol = Channel.fromPath(input_metadata_file) 
        
        baseline_corrected_small_molecule = media_control(
            baseline_corrected_small_molecule, 
            metadata_file_ch_small_mol, 
            baseline_corrected_blank.collect()
        )
    }

    convert_to_mgf(baseline_corrected_small_molecule)

    summarize(baseline_corrected_small_molecule.collect())

    // Nothing to emit
}