#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$baseDir/bin"

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

/* 
 * Simple sanity checks on mzML files that we can warn the users about
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

// Preprocess with MaldiQuant: Sqrt, Smoothing, Baseline Correction, Peak Detection
process baseline {
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

workflow data_prep {
    take:
    input_spectra_folder
    input_small_molecule_folder
    input_media_control_folder
    input_metadata_file

    main:

    // Define input channels based on parameters
    input_mzml_files_ch = Channel.fromPath(input_spectra_folder + "/*.mzML")
    pre_flight_ch = input_mzml_files_ch
    if (input_small_molecule_folder != "") {
        small_mol_ch = Channel.fromPath(input_small_molecule_folder + "/*.mzML")
        pre_flight_ch = pre_flight_ch.concat(small_mol_ch)
    } else {
        small_mol_ch = channel.empty()
    }
    if (input_media_control_folder != "") {
        blank_channel = Channel.fromPath(input_media_control_folder + "/*.mzML")
        pre_flight_ch = pre_flight_ch.concat(blank_channel)
    } else {
        blank_channel = channel.empty()
    }
    if (input_metadata_file != "NO_FILE") {
        metadata_file_ch = Channel.fromPath(input_metadata_file)
        metadataValidation(metadata_file_ch)
    } else {
        metadata_file_ch = channel.empty()
    }

    // Run Pre-flight Check and Error Output
    preFlightFailures = preFlightCheck(pre_flight_ch)
    outputErrors(preFlightFailures.collect())

    // Do a direct merge of protein spectra for downstream plotting
    mergeForPlotting(input_mzml_files_ch)

    // Doing baseline correction on main query spectra
    baseline_query_spectra_ch = baseline(input_mzml_files_ch)

    // Format the metadata for the main analysis
    formatted_metadata_ch = formatMetadata(metadata_file_ch)

    emit:
    input_mzml_files_ch = input_mzml_files_ch // for ML inference in core analysis
    baseline_query_spectra_ch = baseline_query_spectra_ch // main query spectra for merging/searching
    metadata_file_ch = metadata_file_ch // for small molecule and dendrogram
    formatted_metadata_ch = formatted_metadata_ch // for core analysis (if needed)
    small_mol_ch = small_mol_ch // for small molecule analysis (to avoid re-reading)
    blank_channel = blank_channel // for small molecule analysis (to avoid re-reading)
}