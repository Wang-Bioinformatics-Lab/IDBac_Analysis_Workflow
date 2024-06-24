#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra_folder = ""
params.input_small_molecule_folder = ""
params.input_media_control_folder = ""
params.input_metadata_file = ""

params.merge_replicates = "No"
params.distance = "presence"
params.database_search_threshold = "0.7"
params.database_search_mass_range_lower = "2000"
params.database_search_mass_range_upper = "20000"
params.metadata_column = "None"

TOOL_FOLDER = "$baseDir/bin"

/* Simple sanity checks on mzML files that we can warn the users about
 * 1. Ensure each files has scans
 * 2. Ensure each scan has peaks with non-zero intensities (this tend to be a common issue with MALDI data)
 */
process preFlightCheck {
    conda "$TOOL_FOLDER/conda_env.yml"

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

    input:
    file metadata_file

    output:
    path 'errors.csv'

    """
    python3 $TOOL_FOLDER/metadata_validation.py --metadata_table $metadata_file
    """
}

// This process outputs all of the errors from preFlightCheck to nf_output
process outputErrors {
    publishDir "./nf_output", mode: 'copy'

    input:
    path input_files, stageAs: 'input_files/errors_*.csv'

    output:
    file 'errors.csv'

    """
    cat $input_files > errors.csv
    sed -i 1i"Filename,Error" errors.csv   # Add headers (since we are concatenating multiple files)
    """
}

process baselineCorrection {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_maldiquant.yml"

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
    --merge_replicates ${params.merge_replicates} \
    --mass_range_lower ${params.database_search_mass_range_lower} \
    --mass_range_upper ${params.database_search_mass_range_upper}
    """
}

process baselineCorrectionSmallMolecule {
    publishDir "./nf_output/small_molecule/baseline_corrected", mode: 'copy'

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

process baselineCorrectionBlank {
    publishDir "./nf_output/small_molecule/baseline_corrected_blank", mode: 'copy'
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

process small_molecule_media_control {
    publishDir "./nf_output/small_molecule/media_control", mode: 'copy'
    cache false
    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    each small_molecule_file
    file metadata_file
    path media_control_dir, stageAs: "media_control_dir/*"

    output:
    file 'media_controled/*.mzML'

    """
    mkdir media_control
    python $TOOL_FOLDER/media_control.py \
        --small_molecule_file $small_molecule_file \
        --metadata_file $metadata_file \
        --media_control_dir media_control_dir \
        --output_file media_controled/${small_molecule_file.fileName}
    """
}

process summarizeSmallMolecule {
    publishDir "./nf_output/small_molecule/", mode: 'copy'

    cache false

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
    file "output_distance_table.tsv" optional true
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
    output_distance_table.tsv \
    output_histogram_data_directory \
    --merge_replicates ${params.merge_replicates} \
    --distance ${params.distance} \
    --metadata_column "${params.metadata_column}" \
    --mass_range_lower ${params.database_search_mass_range_lower} \
    --mass_range_upper ${params.database_search_mass_range_upper}
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
    file 'db_db_distance.tsv'
    file 'complete_output_results.tsv'

    """
    python $TOOL_FOLDER/database_search.py \
    input_spectra \
    $idbac_database_filtered_json \
    db_results.tsv \
    complete_output_results.tsv \
    db_db_distance.tsv \
    --merge_replicates ${params.merge_replicates} \
    --score_threshold ${params.database_search_threshold} \
    --mass_range_lower ${params.database_search_mass_range_lower} \
    --mass_range_upper ${params.database_search_mass_range_upper}
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
    if (params.input_metadata_file != "") {
        metadata_file_ch = Channel.fromPath(params.input_metadata_file)
        metadataValidation(metadata_file_ch)
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

    // Doing baseline correction
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