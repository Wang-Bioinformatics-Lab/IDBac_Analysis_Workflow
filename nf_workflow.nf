#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_spectra_folder             = ""
params.input_small_molecule_folder      = ""
params.input_media_control_folder       = ""
params.input_metadata_file              = "NO_FILE"
params.publishdir                       = "."

params.merge_replicates                 = "Yes"     // Unused TODO: Remove
params.distance                         = "cosine"  // Options: "cosine", "euclidean", "presence", "reverse_cosine", "reverse_presence"
params.unmatched_peak_penalty           = 1.0       // Only used for "reverse_cosine" and "reverse_presence" distance metrics
params.database_search_threshold        = "0.7"
params.database_search_mass_range_lower = "3000"
params.database_search_mass_range_upper = "20000"
params.ml_search                        = "No"      // If set to "Yes", it will use the ML database search, otherwise it will use the standard database search
params.metadata_column                  = "None"

// Extra peak processing parameters dependent on the MALDI instrument
params.MALDI_instrument = "general" // Default to general instrument settings
params.config = "https://idbac.org/analysis-utils/get_instrument_config" // URL or path to config

params.seed_genera      = ""
params.seed_species     = ""

params.heatmap_bin_size = 10  // Heatmaps use 1.0 Da bins by default
params.search_bin_size  = 10  // Database search uses 10.0 Da bins by default

params.debug_flag = "" // Should be set to either "--debug" or ""

include { data_prep as data_prep }  from "$baseDir/bin/workflows/data_preparation.nf"
include { small_mol as small_mol }  from "$baseDir/bin/workflows/small_molecule.nf"
include { protein as protein     }  from "$baseDir/bin/workflows/protein.nf"

workflow {
    // GNPS2 displays this fixed file above the result links. Register the
    // handler before any validation so even graph-construction errors produce it.
    workflow.onError = {
        def marker = 'IDBAC_USER_ERROR:'
        def report = (workflow.errorReport ?: '').toString()
        def commandErrorStart = report.indexOf('Command error:')
        def commandErrorEnd = report.indexOf('\nWork dir:', commandErrorStart + 1)
        def commandError = commandErrorStart >= 0
            ? report.substring(commandErrorStart, commandErrorEnd >= 0 ? commandErrorEnd : report.size())
            : ''
        def markerSource = "${workflow.errorMessage ?: ''}\n${commandError}"
        def markedMessages = markerSource.readLines()
            .findAll { line -> line.contains(marker) }
            .collect { line -> line.substring(line.indexOf(marker) + marker.size()).trim() }
            .findAll { line -> line }
            .unique()

        def failureMessage = markedMessages
            ? markedMessages.join('\n')
            : (workflow.errorMessage ?: 'The workflow stopped because an unexpected error occurred.').toString()

        def summaryFile = file("${params.publishdir}/nf_output/workflow_summary.txt")
        summaryFile.parent.mkdirs()
        summaryFile.text = "ERROR: Workflow failed.\n\n${failureMessage}\n\nSee the task log for technical details.\n"
        println "Wrote user-facing failure summary to ${summaryFile}"
    }

    // ----------- General data preparation & sanity checks ----------- 
    data_prep(
        params.input_spectra_folder,
        params.input_small_molecule_folder,
        params.input_media_control_folder,
        params.input_metadata_file,
        params.database_search_mass_range_lower,
        params.database_search_mass_range_upper
    )

    input_mzml_files_ch         = data_prep.out.input_mzml_files_ch
    baseline_query_spectra_ch   = data_prep.out.baseline_query_spectra_ch
    small_mol_ch                = data_prep.out.small_mol_ch
    blank_channel               = data_prep.out.blank_channel
    metadata_file_ch            = data_prep.out.metadata_file_ch
    formatted_metadata_ch       = data_prep.out.formatted_metadata_ch

    // ----------- Check there is some valid input data -----------
    input_mzml_files_ch
    .mix( small_mol_ch )
    .ifEmpty { 
        error "IDBAC_USER_ERROR: No valid input spectra or small molecules were found. Ensure files are in .mzML format."
    }


    // ----------- Small Molecule Processing -----------
    small_mol(
        small_mol_ch,
        blank_channel,
        params.input_media_control_folder,
        params.input_metadata_file
    )

    // If ML Search, ensure cosine
    if (params.ml_search == "Yes" && params.distance != "cosine") {
        error "IDBAC_USER_ERROR: ML search requires the distance metric to be set to 'cosine'."
    }
    
    // ----------- Protein Analysis -----------
    def search_args = [
        "database_search_mass_range_lower": params.database_search_mass_range_lower,
        "database_search_mass_range_upper": params.database_search_mass_range_upper,
        "search_bin_size": params.search_bin_size,
        "database_search_threshold": params.database_search_threshold,
        "distance": params.distance,
        "unmatched_peak_penalty": params.unmatched_peak_penalty,
        "heatmap_bin_size": params.heatmap_bin_size,
        "MALDI_instrument": params.MALDI_instrument,
        "config": params.config,
        "merge_replicates": params.merge_replicates,
        "seed_genera": params.seed_genera,
        "seed_species": params.seed_species,
        "ml_search_flag": params.ml_search,
        "debug_flag": params.debug_flag
    ]

    protein(
        baseline_query_spectra_ch,
        input_mzml_files_ch,
        params.input_metadata_file,
        search_args,
    )

}
