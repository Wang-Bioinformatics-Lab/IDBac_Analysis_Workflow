#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.output_dir = "./nf_output"

TOOL_FOLDER = "$moduleDir/bin"

process ml_preprocessing {
    publishDir "${params.output_dir}", mode: 'copy'

    conda "$TOOL_FOLDER/../../conda_maldiquant.yml"

    input:
    file input_json // Single json for the entire DB

    output:
    file 'preprocessed_data/*'

    """
    python $TOOL_FOLDER/preprocess.py \
            --input_json ${input_json} \
            --rscript $TOOL_FOLDER/preprocess_data.R \
            --output_dir ./preprocessed_data/
    """
}

process ml_preprocessing_from_mzml {
    publishDir "${params.output_dir}", mode: 'copy'

    conda "$TOOL_FOLDER/../../conda_maldiquant.yml"

    input:
    path input_mzml, stageAs: 'input_mzml_files/*'

    output:
    file 'preprocessed_data/*'

    """
    python $TOOL_FOLDER/preprocess.py \
            --input_mzml input_mzml_files/ \
            --rscript $TOOL_FOLDER/preprocess_data.R \
            --output_dir ./preprocessed_data/
    """
}

process ml_inference  {
    publishDir "${params.output_dir}", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path preprocessed_spectra, stageAs: 'ml_data_directory/*'

    output:
    path 'ml_vectors.feather'

    """
    python $TOOL_FOLDER/ml_inference.py \
            --ml_data_directory ml_data_directory \
            --output_file ml_vectors.feather \
            --onnx_model $TOOL_FOLDER/models/CLIP_Transformer.onnx
    """
}

process generate_json {
    publishDir "${params.output_dir}", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file idbac_full_spectra_json // Original database with metadata
    file idbac_ml_vectors        // ML vectors from inference
    
    output:
    file 'idbac_ml_db.json'

    """
    python $TOOL_FOLDER/generate_json.py \
            --idbac_full_spectra_json $idbac_full_spectra_json \
            --idbac_ml_vectors $idbac_ml_vectors \
            --output_file 'idbac_ml_db.json'
    """
}

workflow MLInferenceWorkflow {
    take:
    idbac_full_spectra_json

    main:
    preprocessed_data = ml_preprocessing(idbac_full_spectra_json)
    ml_vectors = ml_inference(preprocessed_data.collect())
    idbac_ml_db_json = generate_json(idbac_full_spectra_json, ml_vectors)

}

workflow MLInferenceRawVectorsWorkflow {
    take:
    input_mzml_files_ch

    main:
    preprocessed_data = ml_preprocessing_from_mzml(input_mzml_files_ch.collect())
    ml_vectors = ml_inference(preprocessed_data.collect())
    
    // Output the raw vectors without generating a JSON
    emit:
    ml_vectors
}