#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Process for local CNN training on a single study
process CNN_LOCAL {
    publishDir "results/local_cnn/${study_id}", mode: 'copy'
    tag "local_cnn_${study_id}"
    
    input:
    val study_id
    path study_tidy_file
    path cpeptide_file
    
    output:
    path "models/*", emit: models
    path "figures/**/*", emit: figures
    path "*.csv", emit: results
    
    script:
    """
    mkdir -p models figures/pdf figures/html test_data
    
    # Create a generic local CNN script that works with single study
    python3 ${projectDir}/bin/CNN_Local_Single_Study.py \
        --study_id ${study_id} \
        --study_feature ${study_tidy_file} \
        --study_cpeptide ${cpeptide_file} \
        --output_dir ./
    """
}



// Main workflow - Generic single study local CNN training
workflow {
    // Input validation
    if (!params.study_tidy_file) {
        error "Please provide --study_tidy_file parameter"
    }
    if (!params.cpeptide_file) {
        error "Please provide --cpeptide_file parameter"
    }
    
    // Create channels from input parameters
    study_tidy_ch = Channel.fromPath(params.study_tidy_file, checkIfExists: true)
    cpeptide_ch = Channel.fromPath(params.cpeptide_file, checkIfExists: true)
    
    // Extract study_id from the tidy file name (e.g., SDY1625_tidy.csv -> SDY1625)
    study_id_ch = study_tidy_ch.map { file -> 
        def filename = file.getBaseName()
        def study_id = filename.replaceAll('_tidy$', '')
        return study_id
    }
    
    // Run local CNN analysis for the specified study
    CNN_LOCAL(
        study_id_ch,
        study_tidy_ch,
        cpeptide_ch
    )
}