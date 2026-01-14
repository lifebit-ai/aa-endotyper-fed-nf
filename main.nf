#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Process for local CNN training
process CNN_LOCAL {
    publishDir "results/local_cnn", mode: 'copy'
    tag "local_cnn_analysis"
    
    input:
    path sdy569_tidy
    path sdy569_cpeptide
    path sdy797_tidy
    path sdy797_cpeptide
    path sdy1737_tidy
    path sdy1737_cpeptide
    
    output:
    path "models/*"
    path "figures/**/*"
    
    script:
    """
    mkdir -p models figures/pdf figures/html test_data
    
    CNN_Local_noAutoencoders_3x3.py \
        --sdy569_feature ${sdy569_tidy} \
        --sdy569_cpeptide ${sdy569_cpeptide} \
        --sdy797_feature ${sdy797_tidy} \
        --sdy797_cpeptide ${sdy797_cpeptide} \
        --sdy1737_feature ${sdy1737_tidy} \
        --sdy1737_cpeptide ${sdy1737_cpeptide}
    """
}

// Process for federated CNN training
process CNN_FEDERATED {
    publishDir "results/federated_cnn", mode: 'copy'
    tag "federated_cnn_analysis"
    
    input:
    path sdy569_train
    path sdy569_test
    path sdy797_train
    path sdy797_test
    path sdy1737_train
    path sdy1737_test
    
    output:
    path "models/*"
    path "figures/**/*"
    
    script:
    """
    mkdir -p models figures/pdf figures/html
    
    CNN_Federated_noAutoencoder_3x3.py \
        --sdy569_train ${sdy569_train} \
        --sdy569_test ${sdy569_test} \
        --sdy797_train ${sdy797_train} \
        --sdy797_test ${sdy797_test} \
        --sdy1737_train ${sdy1737_train} \
        --sdy1737_test ${sdy1737_test}
    """
}

// Main workflow
workflow {
    // Create channels for tidy data files (for local CNN)
    sdy569_tidy_ch = Channel.fromPath("${params.data_dir}/SDY569_tidy.csv", checkIfExists: true)
    sdy569_cpeptide_ch = Channel.fromPath("${params.data_dir}/SDY569_cpeptide_auc_tidy.csv", checkIfExists: true)
    sdy797_tidy_ch = Channel.fromPath("${params.data_dir}/SDY797_tidy.csv", checkIfExists: true)
    sdy797_cpeptide_ch = Channel.fromPath("${params.data_dir}/SDY797_cpeptide_auc_tidy.csv", checkIfExists: true)
    sdy1737_tidy_ch = Channel.fromPath("${params.data_dir}/SDY1737_tidy.csv", checkIfExists: true)
    sdy1737_cpeptide_ch = Channel.fromPath("${params.data_dir}/SDY1737_cpeptide_auc_tidy.csv", checkIfExists: true)
    
    // Create channels for cleaned train/test files (for federated CNN)
    sdy569_train_ch = Channel.fromPath("${params.data_dir}/cleaned/SDY569_train.csv", checkIfExists: true)
    sdy569_test_ch = Channel.fromPath("${params.data_dir}/cleaned/SDY569_test.csv", checkIfExists: true)
    sdy797_train_ch = Channel.fromPath("${params.data_dir}/cleaned/SDY797_train.csv", checkIfExists: true)
    sdy797_test_ch = Channel.fromPath("${params.data_dir}/cleaned/SDY797_test.csv", checkIfExists: true)
    sdy1737_train_ch = Channel.fromPath("${params.data_dir}/cleaned/SDY1737_train.csv", checkIfExists: true)
    sdy1737_test_ch = Channel.fromPath("${params.data_dir}/cleaned/SDY1737_test.csv", checkIfExists: true)
    
    // Run local CNN analysis
    CNN_LOCAL(
        sdy569_tidy_ch,
        sdy569_cpeptide_ch,
        sdy797_tidy_ch,
        sdy797_cpeptide_ch,
        sdy1737_tidy_ch,
        sdy1737_cpeptide_ch
    )
    
    // Run federated CNN analysis
    CNN_FEDERATED(
        sdy569_train_ch,
        sdy569_test_ch,
        sdy797_train_ch,
        sdy797_test_ch,
        sdy1737_train_ch,
        sdy1737_test_ch
    )
}