#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Process to read person info TSV file
process RUN_ENDOTYPER {
    publishDir "results", mode: 'copy'
    tag "${person_tsv.simpleName}-${measure_csv.simpleName}"
    
    input:
    path person_tsv
    path measure_csv
    
    output:
    path "*.pdf"
    
    script:
    """
    Rscript ${projectDir}/bin/endotype.R ${person_tsv} ${measure_csv}
    """
}

// Main workflow
workflow {
    // Check if file parameters are provided
    if (!params.case_cohort) {
        error "Please provide case_cohort parameter"
    }
    
    if (!params.measure_info_csv) {
        error "Please provide measure_info_csv parameter"
    }
    
    // Create channels from the file paths
    person_tsv_ch = Channel.fromPath(params.case_cohort, checkIfExists: true)
    measure_csv_ch = Channel.fromPath(params.measure_info_csv, checkIfExists: true)
    
    // Run the processes
    RUN_ENDOTYPER(person_tsv_ch, measure_csv_ch)
}