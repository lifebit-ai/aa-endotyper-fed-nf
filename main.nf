#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Process to read person info TSV file
process READ_PERSON_INFO {
    tag "Reading person info TSV"
    
    input:
    path person_tsv
    
    output:
    stdout
    
    script:
    """
    echo "Reading person info TSV file: ${person_tsv}"
    cat ${person_tsv}
    """
}

// Process to read measure info CSV file
process READ_MEASURE_INFO {
    tag "Reading measure info CSV"
    
    input:
    path measure_csv
    
    output:
    stdout
    
    script:
    """
    echo "Reading measure info CSV file: ${measure_csv}"
    cat ${measure_csv}
    """
}

// Main workflow
workflow {
    // Check if file parameters are provided
    if (!params.person_info_tsv) {
        error "Please provide person_info_tsv parameter"
    }
    
    if (!params.measure_info_csv) {
        error "Please provide measure_info_csv parameter"
    }
    
    // Create channels from the file paths
    person_tsv_ch = Channel.fromPath(params.person_info_tsv, checkIfExists: true)
    measure_csv_ch = Channel.fromPath(params.measure_info_csv, checkIfExists: true)
    
    // Run the processes
    READ_PERSON_INFO(person_tsv_ch)
    READ_MEASURE_INFO(measure_csv_ch)
    
    // Print the outputs
    READ_PERSON_INFO.out.view { "Person Info Content:\n$it" }
    READ_MEASURE_INFO.out.view { "Measure Info Content:\n$it" }
}