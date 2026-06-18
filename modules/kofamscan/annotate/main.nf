process KOFAMSCAN {
    tag "${meta.ID}"
    
    container 'quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2'

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("kofamscan_result.tsv")

    script:
    """
    exec_annotation \\
        -o kofamscan_result.tsv \\
        --cpu ${task.cpus} \\
        --e-value ${params.kofamscan_eval} \\
        --profile ${params.kofam_db}/profiles/prokaryote.hal \\
        --ko-list ${params.kofam_db}/ko_list \\
        -f mapper \\
        ${faa}
    """
}