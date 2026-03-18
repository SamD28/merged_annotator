process KOFAMSCAN {
    container 'quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2'

    input:
    path(all_seqs)

    output:
    path("kofamscan_result.tsv")

    script:
    """
    exec_annotation \\
        -o kofamscan_result.tsv \\
        --cpu ${task.cpus} \\
        --e-value ${params.kofamscan_eval} \\
        --profile ${params.kofam_db}/profiles/prokaryote.hal \\
        --ko-list ${params.kofam_db}/ko_list \\
        -f mapper \\
        ${all_seqs}
    """
}