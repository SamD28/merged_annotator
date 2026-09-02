process GENOMAD {
    tag "${meta.ID}"

    container 'quay.io/biocontainers/genomad:1.12.0--pyhdfd78af_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(genomad_output)

    script:
    genomad_output = "${meta.ID}_result"
    """
    genomad end-to-end \\
        --cleanup \\
        ${fasta} \\
        ${genomad_output} \\
        ${params.genomad_db}
    """
}