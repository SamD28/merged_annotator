process GENOMAD {
    tag "${meta.ID}"

    container 'quay.io/biocontainers/genomad:1.7.0--pyhdfd78af_0'

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path(genomad_output)

    script:
    genomad_output = "${meta.ID}_result"
    """
    genomad end-to-end \\
        --cleanup \\
        ${faa} \\
        ${genomad_output} \\
        ${params.genomad_db}
    """
}