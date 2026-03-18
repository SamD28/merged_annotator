process VFDB {
    tag "${meta.ID}"
    
    container 'quay.io/biocontainers/diamond:2.1.23--hf93d47f_0'

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${meta.ID}_vfdb_result.txt"), emit: vfdb_result

    script:
    """
    diamond blastp -d ${params.vfdb_db} \\
        -q ${faa} \\
        -o ${meta.ID}_vfdb_result.txt \\
        -f 6 \\
        --max-target-seqs 1 \\
        --id ${params.vfdb_identity} \\
        --subject-cover ${params.vfdb_coverage} \\
        -e ${params.vfdb_e_value}
    """
}