process DBCAN {
    container 'quay.io/biocontainers/dbcan:5.2.8--pyhdfd78af_0'

    input:
    path(faa)

    output:
    path("dbcan/dbCAN_hmm_results.tsv")     , emit: results
    path("dbcan/overview.tsv")              , emit: overview
    path "versions.yml"                     , emit: versions

    script:
    """
    run_dbcan \\
        CAZyme_annotation \\
        --input_raw_data ${faa} \\
        --output_dir dbcan \\
        --threads ${task.cpus} \\
        --db_dir ${params.dbcan_db} \\
        --mode protein

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: 5.2.8
        dbcan database: 5.2.5
    END_VERSIONS
    """
}
