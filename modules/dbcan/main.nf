process DBCAN {
    container 'quay.io/biocontainers/dbcan:5.2.8--pyhdfd78af_0'

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("dbcan/dbCAN_hmm_results.tsv")     , emit: results
    tuple val(meta), path("dbcan/overview.tsv")              , emit: overview
    path "versions.yml"                     , emit: versions

    script:
    def eval_override = params.cazyme_hmm_eval ? "--hmmevalue ${params.cazyme_hmm_eval}" : ''
    def cov_override  = params.cazyme_hmm_cov  ? "--hmmcov ${params.cazyme_hmm_cov}"  : ''
    """
    run_dbcan \\
        CAZyme_annotation \\
        --input_raw_data ${faa} \\
        --output_dir dbcan \\
        --threads ${task.cpus} \\
        --db_dir ${params.dbcan_db} \\
        --mode protein \\
        ${eval_override} \\
        ${cov_override}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: 5.2.8
        dbcan database: 5.2.5
    END_VERSIONS
    """
}

process RUNDBCAN_EASYSUBSTRATE {
    tag "$meta.ID"

    label 'process_medium'

    container 'quay.io/biocontainers/dbcan:5.2.8--pyhdfd78af_0'

    input:
    tuple val(meta), path(faa), path(gff)

    output:
    tuple val(meta), path("${meta.ID}_substrate_prediction.tsv")     , emit: dbcan_substrates
    tuple val(meta), path("${meta.ID}_cgc_standard_out*")            , emit: dbcan_predictions
    tuple val(meta), path("${meta.ID}_cgc_standard_out_summary.tsv") , emit: dbcan_predictions_summary
    path  "versions.yml"                                             , emit: versions

    script:
    def eval_override = params.cazyme_hmm_eval ? "--hmmevalue ${params.cazyme_hmm_eval}" : ''
    def cov_override  = params.cazyme_hmm_cov  ? "--hmmcov ${params.cazyme_hmm_cov}"  : ''
    """
    run_dbcan easy_substrate \\
        --mode protein \\
        --db_dir ${params.dbcan_db} \\
        --input_raw_data ${faa} \\
        --output_dir ${meta.ID}_dbcan \\
        --input_gff ${gff} \\
        ${eval_override} \\
        ${cov_override}

    # cp rather than move so it caches
    cp ${meta.ID}_dbcan/substrate_prediction.tsv ${meta.ID}_substrate_prediction.tsv
    cp ${meta.ID}_dbcan/cgc_standard_out.tsv ${meta.ID}_cgc_standard_out.tsv
    cp ${meta.ID}_dbcan/cgc_standard_out_summary.tsv ${meta.ID}_cgc_standard_out_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dbcan: \$(echo \$(run_dbcan version) | cut -f2 -d':' | cut -f2 -d' ')
    END_VERSIONS
    """
}