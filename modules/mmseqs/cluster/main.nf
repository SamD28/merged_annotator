process CLUSTER_SEQS {
    tag "clustering"
    label "clustering"
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    container 'quay.io/d_goryslavets/mmseqs2_pannotator:17-b804f-pannotator.1' 

    input:
    path(seqs_file)

    output:
    tuple path("${params.mmseqs_clusterPrefix}_all_seqs.fasta"), path("${params.mmseqs_clusterPrefix}_cluster.tsv"), path("${params.mmseqs_clusterPrefix}_rep_seq.fasta"), emit: cluster
    path("${params.mmseqs_tmpDir}/latest"), emit: db

    script:
    """
    mmseqs ${params.mmseqs_command} \\
    ${seqs_file} \\
    ${params.mmseqs_clusterPrefix} \\
    ${params.mmseqs_tmpDir} \\
    --min-seq-id ${params.mmseqs_min_id} \\
    -c ${params.mmseqs_min_cov} \\
    --cov-mode ${params.mmseqs_cov_mode} \\
    -a 1 \\
    --remove-tmp-files 0 \\
    --threads ${task.cpus}
    """
}

process ALIGNMENT_SCORE {
    tag "convertalis"

    container 'quay.io/d_goryslavets/mmseqs2_pannotator:17-b804f-pannotator.1'

    input:
    path(db_path)

    output:
    path("${params.mmseqs_clusterPrefix}_aln.tsv")

    script:
    """
    mmseqs align \\
        ${db_path}/input \\
        ${db_path}/input \\
        ${db_path}/clu \\
        ${db_path}/aln_proper \\
        -a \\
        --threads ${task.cpus}

    mmseqs convertalis \\
        ${db_path}/input \\
        ${db_path}/input \\
        ${db_path}/aln_proper \\
        ${params.mmseqs_clusterPrefix}_aln.tsv \\
        --format-mode 4 \\
        --format-output "query,target,qstart,qend,qlen,tstart,tend,tlen,pident,alnlen,mismatch,gapopen,qcov,tcov,evalue,bits,cigar,qseq,tseq,qaln,taln,qheader,theader"
    """
}

