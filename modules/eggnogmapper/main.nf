process EGGNOGMAPPER {
    tag "${meta.ID}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2':
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.emapper.annotations")   , emit: annotations
    tuple val(meta), path("*.emapper.seed_orthologs"), emit: orthologs, optional: true
    tuple val(meta), path("*.emapper.hits")          , emit: hits     , optional: true
    tuple val("${task.process}"), val('eggnog-mapper'), eval("emapper.py --version 2>&1 | grep -o 'emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+' | sed 's/emapper-//'"), topic: versions, emit: versions_eggnogmapper

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix        = task.ext.prefix ?: "${meta.ID}"
    def args          = task.ext.args   ?: ''
    def is_compressed = fasta.extension == '.gz'
    def fasta_name    = is_compressed ? fasta.baseName : "$fasta"
    def search_mode   = params.eggnog_search_mode
    def db            = params.eggnog_data_dir
    def db_flags      = ['diamond': '--dmnd_db', 'novel_fams': '--dmnd_db', 'mmseqs': '--mmseqs_db', 'hmmer': '--database', 'no_search': '--annotate_hits_table', 'cache': '--cache']
    def db_path       = (db instanceof Path && db.isDirectory()) ? "${db}/${db.name}"                    : "$db"
    def db_arg        = db && db_flags[search_mode]              ? "${db_flags[search_mode]} ${db_path}" : ''
    def dbmem         = task.memory.toMega() > 40000             ? '--dbmem'                             : ''
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    emapper.py \\
        $args \\
        --cpu ${task.cpus} \\
        -i ${fasta_name} \\
        --data_dir ${params.eggnog_data_dir} \\
        -m ${params.eggnog_search_mode} \\
        --output ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.ID}"
    """
    echo $args

    touch ${prefix}.emapper.annotations
    touch ${prefix}.emapper.seed_orthologs
    touch ${prefix}.emapper.hits
    """
}