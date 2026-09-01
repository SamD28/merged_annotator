include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { MANIFEST_PARSE   } from './subworkflows/manifest_parse.nf'
include { ANNOTATION       } from './subworkflows/annotation.nf'
include { ARG              } from './subworkflows/arg.nf'
include { FUNC             } from './subworkflows/func.nf'
include { MGE              } from './subworkflows/mge.nf'

workflow {
    main:
    log.info paramsSummaryLog(workflow)

    validateParameters()

    if (params.cluster_proteome)   {
        error "Clustering of proteomes is currently not supported and under development. Please set --cluster_proteome to false."
    }
    
    ch_input_for_annotation =  MANIFEST_PARSE(params.manifest)
    
    ANNOTATION(ch_input_for_annotation)

    ch_new_annotation = ch_input_for_annotation
            .join(ANNOTATION.out.faa)
            .join(ANNOTATION.out.gff)
            .multiMap { meta, fasta, faa, gff ->
                fastas: [meta, fasta]
                faas: [meta, faa]
                gffs: [meta, gff]
            }

    if (params.arg_annotate) {
        ARG(
            ch_new_annotation.fastas,
            ch_new_annotation.faas
        )
    }

    if (params.func_annotate) {
        FUNC(
            ch_new_annotation.faas,
            ch_new_annotation.gffs,
        )
    }

    if (params.mge_annotate) {
        MGE(
            ch_new_annotation.fastas
        )
    }

    workflow.onComplete {
        NextflowTool.summary(workflow, params, log)
    }
}
