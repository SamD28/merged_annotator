include { MANIFEST_PARSE   } from './subworkflows/manifest_parse.nf'
include { ANNOTATION       } from './subworkflows/annotation.nf'
include { ARG              } from './subworkflows/arg.nf'
include { FUNC             } from './subworkflows/func.nf'

workflow {
    main:
    ch_input_for_annotation =  MANIFEST_PARSE(params.manifest)
    
    ANNOTATION(ch_input_for_annotation)

    ch_new_annotation = ch_input_for_annotation
            .join(ANNOTATION.out.faa)
            .join(ANNOTATION.out.gbk)
            .multiMap { meta, fasta, faa, gbk ->
                fastas: [meta, fasta]
                faas: [meta, faa]
                gbks: [meta, gbk]
            }

    if (params.arg_annotate) {
        ARG(
            ch_new_annotation.fastas,
            ch_new_annotation.faas
        )
    }

    if (params.func_annotate) {
        FUNC(ch_new_annotation.faas)
    }
}