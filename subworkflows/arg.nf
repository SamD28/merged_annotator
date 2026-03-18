/*
    Run ARG screening tools
*/

include { ABRICATE_RUN                     } from '../modules/abricate/run/main'
include { AMRFINDERPLUS_UPDATE             } from '../modules/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN                } from '../modules/amrfinderplus/run/main'
include { DEEPARG_DOWNLOADDATA             } from '../modules/deeparg/downloaddata/main'
include { DEEPARG_PREDICT                  } from '../modules/deeparg/predict/main'
include { FARGENE                          } from '../modules/fargene/main'
include { RGI_CARDANNOTATION               } from '../modules/rgi/cardannotation/main'
include { RGI_MAIN                         } from '../modules/rgi/main/main'
include { UNTAR as UNTAR_CARD              } from '../modules/untar/main'
include { TABIX_BGZIP as ARG_TABIX_BGZIP   } from '../modules/tabix/bgzip/main'
include { HAMRONIZATION_RGI                } from '../modules/hamronization/rgi/main'
include { HAMRONIZATION_FARGENE            } from '../modules/hamronization/fargene/main'
include { HAMRONIZATION_SUMMARIZE          } from '../modules/hamronization/summarize/main'
include { HAMRONIZATION_ABRICATE           } from '../modules/hamronization/abricate/main'
include { HAMRONIZATION_DEEPARG            } from '../modules/hamronization/deeparg/main'
include { HAMRONIZATION_AMRFINDERPLUS      } from '../modules/hamronization/amrfinderplus/main'
include { ARGNORM as ARGNORM_DEEPARG       } from '../modules/argnorm/main'
include { ARGNORM as ARGNORM_ABRICATE      } from '../modules/argnorm/main'
include { ARGNORM as ARGNORM_AMRFINDERPLUS } from '../modules/argnorm/main'

workflow ARG {
    take:
    fastas      // tuple val(meta), path(contigs)
    annotations

    main:
    ch_versions = Channel.empty()

    // Prepare HAMRONIZATION reporting channel
    ch_input_to_hamronization_summarize = Channel.empty()

    // AMRfinderplus run
    // Prepare channel for database
    if (!params.arg_skip_amrfinderplus && params.arg_amrfinderplus_db) {
        ch_amrfinderplus_db = Channel
            .fromPath(params.arg_amrfinderplus_db, checkIfExists: true)
            .first()
    }
    else if (!params.arg_skip_amrfinderplus && !params.arg_amrfinderplus_db) {
        AMRFINDERPLUS_UPDATE()
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)
        ch_amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
    }

    if (!params.arg_skip_amrfinderplus) {
        AMRFINDERPLUS_RUN(fastas, ch_amrfinderplus_db)
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)

        // Reporting
        HAMRONIZATION_AMRFINDERPLUS(AMRFINDERPLUS_RUN.out.report, 'tsv', AMRFINDERPLUS_RUN.out.tool_version, AMRFINDERPLUS_RUN.out.db_version)
        ch_versions = ch_versions.mix(HAMRONIZATION_AMRFINDERPLUS.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_AMRFINDERPLUS.out.tsv)

        if (!params.arg_skip_argnorm) {
            ch_input_to_argnorm_amrfinderplus = HAMRONIZATION_AMRFINDERPLUS.out.tsv.filter { meta, file -> !file.isEmpty() }
            ARGNORM_AMRFINDERPLUS(ch_input_to_argnorm_amrfinderplus, 'amrfinderplus', 'ncbi')
            ch_versions = ch_versions.mix(ARGNORM_AMRFINDERPLUS.out.versions)
        }
    }

    // fARGene run
    if (!params.arg_skip_fargene) {
        ch_fargene_classes = Channel.fromList(params.arg_fargene_hmmmodel.tokenize(','))

        ch_fargene_input = fastas
            .combine(ch_fargene_classes)
            .map { meta, fastas, hmm_class ->
                def meta_new = meta.clone()
                meta_new['hmm_class'] = hmm_class
                [meta_new, fastas, hmm_class]
            }
            .multiMap {
                fastas: [it[0], it[1]]
                hmmclass: it[2]
            }

        FARGENE(ch_fargene_input.fastas, ch_fargene_input.hmmclass)
        ch_versions = ch_versions.mix(FARGENE.out.versions)

        // Reporting
        // Note: currently hardcoding versions, has to be updated with every fARGene-update
        HAMRONIZATION_FARGENE(FARGENE.out.hmm_genes.transpose(), 'tsv', '0.1', '0.1')
        ch_versions = ch_versions.mix(HAMRONIZATION_FARGENE.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_FARGENE.out.tsv)
    }

    // RGI run
    if (!params.arg_skip_rgi) {

        if (!params.arg_rgi_db) {

            // Download and untar CARD
            UNTAR_CARD(file('https://card.mcmaster.ca/latest/data', checkIfExists: true))
            ch_versions = ch_versions.mix(UNTAR_CARD.out.versions)
            rgi_db = UNTAR_CARD.out.untar
            RGI_CARDANNOTATION(rgi_db)
            card = RGI_CARDANNOTATION.out.db
            ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions)
        }
        else {

            // Use user-supplied database
            rgi_db = file(params.arg_rgi_db, checkIfExists: true)
            if (!rgi_db.contains("card_database_processed")) {
                RGI_CARDANNOTATION(rgi_db)
                card = RGI_CARDANNOTATION.out.db
                ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions)
            }
            else {
                card = rgi_db
            }
        }

        RGI_MAIN(fastas, card, [])
        ch_versions = ch_versions.mix(RGI_MAIN.out.versions)

        // Reporting
        HAMRONIZATION_RGI(RGI_MAIN.out.tsv, 'tsv', RGI_MAIN.out.tool_version, RGI_MAIN.out.db_version)
        ch_versions = ch_versions.mix(HAMRONIZATION_RGI.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_RGI.out.tsv)
    }

    // DeepARG prepare download
    if (!params.arg_skip_deeparg && params.arg_deeparg_db) {
        ch_deeparg_db = Channel
            .fromPath(params.arg_deeparg_db, checkIfExists: true)
            .first()
    }
    else if (!params.arg_skip_deeparg && !params.arg_deeparg_db) {
        DEEPARG_DOWNLOADDATA()
        ch_versions = ch_versions.mix(DEEPARG_DOWNLOADDATA.out.versions)
        ch_deeparg_db = DEEPARG_DOWNLOADDATA.out.db
    }

    // DeepARG run
    if (!params.arg_skip_deeparg) {

        annotations
            .map { it ->
                def meta = it[0]
                def anno = it[1]
                def model = params.arg_deeparg_model

                [meta, anno, model]
            }
            .set { ch_input_for_deeparg }

        DEEPARG_PREDICT(ch_input_for_deeparg, ch_deeparg_db)
        ch_versions = ch_versions.mix(DEEPARG_PREDICT.out.versions)

        // Reporting
        // Note: currently hardcoding versions as unreported by DeepARG
        // Make sure to update on version bump.
        ch_input_to_hamronization_deeparg = DEEPARG_PREDICT.out.arg.mix(DEEPARG_PREDICT.out.potential_arg)
        HAMRONIZATION_DEEPARG(ch_input_to_hamronization_deeparg, 'tsv', '1.0.4', params.arg_deeparg_db_version)
        ch_versions = ch_versions.mix(HAMRONIZATION_DEEPARG.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_DEEPARG.out.tsv)

        if (!params.arg_skip_argnorm) {
            ch_input_to_argnorm_deeparg = HAMRONIZATION_DEEPARG.out.tsv.filter { meta, file -> !file.isEmpty() }
            ARGNORM_DEEPARG(ch_input_to_argnorm_deeparg, 'deeparg', 'deeparg')
            ch_versions = ch_versions.mix(ARGNORM_DEEPARG.out.versions)
        }
    }

    // ABRicate run
    if (!params.arg_skip_abricate) {
        abricate_dbdir = params.arg_abricate_db ? file(params.arg_abricate_db, checkIfExists: true) : []
        ABRICATE_RUN(fastas, abricate_dbdir)
        ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions)

        HAMRONIZATION_ABRICATE(ABRICATE_RUN.out.report, 'tsv', '1.0.1', '2021-Mar-27')
        ch_versions = ch_versions.mix(HAMRONIZATION_ABRICATE.out.versions)
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_ABRICATE.out.tsv)

        if ((params.arg_abricate_db_id == 'ncbi' || params.arg_abricate_db_id == 'resfinder' || params.arg_abricate_db_id == 'argannot' || params.arg_abricate_db_id == 'megares') && !params.arg_skip_argnorm) {
            ch_input_to_argnorm_abricate = HAMRONIZATION_ABRICATE.out.tsv.filter { meta, file -> !file.isEmpty() }
            ARGNORM_ABRICATE(ch_input_to_argnorm_abricate, 'abricate', params.arg_abricate_db_id)
            ch_versions = ch_versions.mix(ARGNORM_ABRICATE.out.versions)
        }
    }

    ch_input_to_hamronization_summarize
        .map {
            it[1]
        }
        .collect()
        .set { ch_input_for_hamronization_summarize }

    HAMRONIZATION_SUMMARIZE(ch_input_for_hamronization_summarize, params.arg_hamronization_summarizeformat)
    ch_versions = ch_versions.mix(HAMRONIZATION_SUMMARIZE.out.versions)

    emit:
    versions = ch_versions
}
