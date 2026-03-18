#!/usr/bin/env nextflow

include { CLUSTER_SEQS;
          ALIGNMENT_SCORE  } from '../modules/mmseqs/cluster/main'
include { KOFAMSCAN        } from '../modules/kofamscan/annotate/main'
include { EGGNOGMAPPER     } from '../modules/eggnogmapper/main'
include { DBCAN            } from '../modules/dbcan/main'
include { VFDB             } from '../modules/vfdb/main'

workflow FUNC {
    take:
    annotation_results

    main:

    if ( params.cluster_proteome ) {
        collected_proteins = annotation_results
            .map{ _meta, faa -> faa}
            .collectFile( name: 'all_proteins.faa', newLine: true )
            
        (intermediate_cluster_ch, db_path) = CLUSTER_SEQS(collected_proteins)

        cluster_ch = intermediate_cluster_ch
                    .multiMap { all_seqs, clustering_tsv, rep_seq ->
                        all_seqs: all_seqs
                        clustering_tsv: clustering_tsv
                        rep_seq: rep_seq
                    }

        ALIGNMENT_SCORE(db_path)

        annotation_ready = cluster_ch.rep_seq
    } else {
        annotation_ready = annotation_results
                        .map{ _meta, faa -> faa}
    }

    kofamscan_result  = KOFAMSCAN(annotation_ready)

    eggnog_result     = EGGNOGMAPPER(annotation_ready)

    dbcan_result      = DBCAN(annotation_ready)

    vfdb_result       = VFDB(annotation_ready)
    }
