#!/usr/bin/env nextflow

include { CLUSTER_SEQS;
          ALIGNMENT_SCORE               } from '../modules/mmseqs/cluster/main'
include { KOFAMSCAN                     } from '../modules/kofamscan/annotate/main'
include { EGGNOGMAPPER                  } from '../modules/eggnogmapper/main'
include { DBCAN; RUNDBCAN_EASYSUBSTRATE } from '../modules/dbcan/main'
include { VFDB                          } from '../modules/vfdb/main'

workflow FUNC {
    take:
    annotation_faas
    annotation_gffs

    main:

    if ( params.cluster_proteome ) {
        /*
        collected_proteins = annotation_faas
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
        */
    } else {
        annotation_ready = annotation_faas
                        .map{ meta, faa -> faa}
    }

    if (!params.func_skip_eggnogmapper) {
        EGGNOGMAPPER(annotation_faas)
    }

    if (!params.func_skip_kofamscan) {
        KOFAMSCAN(annotation_faas)
    }

    if (!params.func_skip_vfdb) {
        VFDB(annotation_faas)
    }

    if ( params.mode == 'single' || params.mode == 'pyrodigal' || params.mode == 'bakta' || params.mode == 'prokka' && !params.func_skip_dbcan ) {
        linked_annotations = annotation_faas.join(annotation_gffs)

        RUNDBCAN_EASYSUBSTRATE(linked_annotations)
    } else {
        DBCAN(annotation_faas)
    }
    }
