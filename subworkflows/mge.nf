#!/usr/bin/env nextflow

include { GENOMAD } from '../modules/genomad/main'

workflow MGE {
    take:
    annotation_faas

    main:
    GENOMAD(annotation_faas)
}
