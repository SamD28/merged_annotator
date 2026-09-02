#!/usr/bin/env nextflow

include { GENOMAD } from '../modules/genomad/main'

workflow MGE {
    take:
    annotation_fastas

    main:
    GENOMAD(annotation_fastas)
}
