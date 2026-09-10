nextflow.enable.dsl = 2

include { ENSEMBLVEP_VEP } from './modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_DOWNLOAD } from './modules/nf-core/ensemblvep/download/main.nf'

workflow {}
