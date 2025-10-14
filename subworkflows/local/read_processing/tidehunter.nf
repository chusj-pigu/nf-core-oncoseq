/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { TIDEHUNTER_FASTQ } from '../../../modules/local/tidehunter/main.nf'
include { TIDEHUNTER_TAB   } from '../../../modules/local/tidehunter/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TIDEHUNTER_CONCENSUS {

    //TODO Add reports for coverage stats figure ?

    take:
    reads

    main:

    TIDEHUNTER_FASTQ(reads)
    TIDEHUNTER_TAB(reads)

    ch_versions = TIDEHUNTER_FASTQ.out.versions
        .mix(TIDEHUNTER_TAB.out.versions)

    emit:
    fastq               = TIDEHUNTER_FASTQ.out.fq
    summary             = TIDEHUNTER_TAB.out.tsv     // TODO: Quarto report
    versions            = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
