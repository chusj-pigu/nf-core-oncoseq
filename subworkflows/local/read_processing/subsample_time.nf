/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ONTIME_RANGE_FILTER       } from '../../../modules/local/ontime/main.nf' // Module for filtering/splitting BAMs by time
include { ONTIME_RANGE_FILTER_FASTQ } from '../../../modules/local/ontime/main.nf' // Module for filtering/splitting BAMs by time
include { SAMTOOLS_SORT             } from '../../../modules/local/samtools/main.nf' // Module for sorting BAMs
include { SAMTOOLS_INDEX            } from '../../../modules/local/samtools/main.nf' // Module for indexing BAMs

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SUBSAMPLE_TIME {

    //TODO Add reports for coverage stats figure ?

    take:
    input_time            // Channel with bam, bai or fastq, target sequencing start time and end time

    main:

    // When using Marlin or remap to t2t

    ONTIME_RANGE_FILTER(input_time)

    SAMTOOLS_SORT(ONTIME_RANGE_FILTER.out.bam)

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sortedbam)

    ch_versions = ONTIME_RANGE_FILTER.out.versions
        .mix(SAMTOOLS_SORT.out.versions)
        .mix(SAMTOOLS_INDEX.out.versions)

    ch_out_bam = SAMTOOLS_INDEX.out.bamfile_index

    emit:
    bam                 = ch_out_bam
    versions            = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
