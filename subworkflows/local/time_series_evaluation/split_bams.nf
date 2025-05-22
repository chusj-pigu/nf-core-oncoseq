/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import the ONTIME_RANGE_FILTER module, which is responsible for filtering or
    splitting BAM files based on specific criteria.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ONTIME_RANGE_FILTER } from '../../../modules/local/ontime/main.nf'
include { SAMTOOLS_SORT_INDEX } from '../../../modules/local/samtools/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPLIT_BAMS_TIME SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This subworkflow splits BAM files using the ONTIME_RANGE_FILTER module.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPLIT_BAMS_TIME {

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ch_bam: Input channel containing BAM files to be split.
    */
    take:
    ch_bam
    ref

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PROCESS BAM FILES USING ONTIME_RANGE_FILTER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        The ONTIME_RANGE_FILTER module processes the input BAM files and outputs
        the filtered or split BAM files.
    */

    // TODO: move this to an input channel
    timeseries = Channel.from(3, 6, 12, 24, 48)
        .map { x -> tuple(0, "${x}") }

    ch_bam_ts = ch_bam
        .combine(timeseries)
        .map { meta, bam, index, from, to ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_${from}h_${to}h"
            new_meta.ts = "${from}h_${to}h"
            tuple(new_meta, bam, index, from, to)
        }
    // Add time series to the ref data
    ref = ref
        .combine(timeseries)
        .map { meta, refid, ref_fasta, ref_fai, from, to  ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_${from}h_${to}h"
            new_meta.ts = "${from}h_${to}h"
            tuple(new_meta, refid, ref_fasta, ref_fai)
        }

    // Keep only the meta from the ref data
    ref = ref
        .map { meta, refid, ref_fasta, ref_fai ->
            tuple(meta, refid, ref_fasta, ref_fai)
        }

    ONTIME_RANGE_FILTER(
        ch_bam_ts
    )

    SAMTOOLS_SORT_INDEX(
        ONTIME_RANGE_FILTER.out.bam
    )

    // Add FULL to the original bam
    ch_bam = ch_bam
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_FULL"
            new_meta.ts = "FULL"
            tuple(new_meta, bam, bai)
        }

    ch_bam_out = SAMTOOLS_SORT_INDEX.out.sortedbamidx
        .map { meta, bam, bai ->
            tuple(meta, bam, bai)
        }.mix(ch_bam)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSION INFORMATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Define an empty channel for collecting version information. This can be
        populated in the future to track tool versions for reproducibility.
    */
    ch_versions = SAMTOOLS_SORT_INDEX.out.versions
    ch_versions.mix(ONTIME_RANGE_FILTER.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ch_bam_out: Channel containing the processed BAM files.
        ch_versions: Channel for version information (currently empty).
    */
    emit:
    bam              = ch_bam_out
    versions         = ch_versions
    ref
}
