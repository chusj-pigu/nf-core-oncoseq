/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPLIT_BAMS_TIME SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This subworkflow is responsible for splitting BAM files into time series and full
    intervals, and propagating associated reference, chromosome list, model, and clinical
    database information for each interval. It also collects version information for
    reproducibility. All original inputs are printed for debugging.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ONTIME_RANGE_FILTER } from '../../../modules/local/ontime/main.nf' // Module for filtering/splitting BAMs by time
include { SAMTOOLS_SORT_INDEX } from '../../../modules/local/samtools/main.nf' // Module for sorting and indexing BAMs

workflow SPLIT_BAMS_TIME {
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ch_bam:         Channel of BAM files to be split (tuple: meta, bam, bai)
        ref:            Channel of reference data (tuple: meta, refid, ref_fasta, ref_fai)
        clin_database:  Channel of clinical database (tuple: meta, clin_database_val)
    */
    take:
    ch_bam
    ref
    clin_database

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DEFINE TIME SERIES INTERVALS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        The time intervals (in hours) for splitting BAMs and propagating meta info.
        This is currently hardcoded but could be parameterized.
    */
    timeseries = Channel.from(3, 6, 12, 24, 48)
        .map { x -> tuple(0, "${x}") } // Each tuple: (from, to)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPLIT BAM FILES INTO TIME SERIES INTERVALS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        For each BAM, create a new meta for each time interval, updating the id and ts.
        Output: ch_bam_ts (tuple: new_meta, bam, index, from, to)
    */
    ch_bam_ts = ch_bam
        .combine(timeseries)
        .map { meta, bam, index, from, to ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_${from}h_${to}h"
            new_meta.ts = "${from}h_${to}h"
            tuple(new_meta, bam, index, from, to)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PROPAGATE TIME SERIES TO REFERENCE DATA
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        For each ref entry, create a new meta for each time interval, updating id/ts.
        Output: ref (tuple: new_meta, refid, ref_fasta, ref_fai)
    */
    ref = ref
        .combine(timeseries)
        .map { meta, refid, ref_fasta, ref_fai, from, to  ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_${from}h_${to}h"
            new_meta.ts = "${from}h_${to}h"
            tuple(new_meta, refid, ref_fasta, ref_fai)
        }
    // Keep only the meta and reference fields
    ref = ref
        .map { meta, refid, ref_fasta, ref_fai ->
            tuple(meta, refid, ref_fasta, ref_fai)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FILTER/SPLIT BAM FILES BY TIME INTERVALS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Use the ONTIME_RANGE_FILTER module to split/filter BAMs by time intervals.
        Output: ONTIME_RANGE_FILTER.out.bam (filtered BAMs)
    */
    ONTIME_RANGE_FILTER(
        ch_bam_ts
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SORT AND INDEX FILTERED BAM FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Use the SAMTOOLS_SORT_INDEX module to sort and index the filtered BAMs.
        Output: SAMTOOLS_SORT_INDEX.out.sortedbamidx (sorted BAMs with index)
    */
    SAMTOOLS_SORT_INDEX(
        ONTIME_RANGE_FILTER.out.bam
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ADD 'FULL' META ENTRIES FOR ORIGINAL INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        For each input (bam, ref, clin_database), add a 'FULL' meta
        entry to represent the unfiltered/original data. These are mixed into the
        respective channels so downstream steps receive both time series and full data.
    */
    // Add FULL to the original bam
    ch_bam = ch_bam
        .map { meta, bam, bai ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_FULL"
            new_meta.ts = "FULL"
            tuple(new_meta, bam, bai)
        }
    // Add FULL to the original ref
    ref_full = ref
        .map { meta, refid, ref_fasta, ref_fai ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_FULL"
            new_meta.ts = "FULL"
            tuple(new_meta, refid, ref_fasta, ref_fai)
        }
    ref = ref.mix(ref_full)
    // Add FULL to the original clin_database
    clin_database.view()
    clin_database_full = clin_database
        .map { meta, clin_database_val ->
            def new_meta = meta.clone()
            new_meta.id = "${meta.id}_FULL"
            new_meta.ts = "FULL"
            tuple(new_meta, clin_database_val)
        }
    clin_database = clin_database.mix(clin_database_full)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COMBINE PROCESSED AND FULL BAM OUTPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Combine sorted/indexed BAMs from time series with the original FULL BAMs.
        Output: ch_bam_out (tuple: meta, bam, bai)
    */
    ch_bam_out = SAMTOOLS_SORT_INDEX.out.sortedbamidx
        .map { meta, bam, bai ->
            tuple(meta, bam, bai)
        }.mix(ch_bam)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSION INFORMATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collects version information from both SAMTOOLS_SORT_INDEX and
        ONTIME_RANGE_FILTER modules for reproducibility.
    */
    ch_versions = SAMTOOLS_SORT_INDEX.out.versions
    ch_versions.mix(ONTIME_RANGE_FILTER.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        bam:           Channel containing all processed BAM files (time series + full)
        versions:      Channel for version information (currently empty)
        ref:           Channel containing all reference data (time series + full)
        clin_database: Channel containing all clinical database data (time series + full)
    */
    emit:
    bam              = ch_bam_out
    versions         = ch_versions
    ref              = ref
    clin_database    = clin_database
}
