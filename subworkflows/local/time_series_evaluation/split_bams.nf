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
include { SAMTOOLS_SORT       } from '../../../modules/nf-core/samtools/sort/main' // Module for sorting BAMs
include { softwareVersionsToYAML                    } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { modifyMetaId                              } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

workflow SPLIT_BAMS_TIME {
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ch_bam:         Channel of BAM files to be split (tuple: meta, bam, bai)
        ref:            Channel of reference data (tuple: meta, refid, ref_fasta, ref_fai)
    */

    take:
    ch_bam
    bed

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DEFINE TIME SERIES INTERVALS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        The time intervals (in hours) for splitting BAMs and propagating meta info.
        This is currently hardcoded but could be parameterized.
    */
    def createTimeInterval = { x -> tuple(0, "${x}") }

    timeseries = channel
        .from(params.time_points.split(','))
        .map { it -> it.trim() } // optional: remove surrounding spaces
        .map(createTimeInterval) // Each tuple: (from, to)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPLIT BAM FILES INTO TIME SERIES INTERVALS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        For each BAM, create a new meta for each time interval, updating the id and ts.
        Output: ch_bam_ts (tuple: new_meta, bam, index, from, to)
    */
    def createBamTimeSeriesMeta = { meta, bam, index, from, to ->
        def to_stripped = to.replace('"', '')
        def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_${from}h_${to_stripped}h")
        new_meta.ts = "${from}h_${to_stripped}h".toString()
        tuple(new_meta, bam, index, from, to_stripped)
    }

    ch_bam_ts = ch_bam
        .combine(timeseries)
        .map(createBamTimeSeriesMeta)

        /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PROPAGATE TIME SERIES TO BED DATA
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        For each bed entry, create a new meta for each time interval, updating id/ts.
        Output: bed (tuple: new_meta, bedfile, padding, path)
    */

    def createBedTimeSeriesMeta = { meta, bedfile, padding, path, from, to ->
        def to_stripped = to.replace('"', '')
        def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_${from}h_${to_stripped}h")
        new_meta.ts = "${from}h_${to_stripped}h".toString()
        tuple(new_meta, bedfile, padding, path)
    }

    ts_bed = bed
        .combine(timeseries)
        .map(createBedTimeSeriesMeta)

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
    SAMTOOLS_SORT(
        ONTIME_RANGE_FILTER.out.bam,
        [[], [], []],
        'bai'
    )

    ch_bam_sorted = SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_SORT.out.index)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ADD 'FULL' META ENTRIES FOR ORIGINAL INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        For each input (bam, ref, clin_database), add a 'FULL' meta
        entry to represent the unfiltered/original data. These are mixed into the
        respective channels so downstream steps receive both time series and full data.
    */

    if (params.include_full ) {
        // Add FULL to the original bam
        def toFullBam = { meta, bam, bai ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_FULL')
            new_meta.ts = "FULL"
            tuple(new_meta, bam, bai)
        }

        ch_bam = ch_bam
            .map(toFullBam)

        // Add FULL to the original bed
        def toFullBed = { meta, bedfile, padding, path ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_FULL')
            new_meta.ts = "FULL"
            tuple(new_meta, bedfile, padding, path)
        }

        ch_bed = bed
            .map(toFullBed)

        ch_bed = ts_bed.mix(ch_bed)

    } else {
        // If not including full BAMs, just use the time series reference
        ch_bed = ts_bed
    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COMBINE PROCESSED AND FULL BAM OUTPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Combine sorted/indexed BAMs from time series with the original FULL BAMs.
        Output: ch_bam_out (tuple: meta, bam, bai)
    */
    def extractBamFields = { meta, bam, bai ->
        tuple(meta, bam, bai)
    }



    ch_bam_out = ch_bam_sorted
        .map(extractBamFields)
        .mix(ch_bam)


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSION INFORMATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collects version information from both SAMTOOLS_SORT_INDEX and
        ONTIME_RANGE_FILTER modules for reproducibility.
    */
    ch_versions = ONTIME_RANGE_FILTER.out.versions

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        bam:           Channel containing all processed BAM files (time series + full)
        versions:      Channel for version information (currently empty)
        ref:           Channel containing all reference data (time series + full)
    */

    emit:
    bam              = ch_bam_out
    versions         = ch_collated_versions
    bed              = ch_bed
}
