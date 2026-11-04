/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BEDTOOLS_SUBTRACT               } from '../../../modules/nf-core/bedtools/subtract/main'
include { BEDTOOLS_MAKEWINDOWS            } from '../../../modules/nf-core/bedtools/makewindows/main'
include { MOSDEPTH as MOSDEPTH_BACKGROUND } from '../../../modules/nf-core/mosdepth'
include { MOSDEPTH as MOSDEPTH_NOFILT     } from '../../../modules/nf-core/mosdepth'
include { MOSDEPTH as MOSDEPTH_PRIMARY    } from '../../../modules/nf-core/mosdepth'
include { MOSDEPTH as MOSDEPTH_UNIQUE     } from '../../../modules/nf-core/mosdepth'
include { PIGZ_UNCOMPRESS                 } from '../../../modules/nf-core/pigz/uncompress/main'
include { COVERAGE_PLOT                   } from '../../../modules/local/adaptive_specific/main.nf'
// functions
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { modifyMetaId                 } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COVERAGE_SEPARATE {

    //TODO Add reports for coverage stats figure ?

    take:
    bam                     // channel: from mapping workflow, tuple with bam and bai
    bed                     // channel: from path read from params.bed, bed file used for adaptive sampling
    bed_nopad
    ref_fai

    main:


    ch_versions = channel.empty()

    // Create a bed file with window as big as chromosome size to substract ROI from:

    BEDTOOLS_MAKEWINDOWS(
        ref_fai
    )

    // For now, we keep padding in bed file, so that background coverage does not include padding
    ch_bed = bed
        .map { meta,bedfile,_padding,_low_fidelity ->
            tuple(meta,bedfile) }

    ch_in_bedtools = ch_bed
        .combine(BEDTOOLS_MAKEWINDOWS.out.bed.map { _meta, bed_bg -> bed_bg})
        .groupTuple(by:[1])

    BEDTOOLS_SUBTRACT(ch_in_bedtools.map { meta, bedfile, chrom_size -> tuple(meta, bedfile, chrom_size[0])})

    ch_bed_subtracted = BEDTOOLS_SUBTRACT.out.bed
        .transpose()

    ch_bam_bg = bam
        .join(ch_bed_subtracted)
        .map { meta, bamfile, bai, bed_bg ->
        def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_background')
        tuple(new_meta, bamfile, bai, bed_bg) }

    MOSDEPTH_BACKGROUND(
        ch_bam_bg,
        [[],[]],
        []
    )

    ch_in_mosdepth = bam
        .join(bed_nopad)

    // No filters on alignments:

    MOSDEPTH_NOFILT(
        ch_in_mosdepth.map { meta, bamfile, bai, bedfile ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_panel_nofilter')
            tuple(new_meta, bamfile, bai, bedfile)},
        [[],[]],
        []
    )

    // Primary alignments only:

    MOSDEPTH_PRIMARY(
        ch_in_mosdepth.map { meta, bamfile, bai, bedfile ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_panel_primary')
            tuple(new_meta, bamfile, bai, bedfile)},
        [[],[]],
        []
    )

    // mapq60 alignments only:

    MOSDEPTH_UNIQUE(
        ch_in_mosdepth.map { meta, bamfile, bai, bedfile ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_panel_mapq60')
            tuple(new_meta, bamfile, bai, bedfile)},
        [[],[]],
        []
    )

    ch_bed_cov = MOSDEPTH_NOFILT.out.regions_bed
        .mix(MOSDEPTH_PRIMARY.out.regions_bed)
        .mix(MOSDEPTH_UNIQUE.out.regions_bed)

    ch_coverage_bg = MOSDEPTH_BACKGROUND.out.summary_txt
        .map { meta, table ->
        // Read the file content as a list of lines
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_background')
            def lines = table.readLines()
            def coverage = lines.last().tokenize('\t')[3].toFloat()    // Last line and only take mean coverage column (4th)
            tuple(new_meta, coverage)
        }

    // collect each mosdepth adaptive output into it's own channel and join by orinal meta_id (sample_id) to produce plot

    PIGZ_UNCOMPRESS(ch_bed_cov)

    ch_nofilt_bed_out = PIGZ_UNCOMPRESS.out.file
        .filter { meta, bedfile -> meta.id.contains('nofilter') }
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_panel_nofilter')
            tuple(new_meta, bedfile) }

    ch_primary_bed_out = PIGZ_UNCOMPRESS.out.file
        .filter { meta, bedfile -> meta.id.contains('primary') }
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_panel_primary')
            tuple(new_meta, bedfile) }

    ch_mapq60_bed_out = PIGZ_UNCOMPRESS.out.file
        .filter { meta, bedfile -> meta.id.contains('mapq60') }
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_panel_mapq60')
            tuple(new_meta, bedfile) }

    // Now we join all the variables by sample_id (meta)

    ch_bed_joined = ch_nofilt_bed_out
        .join(ch_primary_bed_out)
        .join(ch_mapq60_bed_out)
        .join(ch_coverage_bg)
        .join(bed.map{ meta,_bedfile,_padding,low_fidelity -> tuple(meta,low_fidelity) })

    COVERAGE_PLOT(ch_bed_joined)

    //
    // Collate and save software versions
    //
    ch_versions = COVERAGE_PLOT.out.versions

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



    emit:
    coverage_plot       = COVERAGE_PLOT.out.cov_plot_svg
    coverage_tbl        = COVERAGE_PLOT.out.cov_df
    versions            =  ch_collated_versions

}
