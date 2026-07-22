/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BEDTOOLS_SUBTRACT                             } from '../../../modules/local/bedtools/main.nf'
include { MOSDEPTH_ADAPTIVE as MOSDEPTH_BACKGROUND      } from '../../../modules/local/mosdepth/main.nf'
include { MOSDEPTH_ADAPTIVE as MOSDEPTH_PANEL           } from '../../../modules/local/mosdepth/main.nf'
include { REMOVE_PADDING                                } from '../../../modules/local/adaptive_specific/main.nf'
include { PIGZ_BED                                      } from '../../../modules/local/adaptive_specific/main.nf'
include { COVERAGE_PLOT                                 } from '../../../modules/local/adaptive_specific/main.nf'
include { modifyMetaId                                  } from '../utils_nfcore_oncoseq_pipeline'


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
    ref

    main:


    ch_versions = channel.empty()

    // For now, we keep padding in bed file
    ch_bed = bed
        .map { meta,bedfile,_padding,_low_fidelity ->
            tuple(meta,bedfile) }

    ch_ref_id = ref
        .map { meta, ref_id, _ref_fa, _ref_index ->
            tuple(meta,ref_id)}

    ch_in_bedtools = ch_bed
        .join(ch_ref_id)

    BEDTOOLS_SUBTRACT(ch_in_bedtools)

    ch_bam_panel = bam
        .map { meta, bamfile, bai ->
        def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_panel')
        tuple(new_meta, bamfile, bai) }

    ch_bam_bg = bam
        .join(BEDTOOLS_SUBTRACT.out.bed)
        .map { meta, bamfile, bai, bed_bg ->
        def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_background')
        tuple(new_meta, bamfile, bai, bed_bg, 1796, 0) }

    MOSDEPTH_BACKGROUND(ch_bam_bg)

    // Create channels for running mosdepth with different filters for each sample by creating a new meta variable containing filter type,
        // joining the bed file without the padding and then adding the flag and MAPQ filters to the tuple
    // No filters on alignments:

    ch_nofilt = bam
        .join(bed_nopad)
        .map { meta, bamfile, bai, bedfile ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_panel_nofilter')
            tuple(new_meta, bamfile, bai, bedfile) }
        .map { new_meta, bamfile, bai, bedfile ->
            tuple(new_meta, bamfile, bai, bedfile, 1540, 0) }

    // Primary alignments only:

    ch_primary = bam
        .join(bed_nopad)
        .map { meta, bamfile, bai, bedfile ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_panel_primary')
            tuple(new_meta, bamfile, bai, bedfile) }
        .map { new_meta, bamfile, bai, bedfile ->
            tuple(new_meta, bamfile, bai, bedfile, 1796, 0) }

    // mapq60 alignments only:

    ch_mapq60 = bam
        .join(bed_nopad)
        .map { meta, bamfile, bai, bedfile ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_panel_mapq60')
            tuple(new_meta, bamfile, bai, bedfile) }
        .map { new_meta, bamfile, bai, bedfile ->
            tuple(new_meta, bamfile, bai, bedfile, 1796, 60) }

    mosdepth_in = ch_nofilt
        .mix(ch_primary,ch_mapq60)

    MOSDEPTH_PANEL(mosdepth_in)

    ch_coverage_bg = MOSDEPTH_BACKGROUND.out.summary
        .map { meta, table ->
        // Read the file content as a list of lines
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_background')
            def lines = table.readLines()
            def coverage = lines.last().tokenize('\t')[3].toFloat()    // Last line and only take mean coverage column (4th)
            tuple(new_meta, coverage)
        }
        
    // collect each mosdepth adaptive output into it's own channel and join by orinal meta_id (sample_id) to produce plot

    PIGZ_BED(MOSDEPTH_PANEL.out.bed)

    ch_nofilt_bed_out = PIGZ_BED.out.bed
        .filter { meta, bedfile -> meta.id.contains('nofilter') }
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_panel_nofilter')
            tuple(new_meta, bedfile) }

    ch_primary_bed_out = PIGZ_BED.out.bed
        .filter { meta, bedfile -> meta.id.contains('primary') }
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_panel_primary')
            tuple(new_meta, bedfile) }

    ch_mapq60_bed_out = PIGZ_BED.out.bed
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
    ch_versions = BEDTOOLS_SUBTRACT.out.versions
        .mix(MOSDEPTH_PANEL.out.versions)
        .mix(COVERAGE_PLOT.out.versions)



    emit:
    coverage_plot       = COVERAGE_PLOT.out.cov_plot_svg
    coverage_tbl        = COVERAGE_PLOT.out.cov_df
    split_bed           = bed_nopad
    versions            = ch_versions

}
