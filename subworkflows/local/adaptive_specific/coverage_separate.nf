/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_SPLIT_BY_BED } from '../../../modules/local/samtools/main.nf'
include { CRAMINO_STATS         } from '../../../modules/local/cramino/main.nf'
include { MOSDEPTH_ADAPTIVE     } from '../../../modules/local/mosdepth/main.nf'
include { REMOVE_PADDING        } from '../../../modules/local/adaptive_specific/main.nf'
include { PIGZ_BED              } from '../../../modules/local/adaptive_specific/main.nf'
include { COVERAGE_PLOT         } from '../../../modules/local/adaptive_specific/main.nf'
include { modifyMetaId          } from '../utils_nfcore_oncoseq_pipeline'


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

    main:


    ch_versions = Channel.empty()

    // For now, we keep padding in bed file
    ch_bed = bed
        .map { meta,bedfile,_padding,_low_fidelity ->
            tuple(meta,bedfile) }

    ch_split_in = bam
        .join(ch_bed)
        

    SAMTOOLS_SPLIT_BY_BED(ch_split_in)

    ch_cramino_bg = SAMTOOLS_SPLIT_BY_BED.out.bg
        .map { meta, bamfile, bai ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_background')
            tuple(new_meta, bamfile, bai) }

    ch_cramino_panel = SAMTOOLS_SPLIT_BY_BED.out.panel
        .map { meta, bamfile, bai ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_panel')
            tuple(new_meta, bamfile, bai) }

    ch_cramino_in = ch_cramino_bg
        .mix(ch_cramino_panel)

    CRAMINO_STATS(ch_cramino_in)

    // Remove padding from bed file for further coverage computations
    ch_bed_pad = bed
        .map { meta,bedfile,padding,_low_fidelity ->
            tuple(meta,bedfile,padding) }
    REMOVE_PADDING(ch_bed_pad)

    // Create channels for running mosdepth with different filters for each sample by creating a new meta variable containing filter type,
        // joining the bed file without the padding and then adding the flag and MAPQ filters to the tuple
    // No filters on alignments:

    ch_bed_nopad_nofilt = REMOVE_PADDING.out.bed
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_nofilter')
            tuple(new_meta, bedfile) }

    ch_nofilt = SAMTOOLS_SPLIT_BY_BED.out.panel
        .map { meta, bamfile, bai ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_nofilter')
            tuple(new_meta, bamfile, bai) }
        .join(ch_bed_nopad_nofilt)
        .map { new_meta, bamfile, bai, bedfile ->
            tuple(new_meta, bamfile, bai, bedfile, 1540, 0) }

    // Primary alignments only:
    ch_bed_nopad_primary = REMOVE_PADDING.out.bed
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_primary')
            tuple(new_meta, bedfile) }

    ch_primary = SAMTOOLS_SPLIT_BY_BED.out.panel
        .map { meta, bamfile, bai ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_primary')
            tuple(new_meta, bamfile, bai) }
        .join(ch_bed_nopad_primary)
        .map { new_meta, bamfile, bai, bedfile ->
            tuple(new_meta, bamfile, bai, bedfile, 1796, 0) }

    // mapq60 alignments only:
    ch_bed_nopad_mapq60 = REMOVE_PADDING.out.bed
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_mapq60')
            tuple(new_meta, bedfile) }

    ch_mapq60 = SAMTOOLS_SPLIT_BY_BED.out.panel
        .map { meta, bamfile, bai ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_mapq60')
            tuple(new_meta, bamfile, bai) }
        .join(ch_bed_nopad_mapq60)
        .map { new_meta, bamfile, bai, bedfile ->
            tuple(new_meta, bamfile, bai, bedfile, 1796, 60) }

    mosdepth_in = ch_nofilt
        .mix(ch_primary,ch_mapq60)

    MOSDEPTH_ADAPTIVE(mosdepth_in)

    ch_coverage_bg = CRAMINO_STATS.out.stats
        .filter { meta, table -> meta.id.contains('background') }
        .map { meta, table ->
        // Read the file content as a list of lines
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_background')
            def lines = table.readLines()
            def coverage = lines[5].tokenize('\t')[1].toDouble()    // Last line and only take mean coverage column (4th)
            tuple(new_meta, coverage)
        }
    // collect each mosdepth adaptive output into it's own channel and join by orinal meta_id (sample_id) to produce plot

    PIGZ_BED(MOSDEPTH_ADAPTIVE.out.bed)

    ch_nofilt_bed_out = PIGZ_BED.out.bed
        .filter { meta, bedfile -> meta.id.contains('nofilter') }
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_nofilter')
            tuple(new_meta, bedfile) }

    ch_primary_bed_out = PIGZ_BED.out.bed
        .filter { meta, bedfile -> meta.id.contains('primary') }
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_primary')
            tuple(new_meta, bedfile) }

    ch_mapq60_bed_out = PIGZ_BED.out.bed
        .filter { meta, bedfile -> meta.id.contains('mapq60') }
        .map { meta, bedfile ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_mapq60')
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
    ch_versions = SAMTOOLS_SPLIT_BY_BED.out.versions
        .mix(CRAMINO_STATS.out.versions)
        .mix(MOSDEPTH_ADAPTIVE.out.versions)
        .mix(COVERAGE_PLOT.out.versions)



    emit:
    coverage_separated  = CRAMINO_STATS.out.stats              // TODO: QUARTO REPORT
    coverage_plot       = COVERAGE_PLOT.out.cov_plot_svg       // TODO: QUARTO REPORT
    versions            = ch_versions

}
