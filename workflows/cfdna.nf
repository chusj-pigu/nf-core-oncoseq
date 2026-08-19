include { BASECALL_MULTIPLEX     } from '../subworkflows/local/basecalling/basecall_multiplex'
include { BASECALL_SIMPLEX       } from '../subworkflows/local/basecalling/basecall_simplex'
include { READS_FILTER           } from '../subworkflows/local/read_processing/reads_filter.nf'
include { MAPPING                } from '../subworkflows/local/mapping/mapping'
include { CNV_CALLING            } from '../subworkflows/local/variant_calling/cnv_calling.nf'
include { SV_CALLING             } from  '../subworkflows/local/variant_calling/sv_calling.nf'
include { CLAIRS_TO_CALLING      } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'
include { CLAIR3_CALLING         } from '../subworkflows/local/variant_calling/clair3_calling.nf'
include { VARIANT_PROCESS        } from  '../subworkflows/local/variant_calling/variant_process.nf'
include { ICHORCNA_CALLING       } from '../subworkflows/local/variant_calling/ichor_calling.nf'
include { CLASSY                 } from '../subworkflows/local/methylation_analysis/classy.nf'
include { SUBSAMPLE_TIME as SUBSAMPLE_TIME_BAM  } from '../subworkflows/local/read_processing/subsample_time.nf'
include { SAMTOOLS_COUNT_READS                  } from '../modules/local/samtools/main.nf'
include { COVERAGE_SEPARATE                     } from '../subworkflows/local/adaptive_specific/coverage_separate'
include { SPLIT_BAMS_TIME as SPLIT_BAMS_CLASSIFY } from '../subworkflows/local/time_series_evaluation/split_bams.nf'
include { SPLIT_BAMS_TIME as SPLIT_BAMS_OTHER    } from '../subworkflows/local/time_series_evaluation/split_bams.nf'
include { CRAMINO_STATS          } from '../modules/local/cramino/main.nf'
include { REMOVE_PADDING                                } from '../modules/local/adaptive_specific/main.nf'

// Reporting
include { MIDNIGHT_REPORT      } from '../subworkflows/local/report/final_report.nf'
include { CFDNA_REPORT         } from '../subworkflows/local/report/cfdna.nf'
include { CLASSIFIER_REPORT    } from '../subworkflows/local/report/methylation.nf'
include { FIGENO_REPORT        } from '../subworkflows/local/report/variants.nf'
include { ADAPTIVE_REPORT      } from '../subworkflows/local/report/adaptive.nf'
include { ONTIME_TIME_RANGE    } from '../modules/local/ontime/main.nf'

// Useful functions to handle time parsing
def parseToInstant(str) {
    try {
        return java.time.Instant.parse(str)
    } catch (e) {}
    try {
        return java.time.OffsetDateTime.parse(str).toInstant()
    } catch (e) {}
    try {
        def fmt = java.time.format.DateTimeFormatter
            .ofPattern("yyyy-MM-dd'T'HH:mm:ss.SSSz")
            .withZone(java.time.ZoneId.of("UTC"))
        return java.time.ZonedDateTime.parse(str, fmt).toInstant()
    } catch (e) {}
    throw new RuntimeException("Could not parse timestamp: ${str}")
}

workflow CFDNA {

    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    tumor_type  // channel: tumor type read in from samplesheet
    cfdna_samplesheet   // channel : from demux or samplesheeet
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    max_len
    minqs
    ichor_bin
    mapq_wig
    basecall_model          // channel: model for basecalling
    clairs_model
    bed                     // channel: bed file used for adaptive sampling regions
    targets
    vep_cache

    main:

    //
    // WORKFLOW: Run pipeline
    //

    ch_versions = channel.empty()
    ch_sections = channel.empty()

    // Process bed

    ch_bed_pad = bed
        .map { meta,bedfile,padding,_low_fidelity ->
            tuple(meta,bedfile,padding) }
        .groupTuple(by:[1,2])

    REMOVE_PADDING(ch_bed_pad)

    ch_bed_nopad = REMOVE_PADDING.out.bed
        .transpose()

    if (params.demux) {

        BASECALL_MULTIPLEX (
            samplesheet,
            demux
        )

        ch_versions = BASECALL_MULTIPLEX.out.versions

        READS_FILTER (
            cfdna_samplesheet,
            BASECALL_MULTIPLEX.out.fastq,
            max_len,
            minqs
        )

        MAPPING(
            READS_FILTER.out.reads,
            ref
        )

        if (params.m_bases) {
            ch_to_classify = MAPPING.out.bam
                .map { meta, _bam, _bai -> meta }
            ch_no_methylation_tags = channel.empty()
        } else {
            ch_no_methylation_tags = MAPPING.out.bam
                .map { meta, _bam, _bai -> meta }
            ch_to_classify = channel.empty()
        }

    } else if (params.skip_basecalling || params.skip_mapping) {

        READS_FILTER (
            cfdna_samplesheet,
            samplesheet,
            max_len,
            minqs
        )

        ch_fastq_processed = READS_FILTER.out.reads

        MAPPING(
            ch_fastq_processed,
            ref
        )

        SAMTOOLS_COUNT_READS(MAPPING.out.bam.map { meta, bam, _bai -> tuple(meta, bam)})

        ch_bam_methylation_counts = SAMTOOLS_COUNT_READS.out.txt
            .map { meta, txt ->
                def count = txt.text.trim().toInteger()
                tuple(meta, count)
            }
            .branch { meta, count ->
                pos:  count > 0
                    return meta
                none: true
                    return meta
            }

        ch_bam_methylation_counts.none
            .subscribe { meta ->
                log.warn "Tumor classification will be skipped for ${meta.id} -- no methylation tags found in bam"
            }
        ch_to_classify = ch_bam_methylation_counts.pos
        ch_no_methylation_tags = ch_bam_methylation_counts.none

    } else {

         BASECALL_SIMPLEX (
            samplesheet
        )

        ch_fastq = BASECALL_SIMPLEX.out.fastq

        ch_versions = BASECALL_SIMPLEX.out.versions

        READS_FILTER (
            cfdna_samplesheet,
            ch_fastq,
            max_len,
            minqs
        )

        ch_fastq_processed = READS_FILTER.out.reads

        MAPPING(
            ch_fastq_processed,
            ref
        )

        if (params.m_bases) {
            ch_to_classify = MAPPING.out.bam
                .map { meta, _bam, _bai -> meta }
            ch_no_methylation_tags = channel.empty()
        } else {
            ch_no_methylation_tags = MAPPING.out.bam
                .map { meta, _bam, _bai -> meta }
            ch_to_classify = channel.empty()
        }

    }

    if (params.time_series) {

        ch_classify_to_split = MAPPING.out.bam
            .join(ch_to_classify)

        // Time series mode: Split BAMs into time intervals for temporal analysis
        SPLIT_BAMS_CLASSIFY(
            ch_classify_to_split,
            ref,
            bed
        )

        ch_bam_to_split = MAPPING.out.bam
            .join(ch_no_methylation_tags)
            .map { meta, bam, bai, _count ->
                tuple(meta, bam, bai)}

        SPLIT_BAMS_OTHER(
            ch_bam_to_split,
            ref,
            bed
        )

        // Use time series outputs for downstream variant calling
        ch_bam_for_calling = SPLIT_BAMS_CLASSIFY.out.bam
            .mix(SPLIT_BAMS_OTHER.out.bam)
            .map { meta, bam, bai ->
                tuple(id:meta.id, bam, bai)}
        ch_ref_for_calling = SPLIT_BAMS_CLASSIFY.out.ref
            .mix(SPLIT_BAMS_OTHER.out.ref)
            .map { meta, ref_id, ref_fa, ref_index ->
                tuple(id:meta.id, ref_id, ref_fa, ref_index)}
        ch_bed = SPLIT_BAMS_CLASSIFY.out.bed
            .mix(SPLIT_BAMS_OTHER.out.bed)
            .map { meta,bedfile,padding,low_fidelity ->
                tuple(id:meta.id,bedfile,padding,low_fidelity)}
        ch_bam_to_classify = SPLIT_BAMS_CLASSIFY.out.bam
            .map { meta, bam, bai ->
                tuple(id:meta.id, bam, bai)}

        ch_time_series = channel.of(params.time_points)

        ch_samplesheet_renamed = cfdna_samplesheet
            .join(tumor_type)
            .combine(ch_time_series)
            .map { meta, purity, filter, tumor, times ->
                def time_list = times.tokenize(',')
                tuple(meta, purity, filter, tumor, time_list)}
            .transpose()

        ch_cfdna_ss = ch_samplesheet_renamed
            .map { meta, purity, filter, _tumor, time ->
                def time_stripped = time.replace('"', '')
                def new_meta = meta.id + "_0h_${time_stripped}h"
                tuple(id:new_meta,purity, filter)}

        ch_tumor_ss = ch_samplesheet_renamed
            .map { meta, _purity, _filter, tumor, time ->
                def time_stripped = time.replace('"', '')
                def new_meta = meta.id + "_0h_${time_stripped}h"
                tuple(id:new_meta,tumor)}

        CRAMINO_STATS(ch_bam_for_calling.filter{meta,_bam,_bai -> !meta.id.contains('FULL')})

        ch_coverage_table = CRAMINO_STATS.out.stats
            .mix(MAPPING.out.coverage)

        if (params.include_full) {
            ch_cfdna_full = cfdna_samplesheet
                .map { meta, purity, filter ->
                    def new_meta = meta.id + '_FULL'
                tuple(id:new_meta,purity, filter)}
                .mix(ch_cfdna_ss)
            ch_bam_to_process = ch_bam_for_calling
                .filter{ meta, _bam, _bai -> meta.id.contains('FULL') }
            ch_tumor_type = tumor_type
                .map { meta, tumor ->
                    def new_meta = meta.id + '_FULL'
                    tuple(new_meta, tumor)}
                .mix(ch_tumor_ss)
        } else {
            ch_bam_to_process = ch_bam_for_calling
            ch_cfdna_full     = ch_cfdna_ss
            ch_tumor_type     = ch_tumor_ss
        }

    } else {
        // Standard mode: Use the full BAM directly for variant calling
        ch_bam_to_classify = ch_to_classify
                                .join(MAPPING.out.bam)
        ch_bam_for_calling = MAPPING.out.bam
        ch_bam_to_process  = MAPPING.out.bam
        ch_ref_for_calling = ref
        ch_bed             = bed
        ch_cfdna_full      = cfdna_samplesheet
        ch_coverage_table  = MAPPING.out.coverage
        ch_tumor_type      = tumor_type
    }

    if (params.m_bases || params.skip_basecalling || params.skip_mapping) {

        CLASSY(
            ch_bam_to_classify,
            ch_ref_for_calling,
            ch_tumor_type
        )

        CLASSIFIER_REPORT(
            CLASSY.out.plot,
            CLASSY.out.pred,
            ch_tumor_type
        )

        ch_sections = CLASSIFIER_REPORT.out.sections

        ch_versions = ch_versions
            .mix(CLASSY.out.versions)

    }

    // Analyze coverage separation between target and background regions
    COVERAGE_SEPARATE(
        ch_bam_for_calling,
        ch_bed,
        ch_bed_nopad,
        ch_ref_for_calling
    )
    ADAPTIVE_REPORT(
        COVERAGE_SEPARATE.out.coverage_tbl,
        COVERAGE_SEPARATE.out.coverage_plot
    )

    ch_sections = ch_sections
        .mix(ADAPTIVE_REPORT.out.sections)

    ch_versions = ch_versions
        .mix(COVERAGE_SEPARATE.out.versions)

    ch_bed_for_processing = ch_bed_nopad

    ch_vcf_subchrom = channel.empty()


    // Run snp calling only for higher coverage samples
    ch_coverage = MAPPING.out.coverage
        .map { meta, table ->
            def lines = table.readLines()
            def cov = lines[5].tokenize('\t')[1].toDouble()
            tuple(meta, cov)
        }
        .branch {
            meta, cov ->
            high: cov >= 10
                return meta
            low: cov < 10
                return meta
        }

    ch_high_cov_bam = ch_coverage.high
        .join(ch_bam_to_process)

    CLAIR3_CALLING (
        ch_high_cov_bam,
        ch_ref_for_calling,
        basecall_model,
        ch_bed_for_processing,
        vep_cache
    )

    CLAIRS_TO_CALLING (
        ch_high_cov_bam,
        ch_ref_for_calling,
        clairs_model,
        ch_bed_for_processing,
        vep_cache
    )

    ch_vcf_subchrom = ch_vcf_subchrom
        .mix(CLAIR3_CALLING.out.vcf_snpeff)

    SV_CALLING(
        ch_bam_to_process,
        ch_ref_for_calling,
        CLAIR3_CALLING.out.vcf_vep,
        vep_cache
    )

    CNV_CALLING (
        ch_bam_to_process,
        ch_vcf_subchrom,
        ch_ref_for_calling
    )

    ICHORCNA_CALLING (
        ch_bam_for_calling,
        ch_ref_for_calling,
        ch_cfdna_full,
        ichor_bin,
        mapq_wig
    )

    // Sub reports

    ch_binsize_qdnaseq = channel.of(params.qdnaseq_binsize)
        .map { value ->
        def meta = "qDNAseq"
        tuple(meta, value) }
    ch_binsize_subchrom = channel.of(params.subchrom_binsize)
        .map { value ->
        def meta = "Subchrom"
        tuple(meta, value) }
    ch_binsize_delly = channel.of(params.delly_bin_size)
        .map { value ->
        def meta = "Delly"
        tuple(meta, value) }

    ch_binsizes = ch_binsize_qdnaseq
        .mix(ch_binsize_subchrom)
        .mix(ch_binsize_delly)

    ch_snp_to_process = CLAIR3_CALLING.out.vcf_vep
        .mix(CLAIRS_TO_CALLING.out.vcf_vep)

    VARIANT_PROCESS (
        ch_bam_to_process,
        ch_bed_for_processing,
        SV_CALLING.out.vcf,
        SV_CALLING.out.stellerator,
        CNV_CALLING.out.qdnaseq_bed,
        CNV_CALLING.out.qdnaseq_segs,
        targets,
        CNV_CALLING.out.delly_cov,
        CNV_CALLING.out.delly_segs,
        ch_snp_to_process
    )

    ch_subchrom_plot = CNV_CALLING.out.subchrom_plot_wgs
    ch_subchrom_focal = CNV_CALLING.out.subchrom_gene_plot_wgs

    FIGENO_REPORT(
        VARIANT_PROCESS.out.circos_plot,
        VARIANT_PROCESS.out.panchr_plot,
        ch_binsizes,
        VARIANT_PROCESS.out.sv_plot,
        VARIANT_PROCESS.out.fusion_plot,
        VARIANT_PROCESS.out.targets_plot,
        VARIANT_PROCESS.out.sv_table,
        VARIANT_PROCESS.out.fusion_table,
        ch_subchrom_plot,
        ch_subchrom_focal,
        VARIANT_PROCESS.out.snp_table
    )

    CFDNA_REPORT(
        ch_cfdna_full,
        READS_FILTER.out.stats,
        ch_coverage_table,
        ICHORCNA_CALLING.out.ichorcna_plot
    )

        /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */


    ch_versions = ch_versions
        .mix(MAPPING.out.versions)
        .mix(CNV_CALLING.out.versions)
        .mix(ICHORCNA_CALLING.out.versions)
        .mix(CLAIR3_CALLING.out.versions)
        .mix(CLAIRS_TO_CALLING.out.versions)
        .mix(READS_FILTER.out.versions)
        .mix(SV_CALLING.out.versions)
        .mix(VARIANT_PROCESS.out.versions)
        .mix(FIGENO_REPORT.out.versions)
        .mix(CFDNA_REPORT.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT SECTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Collect sections from all analysis steps
    ch_sections = ch_sections
        .mix(CFDNA_REPORT.out.sections)
        .mix(FIGENO_REPORT.out.sections)

        // channel id containing only meta
    ch_params = ch_ref_for_calling
        .join(ch_bed)

    if (params.realtime) {
        ch_bam = MAPPING.out.bam
            .map { meta, bam, _bai ->
                tuple(meta,bam) }

        ONTIME_TIME_RANGE(ch_bam)

        // Automatically detect the time stamp

        ch_title = ONTIME_TIME_RANGE.out.txt
            .map { meta, file ->
                def lines     = file.readLines()
                def start_str = lines[0].split(':\\s+')[1].trim()
                def end_str   = lines[1].split(':\\s+')[1].trim()

                def start    = parseToInstant(start_str)
                def end      = parseToInstant(end_str)
                def duration = java.time.Duration.between(start, end)
                def hours    = duration.toHours()
                def minutes  = duration.toMinutesPart()

                def time_str = "${hours}h${minutes > 0 ? String.format('%02dm', minutes) : ''}"

                tuple(meta, "OncoSeq cfDNA Report in realtime — ${meta.id} (atfer ${time_str} of sequencing)")
            }
    } else {
        ch_title = ch_params
            .map { meta, refid, ref_fasta, ref_index, bed_file, padding, lowfid ->
            tuple(meta, "OncoSeq cfDNA Report — ${meta.id}")}
    }

    MIDNIGHT_REPORT(
        ch_params,
        ch_cfdna_full,
        ch_sections,
        ch_versions,
        ch_title
    )
}
