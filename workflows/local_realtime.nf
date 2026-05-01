//
// Import subworkflows for the adaptive sequencing pipeline
//

// Basecalling subworkflows
include { BASECALL_SIMPLEX     } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX   } from '../subworkflows/local/basecalling/basecall_multiplex'

// Core analysis subworkflows
include { MAPPING              } from '../subworkflows/local/mapping/mapping'

// Tumor Classifiers
include { CLASSY                } from '../subworkflows/local/methylation_analysis/classy.nf'
include { SAMTOOLS_COUNT_READS                 } from '../modules/local/samtools/main.nf'

// Variant calling subworkflows
include { CLAIR3_CALLING                        } from '../subworkflows/local/variant_calling/clair3_calling.nf'
include { SV_CALLING as SV_UNPHASED             } from  '../subworkflows/local/variant_calling/sv_calling.nf'
include { CNV_CALLING                           } from  '../subworkflows/local/variant_calling/cnv_calling.nf'
include { SUBCHROM_CALL                         } from  '../subworkflows/local/variant_calling/subchrom_call.nf'

// Variant processing and visualization subworkflow
include { VARIANT_PROCESS                       } from  '../subworkflows/local/variant_calling/variant_process.nf'

// Adaptive-specific subworkflows
include { COVERAGE_SEPARATE } from '../subworkflows/local/adaptive_specific/coverage_separate'

//
include { SUBCHROM_PANEL_BIN    } from '../modules/local/subchrom/main.nf'
include { REMOVE_PADDING        } from '../modules/local/adaptive_specific/main.nf'
include { modifyMetaId          } from '../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'

// Reporting
include { ONTIME_TIME_RANGE } from '../modules/local/ontime/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TIME_FASTQ } from '../modules/local/samtools/main.nf'
include { MIDNIGHT_REPORT   } from '../subworkflows/local/report/final_report.nf'
include { CLASSIFIER_REPORT } from '../subworkflows/local/report/methylation.nf'
include { FIGENO_REPORT     } from '../subworkflows/local/report/variants.nf'
include { ADAPTIVE_REPORT   } from '../subworkflows/local/report/adaptive.nf'

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

workflow LOCAL_REALTIME {

    take:
    samplesheet             // channel: samplesheet read in from --input
    demux_samplesheet       // channel: demux samplesheet read in from --demux_samplesheet
    ref                     // channel: reference for mapping, either empty if skipping mapping, or a path
    basecall_model          // channel: model for basecalling
    bed                     // channel: bed file used for adaptive sampling regions
    targets                 // channel : list of genes with their position to represent in Figeno
    vep_cache

    main:

    ch_versions = channel.empty()
    ch_sections = channel.empty()

    if (params.skip_basecalling) {

        MAPPING(
            samplesheet,
            ref
        )

        ch_bam = MAPPING.out.bam
            .map { meta, bam, _bai ->
                tuple(meta,bam) }

        ONTIME_TIME_RANGE(ch_bam)

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
            }

        ch_bam_methylation_counts.none
            .subscribe { meta, count ->
                log.warn "Tumor classification will be skipped for ${meta.id} -- no methylation tags found in bam"
            }
        ch_classy_in = ch_bam_methylation_counts.pos
            .join(MAPPING.out.bam)
            .view()

    } else if (params.skip_mapping) {

        MAPPING(
            samplesheet,
            ref
        )

        // We have to convert to fastq when input is from MinKnow :
        ch_bam = MAPPING.out.bam
            .map { meta, bam, _bai ->
                tuple(meta,bam) }

        SAMTOOLS_TIME_FASTQ(ch_bam)

        ONTIME_TIME_RANGE(SAMTOOLS_TIME_FASTQ.out.fq)

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
            }

        ch_bam_methylation_counts.none
            .subscribe { meta, count ->
                log.warn "Tumor classification will be skipped for ${meta.id} — no methylation tags found in bam"
            }
        ch_classy_in = ch_bam_methylation_counts.pos
            .join(MAPPING.out.bam)

    } else {

        if (params.demux) {

            // Perform multiplex basecalling with demultiplexing
            BASECALL_MULTIPLEX (
                samplesheet,
                demux_samplesheet,
                ref
            )

            // Map basecalled reads to reference
            MAPPING (
                BASECALL_MULTIPLEX.out.fastq,
                BASECALL_MULTIPLEX.out.ref
            )

            ch_versions = ch_versions
                .mix(BASECALL_MULTIPLEX.out.versions)

        } else {
            // Sub-branch 2b: Simplex basecalling (single sample per flow cell)

            // Perform simplex basecalling
            BASECALL_SIMPLEX (
                samplesheet
            )

            // Map basecalled reads to reference
            MAPPING (
                BASECALL_SIMPLEX.out.fastq,
                ref
            )

            ch_versions = ch_versions
                .mix(BASECALL_SIMPLEX.out.versions)
        }

        ch_bam = MAPPING.out.bam
            .map { meta, bam, _bai ->
                tuple(meta,bam) }

        ONTIME_TIME_RANGE(ch_bam)

        ch_classy_in = MAPPING.out.bam
    }

    if (params.realtime < 6) {                 // Before 6h of realtime sequencing, include CNV calling with QDNAseq, SV calling and Marlin

        if (params.m_bases || params.skip_mapping || params.skip_basecalling ) {

            CLASSY(
                ch_classy_in,
                ref
            )

            CLASSIFIER_REPORT(
                CLASSY.out.plot,
                CLASSY.out.pred
            )

            ch_sections = CLASSIFIER_REPORT.out.sections
        }

        CNV_CALLING(
            MAPPING.out.bam,
            ref
        )

        SV_UNPHASED(
            MAPPING.out.bam,
            ref
        )

        ch_clair3_out = channel.empty()

        // Filter variants to visualize :
        VARIANT_PROCESS (
            MAPPING.out.bam,
            SV_UNPHASED.out.vcf,
            CNV_CALLING.out.qdnaseq_bed,
            CNV_CALLING.out.qdnaseq_segs,
            targets,
            CNV_CALLING.out.delly_cov,
            CNV_CALLING.out.delly_segs,
            ch_clair3_out
        )

        COVERAGE_SEPARATE(
            MAPPING.out.bam,
            bed,
            ref
        )

        // Placeholders for report
        ch_subchrom_focal = channel.empty()
        ch_subchrom_plot = channel.empty()

        ch_versions = ch_versions
            .mix(CLASSY.out.versions)

    } else if (params.realtime >=6 & params.realtime < 72 ) {

        CNV_CALLING(
            MAPPING.out.bam,
            ref
        )

        SV_UNPHASED(
            MAPPING.out.bam,
            ref
        )

        COVERAGE_SEPARATE(
            MAPPING.out.bam,
            bed,
            ref
        )
        // Germline variant calling using Clair3 (always uses original mapping output)
        CLAIR3_CALLING (
            MAPPING.out.bam,
            ref,
            basecall_model,
            COVERAGE_SEPARATE.out.split_bed,
            vep_cache
        )

        // Filter variants to visualize :
        VARIANT_PROCESS (
            MAPPING.out.bam,
            SV_UNPHASED.out.vcf,
            CNV_CALLING.out.qdnaseq_bed,
            CNV_CALLING.out.qdnaseq_segs,
            targets,
            CNV_CALLING.out.delly_cov,
            CNV_CALLING.out.delly_segs,
            CLAIR3_CALLING.out.vcf_vep
        )

        // Placeholders for report
        ch_subchrom_focal = channel.empty()
        ch_subchrom_plot = channel.empty()

        ch_versions = ch_versions
            .mix(CLAIR3_CALLING.out.versions)

    } else if (params.realtime == 72) {
        CNV_CALLING(
            MAPPING.out.bam,
            ref
        )

        SV_UNPHASED(
            MAPPING.out.bam,
            ref
        )

        COVERAGE_SEPARATE(
            MAPPING.out.bam,
            bed,
            ref
        )
        // Germline variant calling using Clair3 (always uses original mapping output)
        CLAIR3_CALLING (
            MAPPING.out.bam,
            ref,
            basecall_model,
            COVERAGE_SEPARATE.out.split_bed,
            vep_cache
        )

        // Filter variants to visualize :
        VARIANT_PROCESS (
            MAPPING.out.bam,
            SV_UNPHASED.out.vcf,
            CNV_CALLING.out.qdnaseq_bed,
            CNV_CALLING.out.qdnaseq_segs,
            targets,
            CNV_CALLING.out.delly_cov,
            CNV_CALLING.out.delly_segs,
            CLAIR3_CALLING.out.vcf_vep
        )

        ch_subchrom_panelbin_in = COVERAGE_SEPARATE.out.split_bed
            .map {
                meta, panelbed ->
                tuple(meta, panelbed)
            }
            .join(ref)
            .map {
                meta, panelbed, refid, _ref, _ref_fai ->
                tuple(meta, panelbed, refid, params.subchrom_binsize )
            }

        ch_panel_bin = SUBCHROM_PANEL_BIN(ch_subchrom_panelbin_in).subchrom_panelbin_bed

        SUBCHROM_CALL (
            MAPPING.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf_snpeff,
            ch_panel_bin
        )

        ch_subchrom_plot = SUBCHROM_CALL.out.subchrom_plot_wgs
            //.mix(SUBCHROM_CALL.out.subchrom_plot_panel)

        ch_subchrom_focal = SUBCHROM_CALL.out.subchrom_gene_plot_wgs
            //.mix(SUBCHROM_CALL.out.subchrom_gene_plot_panel)

        ch_versions = ch_versions
            .mix(SUBCHROM_CALL.out.versions)
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT VERSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_versions = ch_versions
        .mix(MAPPING.out.versions)
        .mix(COVERAGE_SEPARATE.out.versions)
        .mix(CNV_CALLING.out.versions)
        .mix(SV_UNPHASED.out.versions)
        .mix(VARIANT_PROCESS.out.versions)

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPILE SECTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ADAPTIVE_REPORT(
        COVERAGE_SEPARATE.out.coverage_tbl,
        COVERAGE_SEPARATE.out.coverage_plot
    )

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

    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SECTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Collect sections from all analysis steps
    ch_sections = ch_sections
        .mix(ADAPTIVE_REPORT.out.sections)
        .mix(FIGENO_REPORT.out.sections)


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

            tuple(meta, "OncoSeq Adaptive Sampling Report in realtime — ${meta.id} (atfer ${time_str} of sequencing)")
        }

     // channel id containing only meta
    ch_id = MAPPING.out.bam
        .map { meta, _bam, _bai ->
        meta }


    MIDNIGHT_REPORT(
        ch_id,
        ch_sections,
        ch_versions,
        ch_title
    )

}
