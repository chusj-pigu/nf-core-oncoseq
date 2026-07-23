//
// Import subworkflows for the adaptive sequencing pipeline
//

// Basecalling subworkflows
include { BASECALL_SIMPLEX   } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'

// Core analysis subworkflows
include { MAPPING            } from '../subworkflows/local/mapping/mapping'

// Variant calling subworkflows
include { CLAIRS_TO_CALLING                     } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'
include { CLAIR3_CALLING                        } from '../subworkflows/local/variant_calling/clair3_calling.nf'
include { PHASING_VARIANTS as PHASING_SOMATIC   } from  '../subworkflows/local/variant_calling/phasing.nf'
include { PHASING_VARIANTS as PHASING_GERMLINE  } from  '../subworkflows/local/variant_calling/phasing.nf'
include { SV_CALLING                            } from  '../subworkflows/local/variant_calling/sv_calling.nf'
include { CNV_CALLING                           } from  '../subworkflows/local/variant_calling/cnv_calling.nf'

// Variant processing and visualization subworkflow
include { VARIANT_PROCESS                       } from  '../subworkflows/local/variant_calling/variant_process.nf'

// Adaptive-specific subworkflows
include { COVERAGE_SEPARATE } from '../subworkflows/local/adaptive_specific/coverage_separate'

// Tumor Classifiers
include { SAMTOOLS_COUNT_READS  } from '../modules/local/samtools/main.nf'
include { CLASSY                } from '../subworkflows/local/methylation_analysis/classy.nf'

// Time series evaluation subworkflows
include { SPLIT_BAMS_TIME                      } from '../subworkflows/local/time_series_evaluation/split_bams.nf'
include { SPLIT_BAMS_TIME_FASTQ                } from '../subworkflows/local/time_series_evaluation/split_bams_fastq.nf'
include { SUBSAMPLE_TIME                       } from '../subworkflows/local/read_processing/subsample_time.nf'

// Reporting
include { MIDNIGHT_REPORT   } from '../subworkflows/local/report/final_report.nf'
include { FIGENO_REPORT     } from '../subworkflows/local/report/variants.nf'
include { ADAPTIVE_REPORT   } from '../subworkflows/local/report/adaptive.nf'
include { CLASSIFIER_REPORT } from '../subworkflows/local/report/methylation.nf'
// Utility functions
include { modifyMetaId          } from '../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'

include { REMOVE_PADDING        } from '../modules/local/adaptive_specific/main.nf'

//
// WORKFLOW: Adaptive sequencing analysis pipeline
//
// This workflow handles two main scenarios:
// 1. Skip basecalling: Start from pre-basecalled FASTQ files
// 2. Full pipeline: Perform basecalling (simplex or multiplex) followed by analysis
//
workflow ADAPTIVE_WGS {

    take:
    samplesheet             // channel: samplesheet read in from --input
    demux_samplesheet       // channel: demux samplesheet read in from --demux_samplesheet
    ref                     // channel: reference for mapping, either empty if skipping mapping, or a path
    clairs_model            // channel: model for ClairS variant calling
    basecall_model          // channel: model for basecalling
    bed                     // channel: bed file used for adaptive sampling regions
    targets                 // channel : list of genes with their position to represent in Figeno
    vep_cache

    main:


    //
    // WORKFLOW: Run pipeline
    //

    ch_versions = channel.empty()

    // Branch 1: Skip basecalling - start from pre-basecalled FASTQ files
    if (params.skip_basecalling || params.skip_mapping) {

        // All samples need to be mapped to hg19 or hg38
        MAPPING(
            samplesheet,
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
            }

        ch_bam_methylation_counts.none
            .subscribe { meta, count ->
                log.warn "Tumor classification will be skipped for ${meta.id} -- no methylation tags found in bam"
            }
        ch_to_classify = ch_bam_methylation_counts.pos
            .join(MAPPING.out.bam)

    } else if (params.demux) {

        // Perform multiplex basecalling with demultiplexing
        BASECALL_MULTIPLEX (
            samplesheet,
            demux_samplesheet
        )

        // Map basecalled reads to reference
        MAPPING (
            BASECALL_MULTIPLEX.out.fastq,
            ref
        )

        ch_fastq = BASECALL_MULTIPLEX.out.fastq

        ch_versions = ch_versions
            .mix(BASECALL_MULTIPLEX.out.versions)

        ch_to_classify = MAPPING.out.bam
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

        ch_fastq = BASECALL_SIMPLEX.out.fastq

        ch_to_classify = MAPPING.out.bam

    }

    if (params.time_series) {
        // Time series mode: Split BAMs into time intervals for temporal analysis
        SPLIT_BAMS_TIME(
            MAPPING.out.bam,
            ref,
            bed
        )

        // Use time series outputs for downstream variant calling
        ch_bam_for_calling = SPLIT_BAMS_TIME.out.bam
            .map { meta, bam, bai ->
                tuple(id:meta.id, bam, bai)}
        ch_ref_for_calling = SPLIT_BAMS_TIME.out.ref
            .map { meta, ref_id, ref_fa, ref_index ->
                tuple(id:meta.id, ref_id, ref_fa, ref_index)}
        ch_bed = SPLIT_BAMS_TIME.out.bed
            .map { meta,bedfile,padding,low_fidelity ->
                tuple(id:meta.id,bedfile,padding,low_fidelity)}

        ch_versions = ch_versions
                .mix(SPLIT_BAMS_TIME.out.versions)

        if (params.m_bases) {
            ch_bam_1h = ch_bam_for_calling
                .filter { meta, _bam, _index ->
                meta.endsWith('_0h_1h') }

            ch_bam_classy = ch_bam_1h

            CLASSY(
                ch_bam_classy,
                ch_ref_for_calling
            )

            CLASSIFIER_REPORT(
                CLASSY.out.plot,
                CLASSY.out.pred
            )

            ch_classy_section = CLASSIFIER_REPORT.out.sections
        }

    } else {
        // Standard mode: Use the full BAM directly for variant calling
        ch_bam_for_calling = MAPPING.out.bam
        ch_bam_full = MAPPING.out.bam
        ch_ref_for_calling = ref
        ch_bed = bed

        // Downsample to 1h to run methylation classification

        if (params.m_bases || params.skip_basecalling || params.skip_mapping) {

            ch_in_subsample_classy = ch_to_classify
                .map { meta, bam, index ->
                tuple(meta, bam, index, 0, 1)
                }

            SUBSAMPLE_TIME(
                ch_in_subsample_classy
            )

            ch_in_classy = SUBSAMPLE_TIME.out.bam

            CLASSY(
                ch_in_classy,
                ch_ref_for_calling
            )

            ch_versions = ch_versions
                .mix(SUBSAMPLE_TIME.out.versions)

            CLASSIFIER_REPORT(
                CLASSY.out.plot,
                CLASSY.out.pred
            )

            ch_classy_section = CLASSIFIER_REPORT.out.sections
        } else {
            ch_classy_section = channel.empty()
        }
    }

    // Analyze coverage separation between target and background regions
    COVERAGE_SEPARATE(
        ch_bam_for_calling,
        ch_bed,
        ch_ref_for_calling
    )

    // Somatic variant calling using ClairS
    CLAIRS_TO_CALLING (
        ch_bam_for_calling,
        ch_ref_for_calling,
        clairs_model,
        COVERAGE_SEPARATE.out.split_bed,
        vep_cache
    )

    // Germline variant calling using Clair3 (always uses original mapping output)
    CLAIR3_CALLING (
        ch_bam_for_calling,
        ch_ref_for_calling,
        basecall_model,
        COVERAGE_SEPARATE.out.split_bed,
        vep_cache
    )

    // // Phase somatic variants (uses original mapping output)
    PHASING_SOMATIC (
        ch_bam_for_calling,
        ch_ref_for_calling,
        CLAIRS_TO_CALLING.out.vcf_snpeff
    )

    // Phase germline variants (can use time series BAM if enabled)
    PHASING_GERMLINE (
        ch_bam_for_calling,
        ch_ref_for_calling,
        CLAIR3_CALLING.out.vcf_snpeff
    )

    // Structural variant calling using phased BAM
    SV_CALLING (
        PHASING_GERMLINE.out.haptag_bam
            .map { meta, bamfile, bai ->
            // Restore original sample ID for output naming
                def meta_restore = modifyMetaId(meta, 'replace', '_germline_snp_snpeff_phased', '', '')
            tuple(meta_restore, bamfile, bai)
            },
        ch_ref_for_calling,
        PHASING_GERMLINE.out.phased_vcf
            .map { meta, vcf ->
                def meta_restore = modifyMetaId(meta, 'replace', '_germline_snp_snpeff_phased', '', '')
            tuple(meta_restore, vcf)
            },
        vep_cache
    )

    // Copy number variant calling
    CNV_CALLING(
        ch_bam_for_calling,
        CLAIR3_CALLING.out.vcf_snpeff,
        ch_ref_for_calling
    )

    ch_snp_to_process = CLAIR3_CALLING.out.vcf_vep
        .mix(CLAIRS_TO_CALLING.out.vcf_vep)

    // Filter variants to visualize :
    VARIANT_PROCESS (
        ch_bam_for_calling,
        COVERAGE_SEPARATE.out.split_bed,
        SV_CALLING.out.vcf,
        SV_CALLING.out.stellerator,
        CNV_CALLING.out.qdnaseq_bed,
        CNV_CALLING.out.qdnaseq_segs,
        targets,
        CNV_CALLING.out.delly_cov,
        CNV_CALLING.out.delly_segs,
        ch_snp_to_process
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT VERSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


    ch_versions = ch_versions
        .mix(MAPPING.out.versions)
        .mix(CLAIRS_TO_CALLING.out.versions)
        .mix(CLAIR3_CALLING.out.versions)
        .mix(PHASING_SOMATIC.out.versions)
        .mix(PHASING_GERMLINE.out.versions)
        .mix(SV_CALLING.out.versions)
        .mix(CNV_CALLING.out.versions)
        .mix(VARIANT_PROCESS.out.versions)

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPILE SECTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ADAPTIVE_REPORT(
        COVERAGE_SEPARATE.out.coverage_tbl,
        COVERAGE_SEPARATE.out.coverage_plot
    )

    ch_subchrom_plot = CNV_CALLING.out.subchrom_plot_wgs

    ch_subchrom_focal = CNV_CALLING.out.subchrom_gene_plot_wgs

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
        SV_CALLING.out.empty_calls,
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
    ch_sections = ADAPTIVE_REPORT.out.sections
        .mix(FIGENO_REPORT.out.sections)
        .mix(ch_classy_section)

    // channel containing parameters used:

    ch_cfdna = channel.empty()

    ch_params = ref
        .join(bed)

    ch_title = ch_params
        .map { meta, _refid, _ref_fasta, _ref_index, _bed_file, _padding, _lowfid ->
        tuple(meta, "OncoSeq Workflow Report — ${meta.id}")}

    MIDNIGHT_REPORT(
        ch_params,
        ch_cfdna,
        ch_sections,
        ch_versions,
        ch_title
    )
}

