//
// Import subworkflows for the adaptive sequencing pipeline
//

// Basecalling subworkflows
include { BASECALL_SIMPLEX  } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'

// Core analysis subworkflows
include { MAPPING           } from '../subworkflows/local/mapping/mapping'

// Variant calling subworkflows
include { CLAIRS_TO_CALLING } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'
include { CLAIR3_CALLING } from '../subworkflows/local/variant_calling/clair3_calling.nf'
include { PHASING_VARIANTS as PHASING_SOMATIC  } from  '../subworkflows/local/variant_calling/phasing.nf'
include { PHASING_VARIANTS as PHASING_GERMLINE  } from  '../subworkflows/local/variant_calling/phasing.nf'
include { SV_CALLING        } from  '../subworkflows/local/variant_calling/sv_calling.nf'
include { CNV_CALLING       } from  '../subworkflows/local/variant_calling/cnv_calling.nf'

// Adaptive-specific subworkflows
include { COVERAGE_SEPARATE } from '../subworkflows/local/adaptive_specific/coverage_separate'

// Time series evaluation subworkflows
include { SPLIT_BAMS_TIME   } from '../subworkflows/local/time_series_evaluation/split_bams.nf'
include { SPLIT_BAMS_TIME_FASTQ   } from '../subworkflows/local/time_series_evaluation/split_bams_fastq.nf'
include { modifyMetaId    } from '../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'

//
// WORKFLOW: Adaptive sequencing analysis pipeline
//
// This workflow handles two main scenarios:
// 1. Skip basecalling: Start from pre-basecalled FASTQ files
// 2. Full pipeline: Perform basecalling (simplex or multiplex) followed by analysis
//
workflow ADAPTIVE {

    take:
    samplesheet             // channel: samplesheet read in from --input
    demux_samplesheet       // channel: demux samplesheet read in from --demux_samplesheet
    ref                     // channel: reference for mapping, either empty if skipping mapping, or a path
    clairs_model            // channel: model for ClairS variant calling
    basecall_model          // channel: model for basecalling
    ch_clin_database        // channel: clinical database for variant annotation
    bed                     // channel: bed file used for adaptive sampling regions

    main:


    //
    // WORKFLOW: Run pipeline
    //

    // Branch 1: Skip basecalling - start from pre-basecalled FASTQ files
    if (params.skip_basecalling) {
        // Map FASTQ reads to reference genome
        MAPPING(
            samplesheet,
            ref
            )

        // Somatic variant calling using ClairS
        CLAIRS_TO_CALLING (
            MAPPING.out.bam,
            ref,
            clairs_model,
            ch_clin_database
        )

        // Germline variant calling using Clair3
        CLAIR3_CALLING (
            MAPPING.out.bam,
            ref,
            basecall_model,
            ch_clin_database
        )

        // Phase somatic variants
        PHASING_SOMATIC (
            MAPPING.out.bam,
            ref,
            CLAIRS_TO_CALLING.out.vcf
        )

        // Phase germline variants
        PHASING_GERMLINE (
            MAPPING.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf
        )

        // Structural variant calling using phased BAM
        SV_CALLING (
                PHASING_GERMLINE.out.haptag_bam
                .map { meta, bamfile, bai ->
                    // Restore original sample ID for output naming
                    def meta_restore = modifyMetaId(meta, 'replace', '_somatic_snp_phased', '', '')
                    meta_restore = modifyMetaId(meta_restore, 'replace', '_germline_snp_phased', '', '')
                    tuple(meta_restore, bamfile, bai)
                },
                ref
        )

        // Copy number variant calling
        CNV_CALLING (
            MAPPING.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf
        )

    } else {
        // Branch 2: Full pipeline - perform basecalling first

        // Sub-branch 2a: Multiplex basecalling (multiple samples per flow cell)
        if (params.demux != null) {

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

        }

        // Conditional processing for time series analysis
        if (params.time_series) {
            // Time series mode: Split BAMs into time intervals for temporal analysis
            SPLIT_BAMS_TIME(
                MAPPING.out.bam,
                ref,
                bed
            )

            // Use time series outputs for downstream variant calling
            ch_bam_for_calling = SPLIT_BAMS_TIME.out.bam
            ch_ref_for_calling = SPLIT_BAMS_TIME.out.ref
            ch_bed = SPLIT_BAMS_TIME.out.bed
        } else {
            // Standard mode: Use the full BAM directly for variant calling
            ch_bam_for_calling = MAPPING.out.bam
            ch_ref_for_calling = ref
            ch_bed = bed
        }

        // Analyze coverage separation between target and background regions
        COVERAGE_SEPARATE(
            ch_bam_for_calling,
            ch_bed
        )

        // Somatic variant calling using ClairS
        CLAIRS_TO_CALLING (
            ch_bam_for_calling,
            ch_ref_for_calling,
            clairs_model,
            ch_clin_database
        )

        // Germline variant calling using Clair3 (always uses original mapping output)
        CLAIR3_CALLING (
            ch_bam_for_calling,
            ch_ref_for_calling,
            basecall_model,
            ch_clin_database
        )

        // // Phase somatic variants (uses original mapping output)
        PHASING_SOMATIC (
            ch_bam_for_calling,
            ch_ref_for_calling,
            CLAIRS_TO_CALLING.out.vcf
        )

        // Phase germline variants (can use time series BAM if enabled)
        PHASING_GERMLINE (
            ch_bam_for_calling,
            ch_ref_for_calling,
            CLAIR3_CALLING.out.vcf
        )

        // Structural variant calling using phased BAM
        // TODO: find
        SV_CALLING (
            PHASING_GERMLINE.out.haptag_bam
                .map { meta, bamfile, bai ->
                    // Restore original sample ID for output naming
                    def meta_restore = modifyMetaId(meta, 'replace', '_somatic_snp_phased', '', '')
                    meta_restore = modifyMetaId(meta_restore, 'replace', '_germline_snp_phased', '', '')
                    tuple(meta_restore, bamfile, bai)
                },
            ch_ref_for_calling
        )

        // Copy number variant calling
        CNV_CALLING(
            ch_bam_for_calling,
            ch_ref_for_calling,
            CLAIR3_CALLING.out.vcf
        )
    }
}
