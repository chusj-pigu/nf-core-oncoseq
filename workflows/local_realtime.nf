//
// Import subworkflows for the adaptive sequencing pipeline
//

// Basecalling subworkflows
include { BASECALL_SIMPLEX     } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX   } from '../subworkflows/local/basecalling/basecall_multiplex'

// Core analysis subworkflows
include { MAPPING as MAPPING_HG  } from '../subworkflows/local/mapping/mapping'
include { MAPPING as MAPPING_T2T } from '../subworkflows/local/mapping/mapping'

// Tumor Classifiers
include { CLASSY                } from '../subworkflows/local/methylation_analysis/marlin.nf'
include { STURGEON              } from '../subworkflows/local/methylation_analysis/sturgeon.nf'

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
include { MIDNIGHT_REPORT } from '../subworkflows/local/report/final_report.nf'
include { CLASSIFIER_REPORT    } from '../subworkflows/local/report/methylation.nf'
include { FIGENO_REPORT  } from '../subworkflows/local/report/variants.nf'
include { ADAPTIVE_REPORT } from '../subworkflows/local/report/adaptive.nf'

workflow LOCAL_REALTIME {

    take:
    samplesheet             // channel: samplesheet read in from --input
    demux_samplesheet       // channel: demux samplesheet read in from --demux_samplesheet
    ref                     // channel: reference for mapping, either empty if skipping mapping, or a path
    tumor_type              // channel: samplesheet read in from --input, contains only tumor type
    ref_t2t                 // channel: Path to T2T reference from params.ref_t2t
    basecall_model          // channel: model for basecalling
    ch_clin_database        // channel: clinical database for variant annotation
    bed                     // channel: bed file used for adaptive sampling regions
    targets                 // channel : list of genes with their position to represent in Figeno

    main:

    ch_versions = channel.empty()
    ch_sections = channel.empty()

    // Branch by tumor type
    ch_tumor_type = tumor_type
        .branch { meta, tumor ->
            leukemia: tumor == "leukemia"
            cns: tumor == "cns"
            other: tumor == "other"
        }

    if (params.skip_basecalling || params.skip_mapping) {

        ch_mapping_t2t = ch_tumor_type.cns
            .join(samplesheet)
            .map { meta, _tumor, input ->
                tuple(meta, input)}

        MAPPING_HG(
            samplesheet,
            ref
        )

        MAPPING_T2T(
            ch_mapping_t2t,
            ref_t2t
        )

    } else {
        if (params.demux) {

            // Perform multiplex basecalling with demultiplexing
            BASECALL_MULTIPLEX (
                samplesheet,
                demux_samplesheet,
                ref
            )

            ch_mapping_t2t = ch_tumor_type.cns
                .join(BASECALL_MULTIPLEX.out.fastq)
                .map { meta, _tumor, input ->
                    tuple(meta, input)}

            // Map basecalled reads to reference
            MAPPING_HG (
                BASECALL_MULTIPLEX.out.fastq,
                BASECALL_MULTIPLEX.out.ref
            )

            // Map basecalled reads from CNS tumor to t2t
            MAPPING_T2T (
                ch_mapping_t2t,
                ref_t2t
            )

            ch_versions = ch_versions
                .mix(BASECALL_MULTIPLEX.out.versions)

        } else {
            // Sub-branch 2b: Simplex basecalling (single sample per flow cell)

            // Perform simplex basecalling
            BASECALL_SIMPLEX (
                samplesheet
            )

            ch_mapping_t2t = ch_tumor_type.cns
                .join(BASECALL_SIMPLEX.out.fastq)
                .map { meta, _tumor, input ->
                    tuple(meta, input)}

            // Map basecalled reads to reference
            MAPPING_HG (
                BASECALL_SIMPLEX.out.fastq,
                ref
            )

            // Map basecalled reads from CNS tumor to t2t
            MAPPING_T2T (
                ch_mapping_t2t,
                ref_t2t
            )

            ch_versions = ch_versions
                .mix(BASECALL_SIMPLEX.out.versions)
        }
    }

    if (params.realtime < 6) {                 // Before 6h of realtime sequencing, include CNV calling with QDNAseq, SV calling and Marlin

        // Run Marlin on leukemia samples only:
        ch_marlin_bam = MAPPING_HG.out.bam
            .join(ch_tumor_type.leukemia)
            .map { meta, bam, bai, _tumor ->
                tuple(meta, bam, bai)}

        CLASSY(
            ch_marlin_bam,
            ref
        )

        STURGEON(
            MAPPING_T2T.out.bam
        )

        CNV_CALLING(
            MAPPING_HG.out.bam,
            ref
        )

        SV_UNPHASED(
            MAPPING_HG.out.bam,
            ref
        )

        // Filter variants to visualize :
        VARIANT_PROCESS (
            MAPPING_HG.out.bam,
            SV_UNPHASED.out.vcf,
            CNV_CALLING.out.qdnaseq_bed,
            CNV_CALLING.out.qdnaseq_segs,
            targets,
            CNV_CALLING.out.delly_cov,
            CNV_CALLING.out.delly_segs
        )

        COVERAGE_SEPARATE(
            MAPPING_HG.out.bam,
            bed,
            ref
        )

        // Placeholders for report
        ch_subchrom_focal = channel.empty()
        ch_subchrom_plot = channel.empty()

        ch_versions = ch_versions
            .mix(MAPPING_T2T.out.versions)
            .mix(CLASSY.out.versions)
            .mix(STURGEON.out.versions)

        ch_classifiers_plots = STURGEON.out.plot
            .mix(CLASSY.out.plot)

        ch_classifiers_pred = STURGEON.out.pred
            .mix(CLASSY.out.pred)

        CLASSIFIER_REPORT(
            ch_classifiers_plots,
            ch_classifiers_pred
        )

        ch_sections = CLASSIFIER_REPORT.out.sections

    } else if (params.realtime >=6 & params.realtime < 72 ) {

        CNV_CALLING(
            MAPPING_HG.out.bam,
            ref
        )

        SV_UNPHASED(
            MAPPING_HG.out.bam,
            ref
        )

        // Filter variants to visualize :
        VARIANT_PROCESS (
            MAPPING_HG.out.bam,
            SV_UNPHASED.out.vcf,
            CNV_CALLING.out.qdnaseq_bed,
            CNV_CALLING.out.qdnaseq_segs,
            targets,
            CNV_CALLING.out.delly_cov,
            CNV_CALLING.out.delly_segs
        )

        COVERAGE_SEPARATE(
            MAPPING_HG.out.bam,
            bed,
            ref
        )
        // Germline variant calling using Clair3 (always uses original mapping output)
        CLAIR3_CALLING (
            MAPPING_HG.out.bam,
            ref,
            basecall_model,
            ch_clin_database,
            bed
        )

        // Placeholders for report
        ch_subchrom_focal = channel.empty()
        ch_subchrom_plot = channel.empty()

        ch_versions = ch_versions
            .mix(CLAIR3_CALLING.out.versions)

    } else if (params.realtime == 72) {
        CNV_CALLING(
            MAPPING_HG.out.bam,
            ref
        )

        SV_UNPHASED(
            MAPPING_HG.out.bam,
            ref
        )

        // Filter variants to visualize :
        VARIANT_PROCESS (
            MAPPING_HG.out.bam,
            SV_UNPHASED.out.vcf,
            CNV_CALLING.out.qdnaseq_bed,
            CNV_CALLING.out.qdnaseq_segs,
            targets,
            CNV_CALLING.out.delly_cov,
            CNV_CALLING.out.delly_segs
        )

        COVERAGE_SEPARATE(
            MAPPING_HG.out.bam,
            bed,
            ref
        )
        // Germline variant calling using Clair3 (always uses original mapping output)
        CLAIR3_CALLING (
            MAPPING_HG.out.bam,
            ref,
            basecall_model,
            ch_clin_database,
            bed
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
            MAPPING_HG.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf,
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
        .mix(MAPPING_HG.out.versions)
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
        ch_subchrom_focal
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

    ch_mode = channel.of("Adaptive Sampling atfer ${params.realtime}h sequencing")

     // channel id containing only meta
    ch_id = MAPPING_HG.out.bam
        .map { meta, _bam, _bai ->
        meta }


    MIDNIGHT_REPORT(
        ch_id,
        ch_sections,
        ch_versions,
        ch_mode
    )

}
