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
include { MARLIN                } from '../subworkflows/local/methylation_analysis/marlin.nf'
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
        }
    }

    if (params.realtime < 6) {                 // Before 6h of realtime sequencing, include CNV calling with QDNAseq, SV calling and Marlin

        MARLIN(
            MAPPING_HG.out.bam,
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
            targets
        )

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
            targets
        )

        COVERAGE_SEPARATE(
            MAPPING_HG.out.bam,
            bed
        )
        // Germline variant calling using Clair3 (always uses original mapping output)
        CLAIR3_CALLING (
            MAPPING_HG.out.bam,
            ref,
            basecall_model,
            ch_clin_database
        )
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
            targets
        )

        COVERAGE_SEPARATE(
            MAPPING_HG.out.bam,
            bed
        )
        // Germline variant calling using Clair3 (always uses original mapping output)
        CLAIR3_CALLING (
            MAPPING_HG.out.bam,
            ref,
            basecall_model,
            ch_clin_database
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
    }
}
