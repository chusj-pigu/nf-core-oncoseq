//
// Import subworkflows for the adaptive sequencing pipeline
//

// Basecalling subworkflows
include { BASECALL_SIMPLEX     } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX   } from '../subworkflows/local/basecalling/basecall_multiplex'

// Core analysis subworkflows
include { MAPPING           } from '../subworkflows/local/mapping/mapping'

// Variant calling subworkflows
include { CLAIR3_CALLING                        } from '../subworkflows/local/variant_calling/clair3_calling.nf'
include { SV_CALLING as SV_UNPHASED             } from  '../subworkflows/local/variant_calling/sv_calling.nf'
include { CNV_CALLING                           } from  '../subworkflows/local/variant_calling/cnv_calling.nf'
include { SUBCHROM_CALL                         } from  '../subworkflows/local/variant_calling/subchrom_call.nf'

// Adaptive-specific subworkflows
include { COVERAGE_SEPARATE } from '../subworkflows/local/adaptive_specific/coverage_separate'

//
include { SUBCHROM_PANEL_BIN    } from '../modules/local/subchrom/main.nf'
include { MARLIN                } from '../subworkflows/local/methylation_analysis/marlin.nf'
include { REMOVE_PADDING        } from '../modules/local/adaptive_specific/main.nf'
include { modifyMetaId          } from '../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'

workflow LOCAL_REALTIME {

    take:
    samplesheet             // channel: samplesheet read in from --input
    demux_samplesheet       // channel: demux samplesheet read in from --demux_samplesheet
    ref                     // channel: reference for mapping, either empty if skipping mapping, or a path
    basecall_model          // channel: model for basecalling
    ch_clin_database        // channel: clinical database for variant annotation
    bed                     // channel: bed file used for adaptive sampling regions

    main:
    if (params.skip_basecalling) {

        MAPPING(
            samplesheet,
            ref
        )
    } else {
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
    }

    if (params.realtime < 6) {                 // Before 10h of realtime sequencing, include CNV calling with QDNAseq, SV calling and Marlin

        MARLIN(
            MAPPING.out.bam,
            ref
        )

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
            bed
        )

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
            bed
        )
        // Germline variant calling using Clair3 (always uses original mapping output)
        CLAIR3_CALLING (
            MAPPING.out.bam,
            ref,
            basecall_model,
            ch_clin_database
        )
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
            bed
        )
        // Germline variant calling using Clair3 (always uses original mapping output)
        CLAIR3_CALLING (
            MAPPING.out.bam,
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
            MAPPING.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf,
            ch_panel_bin
        )
    }
}
