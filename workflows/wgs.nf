// Basecalling subworkflows
include { BASECALL_SIMPLEX   } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'

// Core analysis subworkflows
include { MAPPING            } from '../subworkflows/local/mapping/mapping'

// Variant calling subworkflows
include { CLAIRS_TO_CALLING                    } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'
include { CLAIR3_CALLING                       } from '../subworkflows/local/variant_calling/clair3_calling.nf'
include { PHASING_VARIANTS as PHASING_SOMATIC  } from  '../subworkflows/local/variant_calling/phasing.nf'
include { PHASING_VARIANTS as PHASING_GERMLINE } from  '../subworkflows/local/variant_calling/phasing.nf'
include { SV_CALLING                           } from '../subworkflows/local/variant_calling/sv_calling.nf'
include { CNV_CALLING                          } from '../subworkflows/local/variant_calling/cnv_calling.nf'
include { modifyMetaId                         } from '../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'

workflow WGS {

    take:
    samplesheet             // channel: samplesheet read in from --input
    demux_samplesheet       // channel : demux samplesheet read in from --demux_samplesheet
    ref                     // channel : reference for mapping, either empty if skipping mapping, or a path
    clairs_model
    basecall_model
    ch_clin_database

    main:

    //
    // WORKFLOW: Run pipeline
    //

    if (params.skip_basecalling) {

        MAPPING (
            samplesheet,
            ref
        )

        CLAIRS_TO_CALLING (
            MAPPING.out.bam,
            ref,
            clairs_model,
            ch_clin_database
        )

        CLAIR3_CALLING (
            MAPPING.out.bam,
            ref,
            basecall_model,
            ch_clin_database
        )

        PHASING_SOMATIC (
            MAPPING.out.bam,
            ref,
            CLAIRS_TO_CALLING.out.vcf
        )

        PHASING_GERMLINE (
            MAPPING.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf
        )

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

        CNV_CALLING (
            MAPPING.out.bam,
            ref
        )

    } else {

        if (params.demux != null) {

            BASECALL_MULTIPLEX (
                samplesheet,
                demux_samplesheet
            )

            MAPPING (
                BASECALL_MULTIPLEX.out.fastq,
                ref
            )
        } else {

            BASECALL_SIMPLEX (
                samplesheet
            )

            MAPPING (
                BASECALL_SIMPLEX.out.fastq,
                ref
            )
        }

        CLAIRS_TO_CALLING (
            MAPPING.out.bam,
            ref,
            clairs_model,
            ch_clin_database
        )

        CLAIR3_CALLING (
            MAPPING.out.bam,
            ref,
            basecall_model,
            ch_clin_database
        )

        PHASING_SOMATIC (
            MAPPING.out.bam,
            ref,
            CLAIRS_TO_CALLING.out.vcf
        )

        PHASING_GERMLINE (
            MAPPING.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf
        )

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

        CNV_CALLING (
            MAPPING.out.bam,
            ref
        )
    }
}
