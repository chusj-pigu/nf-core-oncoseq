include { BASECALL_SIMPLEX  } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'
include { MAPPING           } from '../subworkflows/local/mapping/mapping'
include { CLAIRS_TO_CALLING } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'
include { CLAIR3_CALLING } from '../subworkflows/local/variant_calling/clair3_calling.nf'
include { COVERAGE_SEPARATE } from '../subworkflows/local/adaptive_specific/coverage_separate'
include { PHASING_VARIANTS as PHASING_SOMATIC  } from  '../subworkflows/local/variant_calling/phasing.nf'
include { PHASING_VARIANTS as PHASING_GERMLINE  } from  '../subworkflows/local/variant_calling/phasing.nf'
include { SV_CALLING        } from  '../subworkflows/local/variant_calling/sv_calling.nf'
include { CNV_CALLING       } from  '../subworkflows/local/variant_calling/cnv_calling.nf'

workflow ADAPTIVE {

    take:
    samplesheet             // channel: samplesheet read in from --input
    demux_samplesheet       // channel : demux samplesheet read in from --demux_samplesheet
    ref                     // channel : reference for mapping, either empty if skipping mapping, or a path
    clairs_model
    basecall_model
    ch_clin_database
    bed                     // channel: from path read from params.bed, bed file used for adaptive sampling

    main:

    //
    // WORKFLOW: Run pipeline
    //

    if (params.skip_basecalling) {

        MAPPING (
            samplesheet,
            ref
        )

        COVERAGE_SEPARATE (
            MAPPING.out.bam,
            bed
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
            PHASING_GERMLINE.out.haptag_bam,
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

        COVERAGE_SEPARATE (
            MAPPING.out.bam,
            bed
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
            PHASING_GERMLINE.out.haptag_bam,
            ref
        )

        CNV_CALLING (
            MAPPING.out.bam,
            ref
        )

    }
}
