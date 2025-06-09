include { BASECALL_SIMPLEX  } from '../subworkflows/local/basecalling/basecall_simplex'
include { MAPPING           } from '../subworkflows/local/mapping/mapping'
include { CLAIRS_TO_CALLING } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'
include { COVERAGE_SEPARATE } from '../subworkflows/local/adaptive_specific/coverage_separate'
include { PHASING_VARIANTS  } from  '../subworkflows/local/variant_calling/phasing.nf'
include { SV_CALLING        } from  '../subworkflows/local/variant_calling/sv_calling.nf'
include { CNV_CALLING       } from  '../subworkflows/local/variant_calling/cnv_calling.nf'
include { SPLIT_BAMS_TIME   } from '../subworkflows/local/time_series_evaluation/split_bams.nf'
include { SPLIT_BAMS_TIME_FASTQ   } from '../subworkflows/local/time_series_evaluation/split_bams_fastq.nf'


workflow ADAPTIVE {

    take:
    samplesheet // channel: samplesheet read in from --input
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    model
    clin_database
    bed         // channel: from path read from params.bed, bed file used for adaptive sampling

    main:




    //
    // WORKFLOW: Run pipeline
    //

    if (params.skip_basecalling) {
        SPLIT_BAMS_TIME_FASTQ(samplesheet)

        MAPPING(
            SPLIT_BAMS_TIME_FASTQ.out.ch_fastq_out,
            ref
            )

        CLAIRS_TO_CALLING(
            MAPPING.out.bam,
            ref,
            model,
            clin_database
            )

        COVERAGE_SEPARATE(
            MAPPING.out.bam,
            bed
            )

        PHASING_VARIANTS(MAPPING.out.bam,ref,CLAIRS_TO_CALLING.out.vcf)

        SV_CALLING(PHASING_VARIANTS.out.haptag_bam,ref)

        CNV_CALLING(MAPPING.out.bam,ref)

    } else {

        BASECALL_SIMPLEX (
            samplesheet
        )

        MAPPING(BASECALL_SIMPLEX.out.fastq,ref)

        if (params.time_series) {
            // If time series evaluation is enabled, split BAMs into time intervals
            SPLIT_BAMS_TIME(
                MAPPING.out.bam,
                ref,
                bed
            )
            
            // Use time series outputs for variant calling
            ch_bam_for_calling = SPLIT_BAMS_TIME.out.bam
            ch_ref_for_calling = SPLIT_BAMS_TIME.out.ref
            ch_bed = SPLIT_BAMS_TIME.out.bed
        } else {
            // If not, use the full BAM directly
            ch_bam_for_calling = MAPPING.out.bam
            ch_ref_for_calling = ref
            ch_bed = bed
        }

        CLAIRS_TO_CALLING(
            ch_bam_for_calling,
            ch_ref_for_calling,
            model,
            clin_database
            )

        COVERAGE_SEPARATE(ch_bam_for_calling,ch_bed)

        PHASING_VARIANTS(
            ch_bam_for_calling,
            ch_ref_for_calling,
            CLAIRS_TO_CALLING.out.vcf
            )

        SV_CALLING(
            PHASING_VARIANTS.out.haptag_bam,
            ch_ref_for_calling
        )

        CNV_CALLING(
            ch_bam_for_calling,
            ch_ref_for_calling
        )
    }
}
