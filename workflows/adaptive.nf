include { BASECALL_SIMPLEX  } from '../subworkflows/local/basecalling/basecall_simplex'
include { MAPPING           } from '../subworkflows/local/mapping/mapping'
include { CLAIRS_TO_CALLING } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'
include { COVERAGE_SEPARATE } from '../subworkflows/local/adaptive_specific/coverage_separate'
include { PHASING_VARIANTS  } from '../subworkflows/local/variant_calling/phasing.nf'
include { SV_CALLING        } from '../subworkflows/local/variant_calling/sv_calling.nf'
include { SPLIT_BAMS_TIME   } from '../subworkflows/local/time_series_evaluation/split_bams.nf'
include { SPLIT_BAMS_TIME_FASTQ   } from '../subworkflows/local/time_series_evaluation/split_bams_fastq.nf'


workflow ADAPTIVE {

    take:
    samplesheet // channel: samplesheet read in from --input
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    chr_list
    model
    clin_database
    bed         // channel: from path read from params.bed, bed file used for adaptive sampling

    main:

    //
    // WORKFLOW: Run pipeline
    //

    if (params.skip_basecalling) {
        SPLIT_BAMS_TIME_FASTQ(samplesheet)

        MAPPING(SPLIT_BAMS_TIME_FASTQ.out.ch_fastq_out,ref)

        CLAIRS_TO_CALLING(MAPPING.out.bam,ref,chr_list,model,clin_database)

        COVERAGE_SEPARATE(MAPPING.out.bam,bed)

        //PHASING_VARIANTS(MAPPING.out.bam,ref,CLAIRS_TO_CALLING.out.vcf)

        //SV_CALLING(PHASING_VARIANTS.out.haptag_bam,ref)

    } else {

        BASECALL_SIMPLEX (
            samplesheet
        )

        MAPPING(BASECALL_SIMPLEX.out.fastq,ref)

        SPLIT_BAMS_TIME(MAPPING.out.bam, ref, clin_database)

        // SPLIT_BAMS_TIME.out.bam.view()
        // SPLIT_BAMS_TIME.out.ref.view()
        // SPLIT_BAMS_TIME.out.clin_database.view()

        ref = (SPLIT_BAMS_TIME.out.ref)
        clin_database = (SPLIT_BAMS_TIME.out.clin_database)


        CLAIRS_TO_CALLING(SPLIT_BAMS_TIME.out.bam,ref,chr_list,model,clin_database)

        COVERAGE_SEPARATE(SPLIT_BAMS_TIME.out.bam,bed)

        PHASING_VARIANTS(SPLIT_BAMS_TIME.out.bam,ref,CLAIRS_TO_CALLING.out.vcf)

        SV_CALLING(PHASING_VARIANTS.out.haptag_bam,ref)
    }
}
