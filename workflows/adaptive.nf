include { BASECALL_SIMPLEX  } from '../subworkflows/local/basecalling/basecall_simplex'
include { MAPPING           } from '../subworkflows/local/mapping/mapping'
include { CLAIRS_TO_CALLING } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'
include { COVERAGE_SEPARATE } from '../subworkflows/local/adaptive_specific/coverage_separate'
include { PHASING_VARIANTS  } from  '../subworkflows/local/variant_calling/phasing.nf'
include { SV_CALLING        } from  '../subworkflows/local/variant_calling/sv_calling.nf'

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
        MAPPING(samplesheet,ref)

        CLAIRS_TO_CALLING(MAPPING.out.bam,ref,chr_list,model,clin_database)

        COVERAGE_SEPARATE(MAPPING.out.bam,bed)

        // SPLIT_BAM(MAPPING.out.bam)

        //PHASING_VARIANTS(MAPPING.out.bam,ref,CLAIRS_TO_CALLING.out.vcf)

        //SV_CALLING(PHASING_VARIANTS.out.haptag_bam,ref)

    } else {

        BASECALL_SIMPLEX (
            samplesheet
        )

        MAPPING(BASECALL_SIMPLEX.out.fastq,ref)

        CLAIRS_TO_CALLING(MAPPING.out.bam,ref,chr_list,model,clin_database)

        COVERAGE_SEPARATE(MAPPING.out.bam,bed)

        PHASING_VARIANTS(MAPPING.out.bam,ref,CLAIRS_TO_CALLING.out.vcf)

        SV_CALLING(PHASING_VARIANTS.out.haptag_bam,ref)
    }
}
