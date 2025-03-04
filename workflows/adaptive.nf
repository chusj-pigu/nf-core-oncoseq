include { BASECALL_SIMPLEX  } from '../subworkflows/local/basecalling/basecall_simplex'
include { MAPPING           } from '../subworkflows/local/mapping/mapping'
include { COVERAGE_SEPARATE } from '../subworkflows/local/adaptive_specific/coverage_separate'

workflow ADAPTIVE {

    take:
    samplesheet // channel: samplesheet read in from --input
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    bed         // channel: from path read from params.bed, bed file used for adaptive sampling

    main:

    //
    // WORKFLOW: Run pipeline
    //
    BASECALL_SIMPLEX (
        samplesheet
    )

    MAPPING(BASECALL_SIMPLEX.out.fastq,ref)
    COVERAGE_SEPARATE(MAPPING.out.bam,bed)
}
