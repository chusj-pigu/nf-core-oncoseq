include { BASECALL_SIMPLEX } from '../subworkflows/local/basecalling/basecall_simplex'
include { MAPPING          } from '../subworkflows/local/mapping/mapping'

workflow ADAPTIVE {

    take:
    samplesheet // channel: samplesheet read in from --input
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path

    main:

    //
    // WORKFLOW: Run pipeline
    //
    BASECALL_SIMPLEX (
        samplesheet
    )

    MAPPING(BASECALL_SIMPLEX.out.fastq,ref)
}
