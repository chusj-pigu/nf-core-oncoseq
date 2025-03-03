include { BASECALL_SIMPLEX } from '../subworkflows/local/basecalling/basecall_simplex'
include { MAPPING          } from '../subworkflows/local/mapping/mapping'
include { CLAIRS_TO_CALLING } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'

workflow ADAPTIVE {

    take:
    samplesheet // channel: samplesheet read in from --input
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    chr_list
    model

    main:

    //
    // WORKFLOW: Run pipeline
    //
    BASECALL_SIMPLEX (
        samplesheet
    )

    MAPPING(BASECALL_SIMPLEX.out.fastq,ref)

    CLAIRS_TO_CALLING(MAPPING.out.bam,ref,chr_list,model)
}
