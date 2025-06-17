include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'
include { MAPPING            } from '../subworkflows/local/mapping/mapping'
include { CNV_CALLING        } from '../subworkflows/local/variant_calling/cnv_calling.nf'

workflow CFDNA {

    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path

    main:

    //
    // WORKFLOW: Run pipeline
    //
    BASECALL_MULTIPLEX (
        samplesheet,
        demux
    )
    MAPPING (
        BASECALL_MULTIPLEX.out.fastq,
        ref
    )

    // CNV_CALLING (
    //         MAPPING.out.bam,
    //         ref
    //     )
}
