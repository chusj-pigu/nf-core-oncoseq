include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'
include { BASECALL_SIMPLEX   } from '../subworkflows/local/basecalling/basecall_simplex'
include { READS_FILTER       } from '../subworkflows/local/read_processing/chopper'
include { MAPPING            } from '../subworkflows/local/mapping/mapping'
include { CNV_CALLING        } from '../subworkflows/local/variant_calling/cnv_calling.nf'
include { ICHORCNA_CALLING   } from '../subworkflows/local/variant_calling/ichor_calling.nf'

workflow CFDNA {

    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    max_len
    minqs
    ichor_bin
    mapq_wig

    main:

    //
    // WORKFLOW: Run pipeline
    //

    if (params.demux) {
        BASECALL_MULTIPLEX (
            samplesheet,
            demux
        )
        READS_FILTER (
            samplesheet,
            BASECALL_MULTIPLEX.out.fastq,
            max_len,
            minqs
        )
    } else if (params.skip_basecalling) {

        READS_FILTER (
            samplesheet,
            samplesheet,
            max_len,
            minqs
        )
    } else {

         BASECALL_SIMPLEX (
            samplesheet
        )
        READS_FILTER (
            samplesheet,
            BASECALL_SIMPLEX.out.fastq,
            max_len,
            minqs
        )
    }

    MAPPING (
        READS_FILTER.out.reads,
        ref
    )

    CNV_CALLING (
        MAPPING.out.bam,
        ref
    )

    ICHORCNA_CALLING (
        MAPPING.out.bam,
        ref,
        samplesheet,
        ichor_bin,
        mapq_wig
    )

}
