include { BASECALL_MULTIPLEX   } from '../subworkflows/local/basecalling/basecall_multiplex'
include { BASECALL_SIMPLEX     } from '../subworkflows/local/basecalling/basecall_simplex'
include { READS_FILTER         } from '../subworkflows/local/read_processing/chopper'
include { MAPPING              } from '../subworkflows/local/mapping/mapping'
include { TIDEHUNTER_CONCENSUS } from '../subworkflows/local/read_processing/tidehunter'
include { CNV_CALLING          } from '../subworkflows/local/variant_calling/cnv_calling.nf'
include { ICHORCNA_CALLING     } from '../subworkflows/local/variant_calling/ichor_calling.nf'
include { MARLIN               } from '../subworkflows/local/methylation_analysis/marlin.nf'

workflow CFDNA {

    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    cfdna_samplesheet   // channel : from demux or samplesheeet
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
            demux,
            ref
        )

        ch_fastq = BASECALL_SIMPLEX.out.fastq

    } else if (params.skip_basecalling) {

        ch_fastq = samplesheet

    } else {

         BASECALL_SIMPLEX (
            samplesheet
        )

        ch_fastq = BASECALL_SIMPLEX.out.fastq

    }

    if (params.rca) {
        TIDEHUNTER_CONCENSUS(ch_fastq)

        ch_fastq_processed = TIDEHUNTER_CONCENSUS.out.fastq

    } else {
        READS_FILTER (
            cfdna_samplesheet,
            ch_fastq,
            max_len,
            minqs
        )

        ch_fastq_processed = READS_FILTER.out.reads
    }

    MAPPING (
        ch_fastq_processed,
        ref
    )

    MARLIN(
        MAPPING.out.bam,
        ref
    )

    CNV_CALLING (
        MAPPING.out.bam,
        ref
    )

    ICHORCNA_CALLING (
        MAPPING.out.bam,
        cfdna_samplesheet,
        ichor_bin,
        mapq_wig
    )

}
