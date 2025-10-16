include { BASECALL_MULTIPLEX     } from '../subworkflows/local/basecalling/basecall_multiplex'
include { BASECALL_SIMPLEX       } from '../subworkflows/local/basecalling/basecall_simplex'
include { READS_FILTER           } from '../subworkflows/local/read_processing/reads_filter.nf'
include { MAPPING as MAPPING_HG  } from '../subworkflows/local/mapping/mapping'
include { MAPPING as MAPPING_T2T } from '../subworkflows/local/mapping/mapping'
include { TIDEHUNTER_CONCENSUS   } from '../subworkflows/local/read_processing/tidehunter'
include { CNV_CALLING            } from '../subworkflows/local/variant_calling/cnv_calling.nf'
include { ICHORCNA_CALLING       } from '../subworkflows/local/variant_calling/ichor_calling.nf'
include { MARLIN                 } from '../subworkflows/local/methylation_analysis/marlin.nf'
include { STURGEON               } from '../subworkflows/local/methylation_analysis/sturgeon.nf'

workflow CFDNA {

    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    cfdna_samplesheet   // channel : from demux or samplesheeet
    tumor_type      // channel: to decide if marlin or sturgeon is run
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    ref_t2t
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

        ch_fastq = BASECALL_MULTIPLEX.out.fastq

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

        ch_tumor_type = tumor_type
            .branch { meta, tumor ->
                leukemia: tumor == "leukemia"
                cns: tumor == "cns"
                other: tumor == "other"
            }

        ch_mapping_t2t = tumor_type.cns
            .join(ch_fastq_processed)
            .map { meta, _tumor, input ->
                tuple(meta, input)}

        MAPPING_HG(
            ch_fastq_processed,
            ref
        )

        MAPPING_T2T(
            ch_mapping_t2t,
            ref_t2t
        )

        // Run Marlin on leukemia samples only:
        ch_marlin_bam = MAPPING_HG.out.bam
            .join(ch_tumor_type.leukemia)
            .map { meta, bam, bai, _tumor ->
                tuple(meta, bam, bai)}

        MARLIN(
            ch_marlin_bam,
            ref
        )

        STURGEON(
            MAPPING_T2T.out.bam
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
}
