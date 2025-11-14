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

// Reporting
include { MIDNIGHT_REPORT } from '../subworkflows/local/report/final_report.nf'
include { CFNDA_REPORT    } from '../subworkflows/local/report/cfdna.nf'
include { CLASSIFIER_REPORT    } from '../subworkflows/local/report/methylation.nf'

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
    ch_id

    main:

    //
    // WORKFLOW: Run pipeline
    //

    ch_versions = channel.empty()

    if (params.demux) {

        BASECALL_MULTIPLEX (
            samplesheet,
            demux
        )

        ch_fastq = BASECALL_MULTIPLEX.out.fastq

        ch_versions = BASECALL_MULTIPLEX.out.versions

    } else if (params.skip_basecalling) {

        ch_fastq = samplesheet

    } else {

         BASECALL_SIMPLEX (
            samplesheet
        )

        ch_fastq = BASECALL_SIMPLEX.out.fastq

        ch_versions = BASECALL_SIMPLEX.out.versions

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

        ch_mapping_t2t = ch_tumor_type.cns
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
            MAPPING_HG.out.bam,
            ref
        )

        ICHORCNA_CALLING (
            MAPPING_HG.out.bam,
            cfdna_samplesheet,
            ichor_bin,
            mapq_wig
        )

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            COLLECT VERSIONS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */


        ch_versions = ch_versions
            .mix(MAPPING_HG.out.versions)
            .mix(MAPPING_T2T.out.versions)
            .mix(MARLIN.out.versions)
            .mix(STURGEON.out.versions)
            .mix(CNV_CALLING.out.versions)
            .mix(ICHORCNA_CALLING.out.versions)

        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            COMPILE SECTIONS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_classifiers_plots = STURGEON.out.plot
            .mix(MARLIN.out.plot)

        ch_classifiers_pred = STURGEON.out.pred
            .mix(MARLIN.out.pred)

        CFNDA_REPORT(
            ch_id,
            READS_FILTER.out.stats,
            MAPPING_HG.out.coverage,
            ICHORCNA_CALLING.out.ichorcna_plot
        )

        CLASSIFIER_REPORT(
            ch_classifiers_plots,
            ch_classifiers_pred
        )

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            COLLECT SECTIONS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_classifier_out = CLASSIFIER_REPORT.out.sections
            .join(ch_id)
            .map { meta, section, inputs, quatro, project ->
                tuple(id:project, section, inputs, quarto) }

        // Collect sections from all analysis steps
        ch_sections = CFNDA_REPORT.out.sections
            .mix(ch_classifier_out)


        ch_mode = channel.of("cfDNA")

        MIDNIGHT_REPORT(
            ch_id,
            ch_sections,
            ch_versions,
            ch_mode
        )
    }
}
