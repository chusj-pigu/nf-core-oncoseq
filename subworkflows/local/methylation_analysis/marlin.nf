/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MARLIN_PILEUP  } from '../../../modules/local/marlin/main.nf'             // modkit mpileup
include { MARLIN_MERGE   } from '../../../modules/local/marlin/main.nf'
include { MARLIN_PREDICT } from '../../../modules/local/marlin/main.nf'
include { MARLIN_PLOT    } from '../../../modules/local/marlin/main.nf'
include { modifyMetaId   } from '../utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN MAPPING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow MARLIN {
    // Input channels:
    //   fastq_ch: Channel of tuples [meta, reads] (reads can be file or directory)
    //   ref:      Channel of tuples [meta, ref, ref_fasta, ref_fai]
    take:
    bam         // Channel: from mapping workflow
    ref       // Channel: from input samplesheet

    main:

    ch_ref_pileup = ref
        .map { meta, ref_id, ref_fasta, _ref_fai ->
            tuple(meta, ref_id, ref_fasta) }

    ch_pileup_in = bam
        .join(ch_ref_pileup)

    MARLIN_PILEUP(ch_pileup_in)

    ch_ref_merge = ref
        .map { meta, ref_id, _ref_fasta, _ref_fai ->
            tuple(meta, ref_id )}

    ch_merge_in = MARLIN_PILEUP.out.bedmethyl
        .join(ch_ref_merge)

    MARLIN_MERGE(ch_merge_in)

    MARLIN_PREDICT(MARLIN_MERGE.out.merged_bedmethyl)

    MARLIN_PLOT(MARLIN_PREDICT.out.pred)

    // Add marlin to type for report
    ch_type = Channel.of("Marlin")

    ch_marlin_plot = MARLIN_PLOT.out.pred_pdf
        .combine(ch_type)

    ch_marlin_pred = MARLIN_PREDICT.out.pred
        .combine(ch_type)

    // Collect versions from all modules
    ch_versions = MARLIN_PILEUP.out.versions
        .mix(MARLIN_PREDICT.out.versions)

    emit:
    plot     = ch_marlin_plot                 // TODO: Quarto report
    pred     = ch_marlin_pred
    versions = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
