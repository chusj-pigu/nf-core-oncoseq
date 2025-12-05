/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLASSY_CLASSIFY } from '../../../modules/local/classy/main.nf'
include { modifyMetaId    } from '../utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN MAPPING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CLASSY {
    // Input channels:
    //   fastq_ch: Channel of tuples [meta, reads] (reads can be file or directory)
    //   ref:      Channel of tuples [meta, ref, ref_fasta, ref_fai]
    take:
    bam         // Channel: from mapping workflow
    ref       // Channel: from input samplesheet

    main:

    ch_ref = ref
        .map { meta, _ref_id, ref_fasta, _ref_fai ->
            tuple(meta, ref_fasta) }

    ch_classy_in = bam
        .join(ch_ref)

    CLASSY_CLASSIFY(ch_classy_in)

    // Add marlin to type for report
    ch_type = channel.of("Marlin")

    ch_marlin_plot = CLASSY_CLASSIFY.out.svg
        .combine(ch_type)

    ch_marlin_pred = CLASSY_CLASSIFY.out.json
        .combine(ch_type)

    // Collect versions from all modules
    ch_versions = CLASSY_CLASSIFY.out.versions

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
