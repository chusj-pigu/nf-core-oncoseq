/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLASSY_MARLIN           } from '../../../modules/local/classy/main.nf'
include { CLASSY_TUCAN            } from '../../../modules/local/classy/main.nf'
include { CLASSY_STURGEON_GENERAL } from '../../../modules/local/classy/main.nf'
include { CLASSY_CROSSNN_PANCAN   } from '../../../modules/local/classy/main.nf'
include { CLASSY_CROSSNN_CAPER    } from '../../../modules/local/classy/main.nf'

def addTypeToChannel(ch, type) {
    ch = ch.map { meta, data -> tuple(meta, data, type) }
    return ch
}

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
    input         // Channel: from mapping workflow
    ref       // Channel: from input samplesheet

    main:

    ch_versions = channel.empty()

    ch_ref = ref
        .map { meta, ref_id, ref_fasta, _ref_fai ->
            tuple(meta, ref_id, ref_fasta) }

    ch_classy_in = input
        .join(ch_ref)

    CLASSY_MARLIN(ch_classy_in)
    CLASSY_TUCAN(ch_classy_in)
    CLASSY_STURGEON_GENERAL(ch_classy_in)
    CLASSY_CROSSNN_PANCAN(ch_classy_in)
    CLASSY_CROSSNN_CAPER(ch_classy_in)

    // Add classifier types report

    ch_final_plot = addTypeToChannel(CLASSY_MARLIN.out.svg, "Marlin")
        .mix(addTypeToChannel(CLASSY_CROSSNN_PANCAN.out.svg, "CrossNN PanCancer"))
        .mix(addTypeToChannel(CLASSY_CROSSNN_CAPER.out.svg, "CrossNN Caper"))
        .mix(addTypeToChannel(CLASSY_STURGEON_GENERAL.out.svg, "Sturgeon General"))
        .mix(addTypeToChannel(CLASSY_TUCAN.out.svg, "Tucan"))

    ch_final_pred = addTypeToChannel(CLASSY_MARLIN.out.json, "Marlin")
        .mix(addTypeToChannel(CLASSY_CROSSNN_PANCAN.out.json, "CrossNN PanCancer"))
        .mix(addTypeToChannel(CLASSY_CROSSNN_CAPER.out.json, "CrossNN Caper"))
        .mix(addTypeToChannel(CLASSY_STURGEON_GENERAL.out.json, "Sturgeon General"))
        .mix(addTypeToChannel(CLASSY_TUCAN.out.json, "Tucan"))

    // Collect versions from all modules
    ch_versions = CLASSY_MARLIN.out.versions
        .mix(CLASSY_CROSSNN_CAPER.out.versions)
        .mix(CLASSY_CROSSNN_PANCAN.out.versions)
        .mix(CLASSY_STURGEON_GENERAL.out.versions)
        .mix(CLASSY_TUCAN.out.versions)

    emit:
    plot     = ch_final_plot
    pred     = ch_final_pred
    versions = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
