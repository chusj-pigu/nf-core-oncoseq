/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLASSY_MARLIN_PILEUP           } from '../../../modules/local/classy/main.nf'
include { CLASSY_TUCAN_PILEUP            } from '../../../modules/local/classy/main.nf'
include { CLASSY_STURGEON_GENERAL_PILEUP } from '../../../modules/local/classy/main.nf'
include { CLASSY_CROSSNN_PANCAN_PILEUP   } from '../../../modules/local/classy/main.nf'
include { CLASSY_CROSSNN_CAPER_PILEUP    } from '../../../modules/local/classy/main.nf'
include { MODKIT_PILEUP                  } from '../../../modules/local/modkit/main.nf'

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

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LEUKEMIA SAMPLE CLASSIFICATION WITH MARLIN
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    MODKIT_PILEUP(input)

    ch_refid = ref
        .map { meta, ref_id, _ref_fasta, _ref_fai ->
            tuple(meta, ref_id) }

    ch_classy_in = ch_refid
        .join(MODKIT_PILEUP.out.bedmethyl)

    CLASSY_MARLIN_PILEUP(ch_classy_in)
    CLASSY_TUCAN_PILEUP(ch_classy_in)
    CLASSY_STURGEON_GENERAL_PILEUP(ch_classy_in)
    CLASSY_CROSSNN_PANCAN_PILEUP(ch_classy_in)
    CLASSY_CROSSNN_CAPER_PILEUP(ch_classy_in)

    // Add classifier types report

    ch_final_plot = addTypeToChannel(CLASSY_MARLIN_PILEUP.out.svg, "Marlin")
        .mix(addTypeToChannel(CLASSY_CROSSNN_PANCAN_PILEUP.out.svg, "CrossNN PanCancer"))
        .mix(addTypeToChannel(CLASSY_CROSSNN_CAPER_PILEUP.out.svg, "CrossNN Caper"))
        .mix(addTypeToChannel(CLASSY_STURGEON_GENERAL_PILEUP.out.svg, "Sturgeon General"))
        .mix(addTypeToChannel(CLASSY_TUCAN_PILEUP.out.svg, "Tucan"))

    ch_final_pred = addTypeToChannel(CLASSY_MARLIN_PILEUP.out.json, "Marlin")
        .mix(addTypeToChannel(CLASSY_CROSSNN_PANCAN_PILEUP.out.json, "CrossNN PanCancer"))
        .mix(addTypeToChannel(CLASSY_CROSSNN_CAPER_PILEUP.out.json, "CrossNN Caper"))
        .mix(addTypeToChannel(CLASSY_STURGEON_GENERAL_PILEUP.out.json, "Sturgeon General"))
        .mix(addTypeToChannel(CLASSY_TUCAN_PILEUP.out.json, "Tucan"))

    // Collect versions from all modules
    ch_versions = CLASSY_MARLIN_PILEUP.out.versions
        .mix(CLASSY_CROSSNN_CAPER_PILEUP.out.versions)
        .mix(CLASSY_CROSSNN_PANCAN_PILEUP.out.versions)
        .mix(CLASSY_STURGEON_GENERAL_PILEUP.out.versions)
        .mix(CLASSY_TUCAN_PILEUP.out.versions)

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
