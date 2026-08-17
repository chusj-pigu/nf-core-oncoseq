/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLASSY_COMBINED           } from '../../../modules/local/classy/main.nf'
include { PARSE_JSON_COMBINED       } from '../../../modules/local/classifier_process/main.nf'

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
    tumor_types // Channel: from input samplesheet

    main:

    ch_versions = channel.empty()

    ch_ref = ref
        .map { meta, ref_id, ref_fasta, _ref_fai ->
            tuple(meta, ref_id, ref_fasta) }

    ch_classy_in = input
        .join(ch_ref)

    CLASSY_COMBINED(ch_classy_in)

    PARSE_JSON_COMBINED(CLASSY_COMBINED.out.json_combined)

    tumor_types
        .branch { meta, type ->
            blood: type == "blood"
                return meta
            solid: type == "solid"
                return meta
            brain: type == "brain"
                return meta
            other: true
                return meta
        }
        .set { ch_tumor_branched }

    ch_plot_blood = ch_tumor_branched.blood
        .join(addTypeToChannel(CLASSY_COMBINED.out.svg_blood, "Blood"))

    ch_plot_brain = ch_tumor_branched.brain
        .join(addTypeToChannel(CLASSY_COMBINED.out.svg_brain, "Brain"))

    ch_plot_solid = ch_tumor_branched.solid
        .join(addTypeToChannel(CLASSY_COMBINED.out.svg_solid, "Solid"))

    ch_plot_other = ch_tumor_branched.other
        .join(addTypeToChannel(CLASSY_COMBINED.out.svg_blood, "Blood"))
        .join(addTypeToChannel(CLASSY_COMBINED.out.svg_brain, "Brain"))
        .join(addTypeToChannel(CLASSY_COMBINED.out.svg_solid, "Solid"))
        .transpose()

    ch_final_plot = ch_plot_blood
        .mix(ch_plot_brain)
        .mix(ch_plot_solid)
        .mix(ch_plot_other)

    ch_final_pred = addTypeToChannel(PARSE_JSON_COMBINED.out.alma, "Alma")
        .mix(addTypeToChannel(PARSE_JSON_COMBINED.out.capper, "CrossNN Capper"))
        .mix(addTypeToChannel(PARSE_JSON_COMBINED.out.pancan, "CrossNN PanCancer"))
        .mix(addTypeToChannel(PARSE_JSON_COMBINED.out.lamprey, "Lamprey"))
        .mix(addTypeToChannel(PARSE_JSON_COMBINED.out.marlin, "Marlin"))
        .mix(addTypeToChannel(PARSE_JSON_COMBINED.out.mpact, "MPACT"))
        .mix(addTypeToChannel(PARSE_JSON_COMBINED.out.nanomix, "Nanomix"))
        .mix(addTypeToChannel(PARSE_JSON_COMBINED.out.sturgeon_brainstem, "Sturgeon Brainstem"))
        .mix(addTypeToChannel(PARSE_JSON_COMBINED.out.sturgeon_general, "Sturgeon General"))
        .mix(addTypeToChannel(PARSE_JSON_COMBINED.out.tucan, "Tucan"))

    // Collect versions from all modules
    ch_versions = CLASSY_COMBINED.out.versions
        .mix(PARSE_JSON_COMBINED.out.versions)

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
