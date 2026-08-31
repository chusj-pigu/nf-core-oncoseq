/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLASSY_COMBINED     } from '../../../modules/local/classy/main.nf'
include { PARSE_JSON_COMBINED } from '../../../modules/local/classifier_process/main.nf'

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

    ch_classy_in = input
        .join(ref)
        .map { meta, bam, bai, ref_id, ref_fasta, _red_idx ->
            tuple(meta, bam, bai, ref_id, ref_fasta)
            }

    CLASSY_COMBINED(ch_classy_in)

    ch_classy_out_json = CLASSY_COMBINED.out.json_combined

    PARSE_JSON_COMBINED(ch_classy_out_json)

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

    ch_plot_other = ch_tumor_branched.other.join(addTypeToChannel(CLASSY_COMBINED.out.svg_blood, "Blood"))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(CLASSY_COMBINED.out.svg_brain, "Brain")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(CLASSY_COMBINED.out.svg_solid, "Solid")))
        .groupTuple()
        .transpose()

    ch_json_blood = ch_tumor_branched.blood
        .join(addTypeToChannel(PARSE_JSON_COMBINED.out.alma, "Alma"))
        .mix(ch_tumor_branched.blood.join(addTypeToChannel(PARSE_JSON_COMBINED.out.lamprey, "Lamprey")))
        .mix(ch_tumor_branched.blood.join(addTypeToChannel(PARSE_JSON_COMBINED.out.marlin, "Marlin")))
        .groupTuple()
        .transpose()

    ch_json_brain = ch_tumor_branched.brain
        .join(addTypeToChannel(PARSE_JSON_COMBINED.out.capper, "CrossNN Caper"))
        .mix(ch_tumor_branched.brain.join(addTypeToChannel(PARSE_JSON_COMBINED.out.mpact, "MPACT")))
        .mix(ch_tumor_branched.brain.join(addTypeToChannel(PARSE_JSON_COMBINED.out.sturgeon_brainstem, "Sturgeon Brainstem")))
        .mix(ch_tumor_branched.brain.join(addTypeToChannel(PARSE_JSON_COMBINED.out.sturgeon_general, "Sturgeon General")))
        .groupTuple()
        .transpose()

    ch_json_solid = ch_tumor_branched.solid
        .join(addTypeToChannel(PARSE_JSON_COMBINED.out.pancan, "CrossNN PanCancer"))
        .mix(ch_tumor_branched.solid.join(addTypeToChannel(PARSE_JSON_COMBINED.out.tucan, "Tucan")))
        .groupTuple()
        .transpose()

    ch_json_other = ch_tumor_branched.other
        .join(addTypeToChannel(PARSE_JSON_COMBINED.out.alma, "Alma"))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(PARSE_JSON_COMBINED.out.capper, "CrossNN Caper")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(PARSE_JSON_COMBINED.out.pancan, "CrossNN PanCancer")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(PARSE_JSON_COMBINED.out.lamprey, "Lamprey")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(PARSE_JSON_COMBINED.out.marlin, "Marlin")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(PARSE_JSON_COMBINED.out.mpact, "MPACT")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(PARSE_JSON_COMBINED.out.sturgeon_brainstem, "Sturgeon Brainstem")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(PARSE_JSON_COMBINED.out.sturgeon_general, "Sturgeon General")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(PARSE_JSON_COMBINED.out.tucan, "Tucan")))
        .groupTuple()
        .transpose()

    if (params.cfdna) {
        ch_final_plot = ch_plot_blood
            .mix(ch_plot_brain)
            .mix(ch_plot_solid)
            .mix(ch_plot_other)
            .mix(addTypeToChannel(CLASSY_COMBINED.out.svg_nanomix, "Tissue of origin"))

        ch_final_pred = ch_json_blood
            .mix(ch_json_brain)
            .mix(ch_json_solid)
            .mix(ch_json_other)
            .mix(addTypeToChannel(PARSE_JSON_COMBINED.out.nanomix, "Nanomix"))

    } else {
        ch_final_plot = ch_plot_blood
            .mix(ch_plot_brain)
            .mix(ch_plot_solid)
            .mix(ch_plot_other)

        ch_final_pred = ch_json_blood
            .mix(ch_json_brain)
            .mix(ch_json_solid)
            .mix(ch_json_other)
    }

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
