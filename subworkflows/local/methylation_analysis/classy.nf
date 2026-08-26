/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLASSY_COMBINED as CLASSY_ALL         } from '../../../modules/local/classy/main.nf'
include { CLASSY_COMBINED as CLASSY_NO_STURGEON } from '../../../modules/local/classy/main.nf'
include { PARSE_JSON_COMBINED                   } from '../../../modules/local/classifier_process/main.nf'

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

    // Only run sturgeon for samples that have panel coverage of under 10X, otherwise the process is too lengthy.

    ch_classy_in = input
        .join(ref)
        .branch { meta, bam, bai, cov, ref_id, ref_fasta, _ref_fai ->
            low_cov: cov <= 10
                return tuple(meta,bam,bai,ref_id,ref_fasta)
            high_cov: true
                return tuple(meta,bam,bai,ref_id,ref_fasta)
            }

    CLASSY_ALL(ch_classy_in.low_cov)
    CLASSY_NO_STURGEON(ch_classy_in.high_cov)

    ch_classy_out_json = CLASSY_ALL.out.json_combined
        .mix(CLASSY_NO_STURGEON.out.json_combined)

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
        .join(addTypeToChannel(CLASSY_ALL.out.svg_blood, "Blood"))
        .mix(ch_tumor_branched.blood.join(addTypeToChannel(CLASSY_NO_STURGEON.out.svg_blood, "Blood")))

    ch_plot_brain = ch_tumor_branched.brain
        .join(addTypeToChannel(CLASSY_ALL.out.svg_brain, "Brain"))
        .mix(ch_tumor_branched.brain.join(addTypeToChannel(CLASSY_NO_STURGEON.out.svg_brain, "Brain")))

    ch_plot_solid = ch_tumor_branched.solid
        .join(addTypeToChannel(CLASSY_ALL.out.svg_solid, "Solid"))
        .mix(ch_tumor_branched.solid.join(addTypeToChannel(CLASSY_NO_STURGEON.out.svg_solid, "Solid")))

    ch_plot_other = ch_tumor_branched.other.join(addTypeToChannel(CLASSY_ALL.out.svg_blood, "Blood"))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(CLASSY_NO_STURGEON.out.svg_blood, "Blood")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(CLASSY_ALL.out.svg_brain, "Brain")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(CLASSY_NO_STURGEON.out.svg_brain, "Brain")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(CLASSY_ALL.out.svg_solid, "Solid")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(CLASSY_NO_STURGEON.out.svg_solid, "Solid")))
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
        .mix(ch_tumor_branched.other.join(addTypeToChannel(PARSE_JSON_COMBINED.out.tucan, "CrossNN Caper")))
        .mix(ch_tumor_branched.other.join(addTypeToChannel(PARSE_JSON_COMBINED.out.tucan, "CrossNN PanCancer")))
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
            .mix(addTypeToChannel(CLASSY_ALL.out.svg_nanomix, "Tissue of origin"))
            .mix(addTypeToChannel(CLASSY_NO_STURGEON.out.svg_nanomix, "Tissue of origin"))

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
    ch_versions = CLASSY_NO_STURGEON.out.versions
        .mix(CLASSY_ALL.out.versions)
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
