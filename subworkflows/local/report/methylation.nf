/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_SECTION        } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE_TABS    } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TEXT           } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE_TABS     } from '../../../modules/local/quarto/main.nf'

workflow CLASSIFIER_REPORT {

    take:
    ch_methylation_plot
    ch_methylation_call
    tumor_type

    main:

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QUARTO FIGURES with tabs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_classifiers_per_sample = ch_methylation_call
        .map { meta, _tables, classifier -> tuple(meta, classifier) }
        .groupTuple()

    ch_methylation_text_classified = tumor_type
        .join(ch_classifiers_per_sample)
        .map { meta, type, classifiers ->
            def names = classifiers.unique().sort()
            def classifier_list = names.size() > 1
                ? names[0..-2].join(", ") + " and " + names[-1]
                : names[0]

            def section = "Methylation"
            def process = "Methylation-text-${meta.id}"
            tuple(meta,
                "Tumor classifier predictions based on DNA methylation profiles. Top 5 scoring predictions are shown, and the classifiers used for ${type} tumor types are ${classifier_list}.",
                section, process)
        }

    QUARTO_TEXT(
        ch_methylation_text_classified
    )

    ch_quarto_fig = ch_methylation_plot
        .groupTuple()
        .map { meta, plots, types ->

            def tabs = types
            def captions = types.collect { type ->
                def caption_type = type == "Tissue of origin" ?
                "${type} prediction by Nanomix model" :
                "Tumor classifier prediction for ${type} tumors"
                caption_type
            }
            def section = "Methylation"
            def process = "meth-plot-${meta.id}"

            // Return a new tuple with the additional variables
            return [meta, section, process, plots, tabs, captions]
        }

    QUARTO_FIGURE_TABS(
        ch_quarto_fig
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QUARTO TABLES from json outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    def cutoffMap = [
        'Sturgeon General'   : 0.8,
        'Sturgeon Brainstem' : 0.8,
        'Tucan'              : 0.7,
        'Marlin'             : 0.8,
        'CrossNN Caper'      : 0.2,
        'CrossNN PanCancer'  : 0.15,
        'Alma'               : 0.2,
        'Lamprey'            : 0.9,
        'MPACT'              : 0.7,
        'Nanomix'            : 0.8
    ]

    ch_tables_in_quarto = ch_methylation_call
        .groupTuple()
        .map { meta, tables, types ->
            // sort tables and labels together so tabs have a stable order
            def pairs  = [tables, types].transpose().sort { it[1] }
            def tables_sorted = pairs.collect { it[0] }
            def tab           = pairs.collect { it[1] }
            def captions      = tab.collect { captionFor(it, cutoffMap) }

            tuple(meta, "Methylation", "methylation-calls-${meta.id}",
                tables_sorted, tab, captions, "")
        }

    QUARTO_TABLE_TABS(ch_tables_in_quarto)


    ch_section_inputs = QUARTO_TEXT.out.quarto_text
        .mix(QUARTO_FIGURE_TABS.out.quarto_figure)
        .mix(QUARTO_TABLE_TABS.out.quarto_table)

    ch_section_meh = ch_section_inputs
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths, 'Methylation based classification']
        }

    QUARTO_SECTION(
        ch_section_meh
    )

    ch_sections = QUARTO_SECTION.out.quarto_section

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    emit:
    sections = ch_sections

}

def captionFor(name, cutoff_map) {
    if (name == 'Optic') {
        return "Tumor classifier predictions by Optic model. Entries with PASS were accepted by the hierarchical call (score above the cutoff for that depth, 0.5 or 0.8)"
    }
    def cutoff = cutoff_map[name]
    cutoff != null
        ? "Tumor classifier predictions by ${name} model. Entries with PASS are above score cutoff ${cutoff}"
        : "Tumor classifier predictions by ${name} model."
}
