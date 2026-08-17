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

    ch_methylation_text = tumor_type
        .map { meta, type ->
            def section = "Methylation"
            def process = "Methylation-text-${meta.id}"
            def classifier_list = type == "blood" ? "Alma, Lamprey and Marlin" :
                type == "brain" ? "CrossNN Capper, MPACT and Sturgeon" :
                type == "solid" ? "CrossNN PanCancer and Tucan" : "All classifiers"
            tuple(meta, "Tumor classifier predictions based on DNA methylation profiles. Top 5 scoring predictions are shown, and relevant classifiers for ${type} tumor types are ${classifier_list}.", section, process)
        }

    QUARTO_TEXT(
        ch_methylation_text
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
        'Sturgeon General' : 0.8,
        'Tucan'            : 0.7,
        'Marlin'           : 0.8,
        'CrossNN Caper'    : 0.2,
        'CrossNN PanCancer': 0.15,
        'Alma'             : 0.2,
        'Lamprey'          : 0.9,
        'Sturgeon Brainstem' : 0.8,
        'MPACT'            : 0.7,
        'Nanomix'          : 0.8
    ]

    ch_tables_in_quarto = ch_methylation_call
        .groupTuple()
        .map { meta, tables, type ->
            def tab = type
            def captions = type.collect { types ->
                def cutoff = cutoffMap[types]
                "Tumor classifier predictions by ${types} model. Entries with PASS are above score cutoff ${cutoff}"
            }
            def colnames = ""
            def section = "Methylation"
            def process = "methylation-calls-${meta.id}"
            tuple(meta, section, process, tables, tab, captions, colnames)}

    QUARTO_TABLE_TABS(
        ch_tables_in_quarto
    )


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
