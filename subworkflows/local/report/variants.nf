/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_TEXT } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE } from '../../../modules/local/quarto/main.nf'


workflow QDNASEQ_REPORT {

    take:
    ch_qdnaseq_figure

    main:

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SECTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Collect sections from all analysis steps
    // ch_sections = ch_sections.mix(SUMMARIZE_ANALYSIS.out.ch_section)
    ch_section_inputs = Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QDNASEQ FIGURE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_plot_files = ch_qdnaseq_figure
        .map { tuple ->
            // Extract the existing values from the tuple
            def (meta, file) = tuple

            // Transform chrom into two new variables
            def caption = "CNV Plot for ${meta.id}."
            def section = "CNV"
            def process = "qdnaseq-call-${meta.id}"

            // Return a new tuple with the additional variables
            return [meta, file, caption, section, process ]
        }


    QUARTO_FIGURE(
        ch_plot_files
        )

    ch_section_inputs = ch_section_inputs.mix(QUARTO_FIGURE.out.quarto_figure)

    ch_section_inputs = ch_section_inputs
    .groupTuple()
    .map { id, section, filePaths ->
        [id, section[0], filePaths]
    }

    QUARTO_SECTION(
        ch_section_inputs,
        "QDNASEQ CNV Plot."
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    emit:
    sections = QUARTO_SECTION.out.quarto_section
}