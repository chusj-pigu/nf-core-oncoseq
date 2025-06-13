/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_TEXT } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE as QUARTO_TABLE_BACKGROUND } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE as QUARTO_TABLE_PANEL } from '../../../modules/local/quarto/main.nf'
include { modifyMetaId } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'


workflow COVERAGE_REPORT {

    take:
    ch_cramino_stats
    ch_coverage_plot

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
    CONSTRUCT QUARTO TABLE FOR CRAMINO STATS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_cramino_stats_panel = ch_cramino_stats
        .filter { meta, file -> meta.id.contains('_panel') }
        .map { meta, stats ->
            // Restore original sample ID for output naming
            def meta_restore = modifyMetaId(meta, 'replace', '_panel', '', '')
            tuple(meta_restore, stats)
        }

    ch_cramino_stats_background = ch_cramino_stats
        .filter { meta, file -> meta.id.contains('_background') }
        .map { meta, stats ->
            // Restore original sample ID for output naming
            def meta_restore = modifyMetaId(meta, 'replace', '_background', '', '')
            tuple(meta_restore, stats)
        }

    ch_cramino_stats_panel = ch_cramino_stats_panel
        .map { tuple ->
        // Extract the existing values from the tuple
        def (meta, file) = tuple

        // Transform chrom into two new variables
        def caption = "Panel Summary Stats of Cramino Analysis."
        def col_names = "T"
        def section = "Coverage"
        def process = "crammino-panel-stats-${meta.id}"

        // Return a new tuple with the additional variables
        return [meta, file, caption, col_names, section, process ]
    }

    ch_cramino_stats_background = ch_cramino_stats_background
        .map { tuple ->
            // Extract the existing values from the tuple
            def (meta, file) = tuple

            // Transform chrom into two new variables
            def caption = "Background Summary Stats of Cramino Analysis."
            def col_names = "T"
            def section = "Coverage"
            def process = "crammino-background-stats-${meta.id}"

            // Return a new tuple with the additional variables
            return [meta, file, caption, col_names, section, process ]
        }

    QUARTO_TABLE_PANEL(
        ch_cramino_stats_panel
    )

    QUARTO_TABLE_BACKGROUND(
        ch_cramino_stats_background
    )

    ch_section_inputs = QUARTO_TABLE_PANEL.out.quarto_table
        .mix(QUARTO_TABLE_BACKGROUND.out.quarto_table)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSTRUCT QUARTO FIGURE FOR MPGI COVERAGE PLOT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_plot_files = ch_coverage_plot
        .map { tuple ->
            // Extract the existing values from the tuple
            def (meta, file) = tuple

            // Transform chrom into two new variables
            def caption = "Coverage Plot for ${meta.id}."
            def section = "Coverage"
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
        "Coverage analysis with cramino and custom coverage plot."
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    emit:
    sections = QUARTO_SECTION.out.quarto_section

}