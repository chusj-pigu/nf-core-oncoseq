/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_TEXT           } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION        } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE         } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE_COLNAMES } from '../../../modules/local/quarto/main.nf'
include { modifyMetaId          } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'


workflow ADAPTIVE_REPORT {

    take:
    ch_coverage_csv
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
    ch_cov_table = ch_coverage_csv
        .map {meta, table ->
            def lines = table.readLines()
            def panel_cov = lines[1].tokenize(',')[1]
            def bg_cov = lines[1].tokenize(',')[2]
            tuple(meta.id, panel_cov, bg_cov) }
        .collectFile { table ->
            def content = [table[1] + '\t' + table[2] ].join('\n')
            return [ "${table[0]}_adaptive_coverage.tsv", content + '\n' ]
        }
        .map { table ->
            def meta = table.name.replace('_adaptive_coverage.tsv', '')
            tuple(id:meta,table)
            }                                                                   // Add back meta to the new table

    ch_quarto_table = ch_cov_table
        .map { meta, table ->
            def caption = "On target and off target (background) coverage"
            def col_names = "Coverage on Target (X), Coverage off Target (X)"
            def section = "Cov_QC"
            def process = "cov_qc-${meta.id}"
            tuple(meta, table, caption, col_names, section, process)
        }

    QUARTO_TABLE_COLNAMES(
        ch_quarto_table
    )

    ch_section_inputs = QUARTO_TABLE_COLNAMES.out.quarto_table


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
            def process = "coverage-plot-${meta.id}"

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
            [id, section[0], filePaths,
                'Coverage analysis with mosdepth and custom coverage plot']
        }

    QUARTO_SECTION(
        ch_section_inputs
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    emit:
    sections = QUARTO_SECTION.out.quarto_section

}
