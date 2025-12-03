/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_SECTION        } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE         } from '../../../modules/local/quarto/main.nf'
include { CONVERT_PDF_PNG       } from '../../../modules/local/magick/main.nf'
include { QUARTO_TABLE_COLNAMES } from '../../../modules/local/quarto/main.nf'


workflow CLASSIFIER_REPORT {

    take:
    ch_methylation_plot
    ch_methylation_call

    main:

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSTRUCT QUARTO FIGURE AND SECTION FOR MARLIN/STURGEON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_to_convert = ch_methylation_plot
        .filter { _meta, _pdf, type ->
            type == "Sturgeon" }
        .map { meta, file, _type ->
            tuple(meta,file)}
    CONVERT_PDF_PNG(ch_to_convert)

    ch_plots = CONVERT_PDF_PNG.out.png
        .join(ch_methylation_plot)
        .map { meta, png, _pdf, type ->
        tuple(meta,png,type)}
        .mix(ch_methylation_plot.filter { _meta, _pdf, type -> type == "Marlin" })
        .view()

    ch_meth_files = ch_plots
        .map { tuple ->
            // Extract the existing values from the tuple
            def (meta, file, type) = tuple

            // Transform chrom into two new variables
            def caption = "Tumor classifier prediction by ${type} for ${meta.id}"
            def section = "Methylation"
            def process = "meth-plot-${meta.id}"

            // Return a new tuple with the additional variables
            return [meta, file, caption, section, process ]
        }
    QUARTO_FIGURE(
        ch_meth_files
    )

    // Extract Score of the most confident types for both classifiers

    ch_sturgeon_call = ch_methylation_call
        .filter { _meta, _pdf, type ->
            type == "Sturgeon" }
        .splitCsv(elem:1, header:true)
        .map { meta, data, type ->
            // Filter out the 'number_probes' key
            def filtered = data.findAll { k, v -> k != 'number_probes' }

            // Find the entry whose value is the maximum (convert to double if needed)
            def maxEntry = filtered.max { entry -> entry.value.toString().toDouble() }

            // Extract key and value
            def col = maxEntry.key.replaceAll(/\s+/, '')
            def val = maxEntry.value

            tuple(meta.id, col, val, type)
        }
        .collectFile { table ->
            def content = [table[1] + '\t' + table[2] ].join('\n')
            return [ "${table[0]}_score_Sturgeon.tsv", content + '\n' ]
        }
        .map { table ->
            def meta = table.name.replace('_score_Sturgeon.tsv', '')
            def type = "Sturgeon"
            tuple(id:meta,table, type)
        }

    ch_marlin_call = ch_methylation_call
        .filter { _meta, _pdf, type ->
            type == "Marlin" }
        .splitJson(path: 'inference.top_probabilities')
        .groupTuple()
        .map {
            meta, calls, type ->
            def group = calls.flatten().first().group.toString().replace('::','-')
            def class_type = calls.flatten().first().class_name.toString().replace('::','-')
            tuple(meta.id, calls.flatten().first().lineage, group, class_type, calls.flatten().first().probability, type.flatten().first())
        }
        .view()
        .collectFile { table ->
            def content = [ table[1] + '\t' + table[2] + '\t' + table[3] + '\t' + table[4] ].join('\n')
            return [ "${table[0]}_score_Marlin.tsv", content + '\n' ]
        }
        .map { table ->
            def meta = table.name.replace('_score_Marlin.tsv', '')
            def type = "Marlin"
            tuple(id:meta,table, type)
        }

    ch_tables = ch_sturgeon_call
        .mix(ch_marlin_call)

    ch_quarto_table = ch_tables
        .map { meta, table, type ->
            def caption = "Most confident ${type} call"
            def col_names = "Lineage, Group, Class, Score"
            def section = "Methylation"
            def process = "${type}-score-${meta.id}"
            tuple(meta, table, caption, col_names, section, process)
        }

     QUARTO_TABLE_COLNAMES(
        ch_quarto_table
    )

    ch_section_inputs = QUARTO_FIGURE.out.quarto_figure
        .mix(QUARTO_TABLE_COLNAMES.out.quarto_table)

    ch_section_meh = ch_section_inputs
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    QUARTO_SECTION(
        ch_section_meh,
        "Methylation based classification"
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
