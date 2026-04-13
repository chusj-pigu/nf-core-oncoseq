/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_SECTION        } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE_TABS    } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TEXT           } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE_TABS     } from '../../../modules/local/quarto/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   DEFINE FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Reusable closure to filter calls by cutoff or fall back to max
def applyFilter(flatCalls, cutoff, scoreField) {
    def topCalls = flatCalls.findAll { scoreField(it) >= cutoff }
    if (!topCalls) {
        topCalls = [ flatCalls.max { scoreField(it) } ]
    }
    return topCalls
}

// Single shared channel transformation
def makeCallChannel(ch, typeName, modelConfig, cutoffMap) {
    def config = modelConfig[typeName]
    ch
        .filter { meta, json, type -> type.contains(typeName) }
        .splitJson(path: config.path)
        .groupTuple()
        .map { meta, calls, type ->
            def flatCalls = calls.flatten()
            def cutoff    = cutoffMap[type.flatten().first()]
            def topCalls  = applyFilter(flatCalls, cutoff, config.score)

            topCalls.collect { call ->
                def f = config.fields(call)
                tuple(meta.id, f[0], f[1], f[2], f[3], type.flatten().first())
            }
        }
        .flatMap { it }
        .collectFile { table ->
            def type_corr = table[5].replace(' ', '-')
            def content   = "${table[1]}\t${table[2]}\t${table[3]}\t${table[4]}\n"
            return [ "${table[0]}_score_${type_corr}.tsv", content ]
        }
}


workflow CLASSIFIER_REPORT {

    take:
    ch_methylation_plot
    ch_methylation_call

    main:

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QUARTO FIGURES with tabs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_methylation_text = ch_methylation_plot
        .map { meta, _plot, _type ->
            tuple(meta, "Tumor classifier predictions based on DNA methylation profiles. If no predictions pass the score cutoff, the highest scoring prediction is shown regardless of cutoff.") }
        .distinct()
        .map { meta, text ->
            def section = "Methylation"
            def process = "Methylation-text-${meta.id}"
            tuple(meta, text, section, process)
        }

    QUARTO_TEXT(
        ch_methylation_text
    )

    ch_quarto_fig = ch_methylation_plot
        .groupTuple()
        .map { meta, plots, types ->

            def tabs = types
            def captions = types.collect { type ->
                "Tumor classifier prediction by ${type} model"
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
        'Marlin'           : 0.8,
        'CrossNN Caper'    : 0.2,
        'CrossNN PanCancer': 0.15,
        'Sturgeon General' : 0.8,
        'Tucan'            : 0.7
    ]

    // Configuration map: type -> [json_path, score_field, cutoff, parser_closure]
    def modelConfig = [
        'Marlin': [
            path  : 'inference.top_probabilities',
            score : { call -> call.probability },
            fields: { call -> [
                call.lineage.toString().replace(', ', '-'),
                call.group.toString().replace(', ', '-'),
                call.class_name.toString().replace(', ', '-'),
                call.probability
            ]}
        ],
        'Tucan': [
            path  : 'tucan.averaged_scores',
            score : { call -> call.score },
            fields: { call -> [
                call.family.toString(),
                call.class.toString(),
                call.name.toString(),
                call.score
            ]}
        ],
        'CrossNN PanCancer': [
            path  : 'crossnn.votes',
            score : { call -> call.score },
            fields: { call -> [
                call['Methylation.Class.Family'].toString().replace(', ', '-'),
                call['class'].toString().replace(', ', '-'),
                call['Methylation.Class.Name'].toString().replace(', ', '-'),
                call.score
            ]}
        ],
        'CrossNN Caper': [
            path  : 'crossnn.votes',
            score : { call -> call.score },
            fields: { call -> [
                call['Methylation.Class.Family'].toString().replace(',', '').replace(' ', '_'),
                call['class'].toString().replace(',', '').replace(' ', '_'),
                call['Methylation.Class.Name'].toString().replace(',', '').replace(' ', '_'),
                call.score
            ]}
        ],
        'Sturgeon General': [
            path  : 'sturgeon.scores',
            score : { call -> call.score },
            fields: { call ->
                def parts     = call['class'].toString().split(' - ')
                def family    = parts.size() > 0 ? parts[0].trim() : 'Unknown'
                def group     = parts.size() > 1 ? parts[1].trim() : 'Unknown'
                def className = parts.size() > 2 ? parts[2].trim() : call['class'].toString()
                [ family, group, className, call.score ]
            }
        ]
    ]

    ch_marlin_call           = makeCallChannel(ch_methylation_call, 'Marlin', modelConfig, cutoffMap)
    ch_crossnn_call_tucan    = makeCallChannel(ch_methylation_call, 'Tucan', modelConfig, cutoffMap)
    ch_crossnn_call_pancancer = makeCallChannel(ch_methylation_call, 'CrossNN PanCancer', modelConfig, cutoffMap)
    ch_crossnn_call_caper    = makeCallChannel(ch_methylation_call, 'CrossNN Caper', modelConfig, cutoffMap)
    ch_sturgeon_call         = makeCallChannel(ch_methylation_call, 'Sturgeon General', modelConfig, cutoffMap)

    ch_tables = ch_marlin_call
        .mix(ch_crossnn_call_tucan)
        .mix(ch_crossnn_call_pancancer)
        .mix(ch_crossnn_call_caper)
        .mix(ch_sturgeon_call)
        .map { table ->
            def sample = table.name.split('_score_')[0].replace('.tsv', '')
            def type = table.name.split('_score_')[1].replace('.tsv', '').replace('-', ' ')
            tuple(id: sample, table, type) }

    ch_tables_in_quarto = ch_tables
        .groupTuple()
        .map { meta, tables, type ->
            def tab = type
            def captions = type.collect { types ->
                def cutoff = cutoffMap[types]
                "Tumor classifier predictions by ${types} model above score cutoff ${cutoff}"
            }
            def colnames = "Family, Class, Class name, Score"
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
