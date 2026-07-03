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
    flatCalls
        .sort { a, b -> scoreField(a) <=> scoreField(b) }  // ascending first
        .reverse()                                          // then reverse to get descending
        .take(5)
}

// Single shared channel transformation
def makeCallChannel(ch, typeName, modelConfig, cutoffMap) {
    def config = modelConfig[typeName]
    ch
        .filter { meta, json, type -> type.contains(typeName) }
        .splitJson(path: config.path)
        .groupTuple()
        .map { meta, calls, type ->
            def flatCalls  = calls.flatten()
            def cutoff     = cutoffMap[type.flatten().first()]
            def scoreField = config.score
            def topCalls   = applyFilter(flatCalls, cutoff, scoreField)

            def rows = topCalls.collect { call ->
                def f      = config.fields(call)
                def status = scoreField(call) >= cutoff ? "PASS" : "FAIL"
                [meta.id, f[0], f[1], f[2], f[3], status, type.flatten().first()]
            }
            // Return as single tuple with all rows — don't flatMap yet
            tuple(meta.id, type.flatten().first(), rows)
        }
        .collectFile { meta_id, type_name, rows ->
            def header  = "Family\tClass\tClass name\tScore\tCutoff\n"
            def type_corr = type_name.replace(' ', '-')
            def content   = header + rows.collect { r ->
                "${r[1]}\t${r[2]}\t${r[3]}\t${r[4]}\t${r[5]}"
            }.join('\n') + '\n'
            return [ "${meta_id}_score_${type_corr}.tsv", content ]
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
        ch_methylation_text = ch_methylation_plot
        .map { meta, _plot, _type -> meta }   // extract only meta
        .unique { it.id }                      // deduplicate by sample id
        .map { meta ->
            def section = "Methylation"
            def process = "Methylation-text-${meta.id}"
            tuple(meta, "Tumor classifier predictions based on DNA methylation profiles. Top 5 scoring predictions are shown.", section, process)
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
        'Sturgeon General' : 0.8,
        'Tucan'            : 0.7,
        'Marlin'           : 0.8,
        'CrossNN Caper'    : 0.2,
        'CrossNN PanCancer': 0.15
    ]

    // Configuration map: type -> [json_path, score_field, cutoff, parser_closure]
    def modelConfigClass = [
        'Marlin': [
            path  : 'inference.top_probabilities',
            score : { call -> call.probability },
            fields: { call -> [
                call.lineage.toString().replace(', ', '-'),
                call.group.toString().replace(', ', '-'),
                call.class_name.toString().replace(', ', '-'),
                call.probability
            ]},
                lineage : [
                    path  : 'inference.lineage_probabilities',  // different json path
                    score : { call -> call.probability },
                    label : 'Lineage Probability',
                    fields: { call -> [
                        call.label.toString().replace(', ', '-'),
                        call.probability
                    ]}
                ],
                group : [
                    path  : 'inference.group_probabilities',  // different json path
                    score : { call -> call.probability },
                    label : 'Group Probability',
                    fields: { call -> [
                        call.label.toString().replace(', ', '-'),
                        call.probability
                    ]}
            ],
        ],
        'Tucan': [
            path  : 'tucan.averaged_scores',
            score : { call -> call.score },
            fields: { call -> [
                call.family.toString(),
                call.class.toString(),
                call.name.toString(),
                call.score
            ]},
                lineage : [
                    path  : 'tucan.family_scores',  // different json path
                    score : { call -> call.score },
                    label : 'Family Probability',
                    fields: { call -> [
                        call.label.toString(),
                        call.score
                ]}
            ],
        ],
        'CrossNN PanCancer': [
            path  : 'crossnn.votes',
            score : { call -> call.score },
            fields: { call -> [
                call['Methylation.Class.Family'].toString().replace(', ', '_'),
                call['class'].toString().replace(', ', '-'),
                call['Methylation.Class.Name'].toString().replace(', ', '-'),
                call.score
            ]},
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
        ]
    ]

    def modelConfigFamily = [
        'Marlin': [
            path  : 'inference.lineage_probabilities',  // different json path
            score : { call -> call.probability },
            label : 'Lineage Probability',
            fields: { call -> [
                call.label.toString().replace(', ', '-'),
                call.probability
            ]},
        ],
        'Tucan': [
            path  : 'tucan.family_scores',  // different json path
            score : { call -> call.score },
            label : 'Family Probability',
            fields: { call -> [
                call.label.toString(),
                call.score
            ]}
        ]
    ]

    def modelConfigGroup = [
        'Marlin': [
            path  : 'inference.group_probabilities',  // different json path
            score : { call -> call.probability },
            label : 'Group Probability',
            fields: { call -> [
                call.label.toString().replace(', ', '-'),
                call.probability
            ]},
        ]
    ]

    ch_marlin_call_top_lineage = ch_methylation_call
        .filter { meta, json, type -> type.contains('Marlin') }
        .splitJson(path:   modelConfigClass['Marlin'].path)
        .map { meta, calls, type ->
            def unique_lineage = calls.lineage
            tuple(meta, unique_lineage, type)
        }
        .groupTuple()
        .map { meta, lineages, type ->
            def unique_lineage = lineages.unique()
            tuple(meta, unique_lineage)
        }
        .transpose()

    ch_marlin_call_top_group = ch_methylation_call
        .filter { meta, json, type -> type.contains('Marlin') }
        .splitJson(path:   modelConfigClass['Marlin'].path)
        .map { meta, calls, type ->
            def unique_group = calls.group
            tuple(meta, unique_group, type)
        }
        .groupTuple()
        .map { meta, groups, type ->
            def unique_group = groups.unique()
            tuple(meta, unique_group)
        }
        .transpose()

    ch_marlin_call_lineage = ch_methylation_call
        .filter { meta, json, type -> type.contains('Marlin') }
        .splitJson(path:   modelConfigFamily['Marlin'].path)
        .map { meta, calls, type ->
            def unique_lineage = calls.label
            def score = calls.probability
            tuple(meta, unique_lineage, [lineage: score])
        }
        .join(ch_marlin_call_top_lineage, by: [0,1])

    ch_marlin_call_group = ch_methylation_call
        .filter { meta, json, type -> type.contains('Marlin') }
        .splitJson(path:   modelConfigGroup['Marlin'].path)
        .map { meta, calls, type ->
            def unique_group = calls.label
            def score = calls.probability
            tuple(meta, unique_group, [group: score])
        }
        .join(ch_marlin_call_top_group, by: [0,1])

    ch_marlin_call_class = ch_methylation_call
        .filter { meta, json, type -> type.contains('Marlin') }
        .splitJson(path:   modelConfigClass['Marlin'].path)
        .map { meta, calls, type ->
            def unique_class = calls.class_name
            def unique_lineage = calls.lineage
            def unique_group = calls.group
            def score = calls.probability
            tuple(meta, unique_lineage, unique_group, unique_class, [classes: score], type)
        }
        .combine(ch_marlin_call_lineage, by:[0,1])
        .map { meta, lineage, group, classes, class_score, type, lineage_score ->
            tuple(meta, group, lineage, lineage_score, classes, class_score, type)
        }
        .combine(ch_marlin_call_group, by:[0,1])
        .map { meta, group, lineage, lineage_score, classes, class_score, type, group_score ->
            tuple(meta, lineage, lineage_score, group, group_score, classes, class_score, type)
        }
        .groupTuple()

    ch_marlin_table = ch_marlin_call_class
        .map { meta, lineages, lineage_scores, groups, group_scores, classes, class_scores, type ->
            def typeStr = type.flatten().first()

            // Zip class-level data together and take top 5 by class score
            def rows = [lineages, lineage_scores, groups, group_scores, classes, class_scores]
                .transpose()
                .sort { a, b -> b[5]['classes'] <=> a[5]['classes'] }
                .take(5)
                .collect { lineage, lin_score, group, grp_score, cls, cls_score ->
                    def status = cls_score['classes'] >= cutoffMap['Marlin'] ? "PASS" : "FAIL"
                    [lineage, lin_score['lineage'], group, grp_score['group'], cls, cls_score['classes'], status]
                }
            tuple(meta.id, typeStr, rows)
        }
        .collectFile { meta_id, type_name, rows ->
            def header  = "Lineage\tLineage_Score\tGroup\tGroup_Score\tClass\tClass_Score\tCutoff\n"
            def content = rows.collect { r ->
                "${r[0]}\t${r[1]}\t${r[2]}\t${r[3]}\t${r[4]}\t${r[5]}\t${r[6]}"
            }.join('\n') + '\n'
            ["${meta_id}_score_Marlin.tsv", header + content]
        }

    // Tucan

    ch_tucan_top_family = ch_methylation_call
        .filter { meta, json, type -> type.contains('Tucan') }
        .splitJson(path:   modelConfigClass['Tucan'].path)
        .take(5)
        .map { meta, calls, type ->
            def unique_family = calls.family
            tuple(meta, unique_family, type)
        }
        .groupTuple()
        .map { meta, families, _type ->
            def unique_family = families.unique()
            tuple(meta, unique_family)
        }
        .transpose()

    ch_tucan_call_family = ch_methylation_call
        .filter { meta, json, type -> type.contains('Tucan') }
        .splitJson(path:   modelConfigFamily['Tucan'].path)
        .map { meta, calls, type ->
            def unique_family = calls.family
            def score = calls.score
            tuple(meta, unique_family, [family: score])
        }
        .join(ch_tucan_top_family, by: [0,1])

    ch_tucan_call_class = ch_methylation_call
        .filter { meta, json, type -> type.contains('Tucan') }
        .splitJson(path:   modelConfigClass['Tucan'].path)
        .map { meta, calls, type ->
            def unique_class = calls.name
            def unique_family = calls.family
            def score = calls.score
            tuple(meta, unique_family, unique_class, [classes: score], type)
        }
        .take(5)
        .combine(ch_tucan_call_family, by:[0,1])
        .map { meta, family, classes, class_score, type, family_score ->
            tuple(meta, family, family_score, classes, class_score, type)
        }
        .groupTuple()

    ch_tucan_table = ch_tucan_call_class
        .map { meta, families, family_scores, classes, class_scores, type ->
            def typeStr = type.flatten().first()

            def rows = [families, family_scores, classes, class_scores]
                .transpose()
                .sort { a, b -> b[3]['classes'] <=> a[3]['classes'] }
                .take(5)
                .collect { family, fam_score, cls, cls_score ->
                    def status = cls_score['classes'] >= cutoffMap['Tucan'] ? "PASS" : "FAIL"
                    [family, fam_score['family'], cls, cls_score['classes'], status]
                }
            tuple(meta.id, typeStr, rows)
        }
        .collectFile { meta_id, type_name, rows ->
            def header  = "Family\tFamily_Score\tClass\tClass_Score\tCutoff\n"
            def content = rows.collect { r ->
                "${r[0]}\t${r[1]}\t${r[2]}\t${r[3]}\t${r[4]}"
            }.join('\n') + '\n'
            ["${meta_id}_score_Tucan.tsv", header + content]
        }

    ch_crossnn_call_pancancer = makeCallChannel(ch_methylation_call, 'CrossNN PanCancer', modelConfigClass, cutoffMap)
    ch_crossnn_call_caper    = makeCallChannel(ch_methylation_call, 'CrossNN Caper', modelConfigClass, cutoffMap)
    ch_sturgeon_call         = makeCallChannel(ch_methylation_call, 'Sturgeon General', modelConfigClass, cutoffMap)

    ch_tables = ch_marlin_table
        .mix(ch_tucan_table)
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
