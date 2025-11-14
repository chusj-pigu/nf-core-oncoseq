/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_TEXT                                   } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_CNV_SECTION          } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_SV_SECTION           } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_FUSION_SECTION       } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_TARGETS_SECTION      } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_CNV            } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_SV             } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_FUSION         } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_CNLOH          } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_FOCAL          } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_TARGETS        } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE_COLNAMES as QUARTO_FUSION_TABLES } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE_COLNAMES as QUARTO_SV_TABLES     } from '../../../modules/local/quarto/main.nf'
include { modifyMetaId                                  } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'


workflow FIGENO_REPORT {

    take:
    ch_circos_figure
    ch_bin_sizes          // from params qdnaseq_binsize and subchrom_binsize
    ch_sv_figures
    ch_fusion_figures
    ch_targets_figures
    ch_sv_tables
    ch_fusion_tables
    ch_subchrom_figure
    ch_subchrom_focal

    main:

    // channels for type of figures:
    ch_qdnaseq = Channel.of("CNV")
    ch_subchrom = Channel.of("CNLOH")
    ch_focal = Channel.of("CNLOH-Focal")
    ch_fusion = Channel.of("Fusions")
    ch_sv = Channel.of("Structural-Variants")
    ch_targets = Channel.of("Important-Genes")
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CIRCOS FIGURE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Channel for figure subtitle

    ch_subtitle_cnv = ch_bin_sizes
        .filter { meta, _value ->
        meta == "qDNAseq" }
        .map { meta, value ->
        "Circos plot showing ${meta} calls with bin size of ${value} kb" }

    ch_in_circos = ch_circos_figure
        .combine(ch_qdnaseq)
        .combine(ch_subtitle_cnv)

    ch_circos_files = ch_in_circos
        .map { meta, file, type, text ->
            CreateFigureCNVInput(meta, file, type, text)
        }

    QUARTO_FIGURE_CNV(
        ch_circos_files
    )

    ch_section_cnv = QUARTO_FIGURE_CNV.out.quarto_figure
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CNLOH FIGURES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    if (params.realtime == null || params.realtime == 72) {

        ch_subtitle_cnloh = ch_bin_sizes
            .filter { meta, _value ->
            meta == "Subchrom" }
            .map { meta, value ->
            "CNLOH calls by ${meta} with a bin size of ${value} kb" }

        ch_cnloh_in = ch_subchrom_figure
            .combine(ch_subchrom)
            .combine(ch_subtitle_cnloh)

        ch_subchrom_files = ch_cnloh_in
            .map { meta, file, type, text ->
            CreateFigureCNVInput(meta, file, type, text)
            }

        QUARTO_FIGURE_CNLOH(
            ch_subchrom_files
        )

        ch_subtitle_focal = ch_bin_sizes
            .filter { meta, _value ->
            meta == "Subchrom" }
            .map { meta, value ->
            "Gene Focal view by ${meta} with a bin size of ${value} kb" }

        ch_focal_in = ch_subchrom_focal
            .combine(ch_focal)
            .combine(ch_subtitle_focal)

        ch_focal_files = ch_focal_in
            .map { meta, file, type, text ->
            CreateFigureCNVInput(meta, file, type, text)
            }

        QUARTO_FIGURE_FOCAL(
            ch_focal_files
        )

        ch_section_cnv = QUARTO_FIGURE_CNV.out.quarto_figure
            .mix(QUARTO_FIGURE_CNLOH.out.quarto_figure)
            .mix(QUARTO_FIGURE_FOCAL.out.quarto_figure)
            .groupTuple()
            .map { id, section, filePaths ->
                [id, section[0], filePaths]
            }
    } else {
        ch_section_cnv = QUARTO_FIGURE_CNV.out.quarto_figure
    }

    QUARTO_CNV_SECTION(
        ch_section_cnv,
        "CNV Plots"
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SV SECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_sv_table_in = ch_sv_tables
        .combine(ch_sv)
        .map { meta, table, type ->
            def type_lower = type.toLowerCase()
            def new_meta = modifyMetaId(meta, 'replace', '_sv', '', '')
            def caption = "Summary of ${type_lower} detected"
            def col_names = "CHR, GENE, TYPE, SUPPORT, START, END, LEN"
            def section = type_lower
            def process = "${type_lower}-stats-${meta.id}"
            tuple(new_meta, table, caption, col_names, section, process)
        }

    QUARTO_SV_TABLES(ch_sv_table_in)

    // Process tables to extract necessary information to put in legend

    ch_sv_stats = ch_sv_tables
        .flatMap { meta, table ->
            def tuples = []
            table.readLines().each { line ->
                def cols = line.split('\t')
                def gene = cols[1]
                def support = cols[3]
                def type = cols[2]
                def new_meta = meta.id + '_' + gene
                tuples << tuple(new_meta, support, type)
            }
            return tuples
        }

    ch_sv_files = ch_sv_figures
        .combine(ch_sv)
        .transpose()
        .map { meta, figure, type ->
            def meta_refined = figure.baseName
            tuple(meta_refined, meta, figure, type)
        }
        .join(ch_sv_stats)
        .map { _meta, old_meta, file, type, support, sv ->
        CreateSVInput(old_meta, file, type, support, sv)
        }

    QUARTO_FIGURE_SV(
        ch_sv_files
        )

    ch_section_sv = QUARTO_SV_TABLES.out.quarto_table
        .mix(QUARTO_FIGURE_SV.out.quarto_figure)

    ch_section_sv = ch_section_sv
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    QUARTO_SV_SECTION(
        ch_section_sv,
        "Other structural variants called by Sniffles2 with > 4 reads support and at least one annotation of high or moderate impact by SnpEff"
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUSIONS SECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_fusion_table_in = ch_fusion_tables
        .combine(ch_fusion)
        .map { meta, table, type ->
            def type_lower = type.toLowerCase()
            def new_meta = modifyMetaId(meta, 'replace', '_sv', '', '')
            def caption = "Summary of ${type_lower} detected"
            def col_names = "FUSION, CHR1, BREAKPOINT1, CHR2, BREAKPOINT2, TYPE, DIRECTION, SUPPORT"
            def section = type_lower
            def process = "${type_lower}-stats-${meta.id}"
            tuple(new_meta, table, caption, col_names, section, process)
        }

    QUARTO_FUSION_TABLES(ch_fusion_table_in)

    ch_fusion_stats = ch_fusion_tables
        .flatMap { meta, table ->
            def tuples = []
            table.readLines().each { line ->
                def cols = line.split('\t')
                def gene = cols[0]
                def support = cols[7]
                def type = cols[5].replaceAll('_', ' ')
                def new_meta = meta.id + '_' + gene
                tuples << tuple(new_meta, support, type)
            }
            return tuples
        }

    ch_fusion_files = ch_fusion_figures
        .combine(ch_fusion)
        .transpose()
        .map { meta, figure, type ->
            def meta_refined = figure.baseName
            tuple(meta_refined, meta, figure, type)
        }
        .join(ch_fusion_stats)
        .map { _meta, old_meta, file, type, support, sv ->
        CreateSVInput(old_meta, file, type, support, sv)
        }

    QUARTO_FIGURE_FUSION(
        ch_fusion_files
        )

    ch_section_fusion = QUARTO_FUSION_TABLES.out.quarto_table
        .mix(QUARTO_FIGURE_FUSION.out.quarto_figure)

    ch_section_fusion = ch_section_fusion
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    QUARTO_FUSION_SECTION(
       ch_section_fusion,
        "Gene fusions called by Sniffles2 with > 4 reads support and at least one annotation of high or moderate impact by SnpEff"
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENES OF INTEREST SECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_targets_in = ch_targets_figures
        .combine(ch_targets)

    ch_target_files = ch_targets_in
        .transpose()
        .map { meta, file, type ->
            def support_placeholder = ""
            def sv_placeholder = ""
            tuple(meta, file, type, support_placeholder, sv_placeholder) }
        .map { meta, file, type, support_placeholder, sv_placeholder ->
            CreateSVInput(meta, file, type, support_placeholder, sv_placeholder)
            }

    QUARTO_FIGURE_TARGETS(
        ch_target_files
        )

    ch_section_targets = QUARTO_FIGURE_TARGETS.out.quarto_figure

    ch_section_targets = ch_section_targets
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    QUARTO_TARGETS_SECTION (
       ch_section_targets,
        "Genes of interest shown in figeno WITHOUT any Sniffles2 calls with > 4 reads support and at least one annotation of high or moderate impact by SnpEff. SV with length > 1mb is not included in figeno plots."
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_sections = QUARTO_CNV_SECTION.out.quarto_section
        .mix(QUARTO_FUSION_SECTION.out.quarto_section)
        .mix(QUARTO_SV_SECTION.out.quarto_section)
        .mix(QUARTO_TARGETS_SECTION.out.quarto_section)

    emit:
    sections = ch_sections
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def CreateFigureCNVInput(meta, file, type, text) {
    // Transform chrom into two new variables
    def caption = text
    def section = "CNV"
    def process = "${type}-call-${meta.id}"

    // Return a new tuple with the additional variables
    return [meta, file, caption, section, process ]
}

def CreateSVInput(meta, file, type, support, sv) {
    def new_meta = modifyMetaId(meta, 'replace', '_sv', '', '')
    def sv_type = file.name.toString().replaceAll("${meta.id}_", "")
    def sv_type_noext = sv_type.replaceAll(".png", "")
    def type_red = type.toLowerCase().replaceFirst(/s$/, '').replace('-', ' ')

    // Transform chrom into two new variables
    def caption = type_red == "important gene" ?
        "Figeno Plot showing ${sv_type_noext}" :
        "Figeno Plot showing ${sv} in ${sv_type_noext} with ${support} reads support"
    def section = "${type}"
    def process = "figeno-${type}-${new_meta.id}-${sv_type_noext}"

    // Return a new tuple with the additional variables
    return [new_meta, file, caption, section, process ]
}
