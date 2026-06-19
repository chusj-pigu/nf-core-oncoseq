/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_TEXT                                   } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_CNV_SECTION          } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_DELLY_SECTION        } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_SV_SECTION           } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_FUSION_SECTION       } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_TARGETS_SECTION      } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_SNP_SECTION          } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_QDNASEQ        } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_DELLY          } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_SV             } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_FUSION         } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_CNLOH          } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_FOCAL          } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_TARGETS        } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE_TABS as QUARTO_FUSION_TABLES     } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE_TABS as QUARTO_SV_TABLES         } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE_TABS as QUARTO_SNP_TABLES        } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TEXT as QUARTO_TEXT_SV                 } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TEXT as QUARTO_TEXT_FUSION             } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TEXT as QUARTO_TEXT_SNP                } from '../../../modules/local/quarto/main.nf'
include { modifyMetaId                                  } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'


workflow FIGENO_REPORT {

    take:
    ch_circos_figure
    ch_delly_figure
    ch_bin_sizes          // from params qdnaseq_binsize and subchrom_binsize
    ch_sv_figures
    ch_fusion_figures
    ch_targets_figures
    ch_sv_tables
    ch_fusion_tables
    ch_subchrom_figure
    ch_subchrom_focal
    ch_snp_table

    main:

    // channels for type of figures:
    ch_qdnaseq = channel.of("CNV")
    ch_delly = channel.of("PanChr")
    ch_subchrom = channel.of("CNLOH")
    ch_focal = channel.of("CNLOH-Focal")
    ch_fusion = channel.of("Fusions")
    ch_sv = channel.of("Structural-Variants")
    ch_targets = channel.of("Important-Genes")
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CIRCOS FIGURE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // channel for figure subtitle

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

    QUARTO_FIGURE_QDNASEQ(
        ch_circos_files
    )

    ch_section_cnv = QUARTO_FIGURE_QDNASEQ.out.quarto_figure
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DELLY FIGURE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // channel for figure subtitle

    ch_in_panchr = ch_delly_figure
        .transpose()
        .combine(ch_delly)
        .map { meta, figure, type ->
            // Add chromosome to meta to avoid file collision
            def chr = (figure.baseName =~ /.*_(chr\w+)$/)[0][1]
            def new_meta = [
                id      : meta.id + '_' + chr,
                sample  : meta.id     // store original once
            ]
            def subtitle = chr
            tuple(new_meta, figure, type, chr)
            }

    ch_delly_files = ch_in_panchr
        .map { meta, file, type, text ->
            // Transform chrom into two new variables
            def caption = text
            def section = "Pan-Chromosome"
            def process = "${type}-call-${meta.id}"

            // Return a new tuple with the additional variables
            return [meta, file, caption, section, process ]
        }

    QUARTO_FIGURE_DELLY(
        ch_delly_files
    )

    ch_section_delly = QUARTO_FIGURE_DELLY.out.quarto_figure
        .map { meta, section, figure ->
            tuple(id:meta.sample, section, figure) }
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths,
                "Pan chromosome plots showing Delly calls with bin size of ${params.delly_bin_size} kb"]
        }

    QUARTO_DELLY_SECTION(
        ch_section_delly
    )




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

        ch_section_cnv = QUARTO_FIGURE_QDNASEQ.out.quarto_figure
            .mix(QUARTO_FIGURE_CNLOH.out.quarto_figure)
            .mix(QUARTO_FIGURE_FOCAL.out.quarto_figure)
            .groupTuple()
            .map { id, section, filePaths ->
                [id, section[0], filePaths, 'CNV Plots']
            }
    } else {
        ch_section_cnv = QUARTO_FIGURE_QDNASEQ.out.quarto_figure
            .map { id, section, filePaths ->
                [id, section, filePaths, 'CNV Plots']
            }
    }

    QUARTO_CNV_SECTION(
        ch_section_cnv
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SV SECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // If table is empty, skip quarto table

    ch_sv_tables = ch_sv_tables
        .map { meta, table ->
            def type = table.name.contains('sniffles') ? 'sniffles' : 'severus'
            tuple(meta, table, type) }
        .branch { meta, table, type ->
            pos: table.size() > 0
            empty: true }

    ch_sv_table_in = ch_sv_tables.pos
        .groupTuple()
        .combine(ch_sv)
        .map { meta, tables, type, variant ->
            def variant_lower = variant.toLowerCase()
            def tab = type
            def caption = type.collect { types ->
                "Summary of ${variant_lower} detected by ${types}"
            }
            def col_names = "CHR, GENE, TYPE, START, END, LEN, SUPPORT"
            def section = variant_lower
            def process = "${variant_lower}-stats-${meta.id}"
            tuple(meta, section, process, tables, tab, caption, col_names)}

    QUARTO_SV_TABLES(ch_sv_table_in)

    // Process tables to extract necessary information to put in legend

    ch_sv_stats = ch_sv_tables.pos
        .flatMap { meta, table, program ->
            def tuples = []
            table.readLines().each { line ->
                def cols = line.split('\t')
                def gene = cols[1]
                def support = cols[6]
                def type = cols[2]
                def new_meta = meta.id + '_' + gene
                tuples << tuple(new_meta, support, type, program)
            }
            return tuples
        }
        .groupTuple()

    ch_sv_files = ch_sv_figures
        .combine(ch_sv)
        .transpose()
        .map { meta, figure, type ->
            def meta_refined = figure.baseName
            tuple(meta_refined, meta, figure, type)
        }
        .join(ch_sv_stats)
        .map { _meta, old_meta, file, type, support, sv, program ->
        CreateSVInput(old_meta, file, type, support, sv, program)
        }

    QUARTO_FIGURE_SV(
        ch_sv_files
        )

    // Empty table :
    ch_empty = ch_sv_tables.empty
        .combine(ch_sv)
        .groupTuple()
        .map { meta, _tables, types, variants ->
            def variant_uniqe = variants[0]
            def variant_lower = variant_uniqe.toLowerCase()
            def types_empty = types.size() > 1 ? "${types[0]} or ${types[1]}" : "${types[0]}"
            def text = "No Structural variants called by ${types_empty} remaining after applying filters"  // Placeholder text for empty table
            def section = variant_lower
            def type_multi = types_empty.replace(' ', '-')
            def process = "${variant_lower}-${type_multi}-empty-${meta.id}"
            tuple(meta, text, section, process)
        }

    QUARTO_TEXT_SV(
        ch_empty
    )

    ch_section_sv = QUARTO_SV_TABLES.out.quarto_table
        .mix(QUARTO_FIGURE_SV.out.quarto_figure)
        .mix(QUARTO_TEXT_SV.out.quarto_text)

    ch_section_sv = ch_section_sv
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths,
                'Other structural variants with > 4 reads support, at least one annotation of high or moderate impact, and affecting at least one gene in the panel']
        }

    QUARTO_SV_SECTION(
        ch_section_sv
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUSIONS SECTION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_fusion_tables = ch_fusion_tables
        .map { meta, table ->
            def type = table.name.contains('sniffles') ? 'sniffles' : 'severus'
            tuple(meta, table, type)}
        .branch { meta, table, type ->
            pos: table.size() > 0
            empty: true }

    ch_fusion_table_in = ch_fusion_tables.pos
        .groupTuple()
        .combine(ch_fusion)
        .map { meta, tables, type, variant ->
            def variant_lower = variant.toLowerCase()
            def tab = type
            def caption = type.collect { types ->
                "Summary of ${variant_lower} detected by ${types}"
            }
            def col_names = "FUSION, CHR1, BREAKPOINT1, CHR2, BREAKPOINT2, TYPE, DIRECTION, SUPPORT"
            def section = variant_lower
            def process = "${variant_lower}-stats-${meta.id}"
            tuple(meta, section, process, tables, tab, caption, col_names)}

    QUARTO_FUSION_TABLES(ch_fusion_table_in)

    ch_fusion_stats = ch_fusion_tables.pos
        .flatMap { meta, table, program ->
            def tuples = []
            table.readLines().each { line ->
                def cols = line.split('\t')
                def gene = cols[0]
                def support = cols[7]
                def type = cols[5].replaceAll('_', ' ')
                def type_norm = program == "severus" ? type + 'fusion' : type
                def new_meta = meta.id + '_' + gene
                tuples << tuple(new_meta, support, type_norm, program)
            }
            return tuples
        }
        .groupTuple()

    ch_fusion_files = ch_fusion_figures
        .combine(ch_fusion)
        .transpose()
        .map { meta, figure, type ->
            def meta_refined = figure.baseName
            tuple(meta_refined, meta, figure, type)
        }
        .join(ch_fusion_stats)
        .map { _meta, old_meta, file, type, support, sv, program ->
        CreateSVInput(old_meta, file, type, support, sv, program)
        }

    QUARTO_FIGURE_FUSION(
        ch_fusion_files
        )

    // Empty table :

    ch_empty_fusion = ch_fusion_tables.empty
        .combine(ch_fusion)
        .groupTuple()
        .map { meta, _tables, types, variants ->
            def variant_uniqe = variants[0]
            def variant_lower = variant_uniqe.toLowerCase()
            def types_empty = types.size() > 1 ? "${types[0]} or ${types[1]}" : "${types[0]}"
            def text = "No gene fusions called by ${types_empty} remaining after applying filters"  // Placeholder text for empty table
            def section = variant_lower
            def type_multi = types_empty.replace(' ', '-')
            def process = "${variant_lower}-${type_multi}-empty-${meta.id}"
            tuple(meta, text, section, process)
        }

    QUARTO_TEXT_FUSION(
        ch_empty_fusion
    )

    ch_section_fusion = QUARTO_FUSION_TABLES.out.quarto_table
        .mix(QUARTO_FIGURE_FUSION.out.quarto_figure)
        .mix(QUARTO_TEXT_FUSION.out.quarto_text)

    ch_section_fusion = ch_section_fusion
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths,
                'Gene fusions with > 4 reads support, at least one annotation of high or moderate impact, and affecting at least one gene in the panel']
        }

    QUARTO_FUSION_SECTION(
       ch_section_fusion
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
            def program_placeholder = ""
            tuple(meta, file, type, support_placeholder, sv_placeholder, program_placeholder) }
        .map { meta, file, type, support_placeholder, sv_placeholder, program_placeholder ->
            CreateSVInput(meta, file, type, support_placeholder, sv_placeholder, program_placeholder)
            }

    QUARTO_FIGURE_TARGETS(
        ch_target_files
        )

    ch_section_targets = QUARTO_FIGURE_TARGETS.out.quarto_figure

    ch_section_targets = ch_section_targets
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths,
                'Genes of interest shown in figeno WITHOUT any calls with > 4 reads support and at least one annotation of high or moderate impact. SV with length > 1mb is not included in figeno plots.']
        }

    QUARTO_TARGETS_SECTION (
       ch_section_targets
    )


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FILTERED SNP TABLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_snp_table = ch_snp_table
        .map { meta,table ->
            def type = meta.id.contains('somatic') ? 'ClairS-TO' : 'Clair3'
            def new_meta = modifyMetaId(meta, 'replace', '_somatic_snp_vep', '', '')
            def meta_final = modifyMetaId(new_meta, 'replace', '_germline_snp_vep', '', '')
            tuple(meta_final, table, type)}
        .branch { meta, table, type ->
            pos: table.readLines().size() > 1
            empty: true }

    ch_snp_table_in = ch_snp_table.pos
        .groupTuple()
        .map { meta, tables, type ->
            def tab = type
            def caption = type.collect { types ->
                "SNPs called by $types and filtered for regions included in the panels"
            }
            def col_names = ""
            def section = "SNPs"
            def process = "snp-${type}-stats-${meta.id}"
            tuple(meta, section, process, tables, tab, caption, col_names)
        }

    ch_empty_snp = ch_snp_table.empty
        .groupTuple()
        .map { meta, _tables, types ->
            def types_empty = types.size() > 1 ? "${types[0]} or ${types[1]}" : "${types[0]}"
            def text = "No SNPs called by ${types_empty} remaining after applying filter in targeted regions."
            def section = "SNPs"
            def type_multi = types_empty.replace(' ', '-')
            def process = "snp-${type_multi}-empty-${meta.id}"
            tuple(meta, text, section, process)
        }

    QUARTO_TEXT_SNP(ch_empty_snp)

    QUARTO_SNP_TABLES(ch_snp_table_in)

    ch_section_snp = QUARTO_SNP_TABLES.out.quarto_table
        .mix(QUARTO_TEXT_SNP.out.quarto_text)
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths, "SNPs filtered with EnsemblVep using filters :'$params.filtervep_expression'"]
        }

    QUARTO_SNP_SECTION(ch_section_snp)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_sections = QUARTO_CNV_SECTION.out.quarto_section
        .mix(QUARTO_FUSION_SECTION.out.quarto_section)
        .mix(QUARTO_SV_SECTION.out.quarto_section)
        .mix(QUARTO_TARGETS_SECTION.out.quarto_section)
        .mix(QUARTO_DELLY_SECTION.out.quarto_section)
        .mix(QUARTO_SNP_SECTION.out.quarto_section)

    ch_versions = QUARTO_CNV_SECTION.out.versions
        .mix(QUARTO_FUSION_SECTION.out.versions)
        .mix(QUARTO_SV_SECTION.out.versions)
        .mix(QUARTO_TARGETS_SECTION.out.versions)
        .mix(QUARTO_DELLY_SECTION.out.versions)
        .mix(QUARTO_SNP_SECTION.out.versions)

    emit:
    sections = ch_sections
    versions = ch_versions
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

def CreateSVInput(meta, file, type, support, sv, program) {
    def new_meta = modifyMetaId(meta, 'replace', '_sv', '', '')
    def sv_unique = sv instanceof List ? sv[0] : sv
    def sv_type = file.name.toString().replaceAll("${meta.id}_", "")
    def sv_type_noext = sv_type.replaceAll(".png", "")
    def type_red = type.toLowerCase().replaceFirst(/s$/, '').replace('-', ' ')
    def multi_program_support = program.size() > 1 ? "${support[0]} (${program[0]}) and ${support[1]} (${program[1]})" :
        (support instanceof List ? "${support[0]} (${program})" : "${support} (${program})")
    // Transform chrom into two new variables
    def caption = type_red == "important gene" ?
        "Figeno Plot showing ${sv_type_noext}" :
        "Figeno Plot showing ${sv_unique} in ${sv_type_noext} with ${multi_program_support} reads support"
    def section = "${type}"
    def process = "figeno-${type}-${new_meta.id}-${sv_type_noext}"

    // Return a new tuple with the additional variables
    return [new_meta, file, caption, section, process ]
}

