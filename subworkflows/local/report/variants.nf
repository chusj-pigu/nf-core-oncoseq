/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_TEXT                             } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_CNV_SECTION    } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_SV_SECTION     } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_FUSION_SECTION } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_CNV      } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_SV       } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_FUSION   } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_CNLOH    } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_FOCAL    } from '../../../modules/local/quarto/main.nf'
include { modifyMetaId                            } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'


workflow FIGENO_REPORT {

    take:
    ch_circos_figure
    ch_sv_figures
    ch_fusion_figures
    ch_subchrom_figure
    ch_subchrom_focal

    main:

    // channels for type of figures:
    ch_qdnaseq = Channel.of("CNV")
    ch_subchrom = Channel.of("CNLOH")
    ch_focal = Channel.of("CNLOH-Focal")
    ch_fusion = Channel.of("Fusion")
    ch_sv = Channel.of("Structural-Variants")
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CIRCOS FIGURE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_in_circos = ch_circos_figure
        .combine(ch_qdnaseq)

    ch_circos_files = ch_in_circos
        .map { meta, file, type ->
        CreateFigureCNVInput(meta, file, type)
        }

    QUARTO_FIGURE_CNV(
        ch_circos_files
    )

    ch_section_qdnaseq = QUARTO_FIGURE_CNV.out.quarto_figure
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CNLOH FIGURES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_cnloh_in = ch_subchrom_figure
        .combine(ch_subchrom)

    ch_subchrom_files = ch_cnloh_in
        .map { meta, file, type ->
        CreateFigureCNVInput(meta, file, type)
        }

    QUARTO_FIGURE_CNLOH(
        ch_subchrom_files
        )

    ch_section_cnloh = QUARTO_FIGURE_CNLOH.out.quarto_figure
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    ch_focal_in = ch_subchrom_focal
        .combine(ch_focal)

    ch_focal_files = ch_focal_in
        .map { meta, file, type ->
        CreateFigureCNVInput(meta, file, type)
        }

    QUARTO_FIGURE_FOCAL(
        ch_focal_files
    )

    ch_section_focal = QUARTO_FIGURE_FOCAL.out.quarto_figure
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    ch_section_cnv = ch_section_qdnaseq
        .mix(ch_section_cnloh)
        .mix(ch_section_focal)

    QUARTO_CNV_SECTION(
        ch_section_cnv,
        "CNV Plots"
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FIGENO FOR SV
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_sv_in = ch_sv_figures
        .combine(ch_sv)

    ch_sv_files = ch_sv_in
        .map { meta, file, type ->
        CreateSVInput(meta, file, type)
        }

    QUARTO_FIGURE_SV(
        ch_sv_files
        )

    ch_section_sv = QUARTO_FIGURE_SV.out.quarto_figure

    ch_section_sv = ch_section_sv
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    QUARTO_SV_SECTION(
        ch_section_sv,
        "Figeno SV plots"
    )

    ch_fusion_in = ch_fusion_figures
        .combine(ch_fusion)

    ch_fusion_files = ch_fusion_in
        .map { meta, file, type ->
            CreateSVInput(meta, file, type)
            }

    QUARTO_FIGURE_FUSION(
        ch_fusion_files
        )

    ch_section_fusion = QUARTO_FIGURE_FUSION.out.quarto_figure

    ch_section_fusion = ch_section_fusion
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    QUARTO_FUSION_SECTION(
       ch_section_fusion,
        "Figeno Gene Fusion plots"
    )

    ch_sections = QUARTO_CNV_SECTION.out.quarto_section
        .mix(QUARTO_FUSION_SECTION.out.quarto_section)
        .mix(QUARTO_SV_SECTION.out.quarto_section)
        .mix(QUARTO_CNLOH_SECTION.out.quarto_section)
        .mix(QUARTO_FOCAL_SECTION.out.quarto_section)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    emit:
    sections = ch_sections
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def CreateFigureCNVInput(meta, file, type) {
    // Transform chrom into two new variables
    def caption = "${type} Plot for ${meta.id}"
    def section = "CNV"
    def process = "${type}-call-${meta.id}"

    // Return a new tuple with the additional variables
    return [meta, file, caption, section, process ]
}

def CreateSVInput(meta, file, type) {
    def new_meta = modifyMetaId(meta, 'replace', '_sv', '', '')
    def sv_type = file.name.replace("${meta.id}_", "")
    def sv_type_noext = sv_type.replace(".png", "")

    // Transform chrom into two new variables
    def caption = "Figeno Plot for ${new_meta.id} in ${sv_type_noext}"
    def section = "${type}"
    def process = "figeno-${type}-${new_meta.id}-${sv_type_noext}"

    // Return a new tuple with the additional variables
    return [new_meta, file, caption, section, process ]
}
