/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_TEXT } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE } from '../../../modules/local/quarto/main.nf'
include { QUARTO_REPORT } from '../../../modules/local/quarto/main.nf'
include { softwareVersionsToYAML } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'

workflow MIDNIGHT_REPORT {

    take:
    ch_params
    ch_cfdna
    ch_sections
    ch_versions
    ch_title

    main:


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOFTWARE VERSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_id = ch_params
        .map { meta, refid, ref_fasta, ref_index, bed_file, padding, lowfid ->
            meta }
        .unique()

    // Extract all versions into a single channel of values
    versions = softwareVersionsToYAML(ch_versions)
    // Collapse the channel of versions into a single value
    versions = versions.collect().map { it.join('\n\n') }
    versions = ch_id.combine(versions)

    // Give it an ID of versions
    versions = versions
        .map {
            versions_out ->
            def section = "Versions"
            def process = "versions"

            [versions_out[0], versions_out[1]] + [section, process]
            }

    QUARTO_TEXT(
        versions
        )

    ch_section_inputs = QUARTO_TEXT.out.quarto_text

    ch_section_inputs = ch_section_inputs
        .map { meta, section, section_inputs ->
            [meta, section, section_inputs, 'Software Versions']
        }

    QUARTO_SECTION(
        ch_section_inputs
    )


    // // Add the versions to the channel of sections for every report

    ch_sections = ch_sections.mix(QUARTO_SECTION.out.quarto_section)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INDIVIDUAL REPORTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_report_sections = ch_sections
        .groupTuple()
        .map { id, section, filePaths, reports ->
            [id, section, filePaths, reports]
        }

    ch_title = ch_id
        .join(ch_title)
        .view()

    if (params.cfdna) {

        ch_subtitle = ch_params
            .join(ch_cfdna)
            .map { meta, refid, _ref_fasta, _ref_index, bed_file, _padding, _lowfid, _purity, filter ->
                def bed = bed_file.name
                def skip_modes = params.skip_mapping ? "and bam files as input, skipping mapping." :
                    (params.skip_basecalling ? "and fastq files as input, skipping basecalling." :
                        "using pod5 files as input, basecalling model ${params.basecall_model}")
                def modif = params.m_bases ? " and modified basecalling model ${params.m_bases}" : ""
                def filt = filter == "yes" ? "Reads shorter than ${params.max_length} bp were filtered out before mapping." : ""
                tuple(meta, "Ouputs summary for oncoseq workflow, using reference genome ${refid}, targeted panel bed file ${bed}, ${skip_modes}${modif}. ${filt}")}

    } else {
        ch_subtitle = ch_params
            .map { meta, refid, _ref_fasta, _ref_index, bed_file, _padding, _lowfid ->
                def bed = bed_file.name
                def skip_modes = params.skip_mapping ? "and bam files as input, skipping mapping." :
                    (params.skip_basecalling ? "and fastq files as input, skipping basecalling." :
                        "using pod5 files as input, basecalling model ${params.basecall_model}")
                def modif = params.m_bases ? " and modified basecalling model ${params.m_bases}" : ""
                tuple(meta, "Ouputs summary for oncoseq workflow, using reference genome ${refid}, targeted panel bed file ${bed}, ${skip_modes}${modif}")}
    }

    ch_template = channel.fromPath(params.report_template)

    ch_report_in = ch_report_sections
        .join(ch_title)
        .join(ch_subtitle)
        .combine(ch_template)

    QUARTO_REPORT(
        ch_report_in
        )

    ch_report = QUARTO_REPORT.out.report


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    emit:
    ch_report
}
