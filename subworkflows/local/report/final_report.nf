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
    ch_id
    ch_sections
    ch_versions
    ch_mode

    main:


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOFTWARE VERSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

    QUARTO_SECTION(
        ch_section_inputs,
        "Software Versions"
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
        .combine(ch_mode)
        .map { meta, mode ->
            def title = "Oncoseq Pipeline ${mode} Report for ${meta.id}"
            tuple(meta, title) }

    ch_subtitle = ch_id
        .combine(ch_mode)
        .map { meta, mode ->
            def subtitles = "Outputs for the ${mode} branch of the Oncoseq pipeline"
            tuple(meta, subtitles)}

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
