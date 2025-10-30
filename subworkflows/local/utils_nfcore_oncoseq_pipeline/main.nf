//
// Subworkflow with functionality specific to the nf-core/oncoseq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input_sheet             //  string: Path to input samplesheet
    ubam_samplesheet  // string: Path to ubam samplesheet
    demux_samplesheet // string: Path to demux samplesheet
    adaptive_samplesheet // string: Path to adaptive samplesheet ( not null if different for samples)
    input_bed                 // string: path to input bed file (not null if the same for all sample)
    input_padding               // Padding around ROI (parameter padding)
    list_low_fidelity                // List of low fidelity genes to discard for coverage calculations (parameter low_fidelity)

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    // channels from parameters:

    ch_ubam         = channel.fromPath(params.ubam)
    ch_bed          = channel.fromPath(params.bed)
    ch_low_fidelity = channel.fromPath(params.low_fidelity)
    ch_padding      = channel.of(params.padding)
    ch_ref_set      = channel.fromPath(params.ref)
                        .toSortedList()
                        .branch {
                            empty: (it[0].name == "NOFILE")
                                return [[], []]
                            other:true
                                return it
                        }
    ch_ref          = ch_ref_set.empty
                        .mix(ch_ref_set.other)
    ch_ref_id       = channel.of(params.ref_id)

    // Initialize empty channels to fall back on:

    ch_adaptive = channel.empty()
    ch_cfdna    = channel.empty()
    ch_tumor    = channel.empty()

    // Main input channel:

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .combine(ch_ubam)
        .branch {
            meta, input, ubam, _ref, _ref_path ,kit, barcode, id, _tumor_type, _bed, _padding, _low_fidelity, _purity, _filter, ubam_ph ->
            reg: (!ubam)
                return [ meta, file(input), ubam_ph ]
            resume: (ubam)
                return [ meta, file(input), file(ubam) ]
            demux: (barcode && id && kit)
                return [ meta, id, barcode, kit ]
        }
        .set { ch_samplesheet_branched }

    ch_samplesheet_branched.reg
        .mix(ch_samplesheet_branched.resume)
        .set { ch_in_samplesheet }

    ch_samplesheet_branched.demux
        .set { ch_demux }

    if (params.cfdna) {

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .branch {
                meta, _input, _ubam, _ref, _ref_path ,kit, barcode, id, _tumor_type, _bed, _padding, _low_fidelity, purity, filter ->
                simplex: (!kit && !barcode && !id)
                    return [ meta, purity, filter ]
                multiplex: (kit && barcode && id)
                    return tuple(id:id, purity, filter)
            }
            .set { ch_cfdna_branched }

        ch_cfdna_branched.simplex
            .mix(ch_cfdna_branched.multiplex)
            .set { ch_cfdna }

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .branch {
                meta, _input, _ubam, _ref, _ref_path ,_kit, _barcode, id, tumor_type, _bed, _padding, _low_fidelity, _purity, f_ilter ->
                simplex_default: (!id && !tumor_type)
                    return tuple(meta, "leukemia")
                simplex: (!id && tumor_type)
                    return [ meta, tumor_type ]
                multiplex_default: (id && !tumor_type)
                    return tuple(id:id, "leukemia")
                multiplex: (id && tumor_type)
                    return tuple(id:id, tumor_type)
            }
            .set { ch_tumor_branched }

        ch_tumor_branched.simplex_default
            .mix(ch_tumor_branched.simplex)
            .mix(ch_tumor_branched.multiplex_default)
            .mix(ch_tumor_branched.multiplex)
            .set { ch_tumor }


    } else if (params.adaptive) {

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .combine(ch_bed)
            .combine(ch_padding)
            .combine(ch_low_fidelity)
            .branch {
                meta, _input, _ubam, _ref, _ref_path ,kit, barcode, id, _tumor_type, bed, padding, lf, _purity, _filter, bed_c, padding_c, lf_c ->
                common: (!bed && !padding && !lf && !id)
                    return [ meta, bed_c, padding_c, lf_c ]
                common_demux: (!bed && !padding && !lf && kit && barcode && id)
                    return tuple(id:id, bed_c, padding_c, lf_c )
                bed: (bed && !paddng && !lf && !id)
                    return [ meta, file(bed), padding_c, lf_c  ]
                bed_demux: (bed && !padding && !lf && kit && barcode && id)
                    return tuple(id:id, file(bed), padding_c, lf_c )
                padding: (bed && padding && !lf && !id)
                    return [ meta, file(bed), padding, lf_c  ]
                padding_demux: (bed && padding && !lf && kit && barcode && id)
                    return tuple(id:id, file(bed), padding, lf_c)
                bed_diff: (bed && padding && lf && !id)
                    return [ meta, file(bed), padding, file(lf) ]
                bed_diff_demux: (bed && padding && lf && kit && barcode && id)
                    return tuple(id:id, file(bed), padding, file(lf))
            }
            .set { ch_adaptive_branched }

        ch_adaptive_branched.common
            .mix(ch_adaptive_branched.common_demux)
            .mix(ch_adaptive_branched.bed)
            .mix(ch_adaptive_branched.bed_demux)
            .mix(ch_adaptive_branched.padding)
            .mix(ch_adaptive_branched.padding_demux)
            .mix(ch_adaptive_branched.bed_diff)
            .mix(ch_adaptive_branched.bed_diff_demux)
            .set { ch_adaptive }


        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .branch {
                meta, _input, _ubam, _ref, _ref_path ,_kit, _barcode, id, tumor_type, _bed, _padding, _low_fidelity, _purity, _filter ->
                simplex_default: (!id && !tumor_type)
                    return tuple(meta, "leukemia")
                simplex: (!id && tumor_type)
                    return [ meta, tumor_type ]
                multiplex_default: (id && !tumor_type)
                    return tuple(id:id, "leukemia")
                multiplex: (id && tumor_type)
                    return tuple(id:id, tumor_type)
            }
            .set { ch_tumor_branched }

        ch_tumor_branched.simplex_default
            .mix(ch_tumor_branched.simplex)
            .mix(ch_tumor_branched.multiplex_default)
            .mix(ch_tumor_branched.multiplex)
            .set { ch_tumor }

    }

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .combine(ch_ref)
        .combine(ch_ref_id)
        .branch {
            meta, _input, _ubam, ref, ref_path ,kit, _barcode, id, _tumor_type, _bed, _padding, _low_fidelity, _purity, _filter, ref_c, ref_index, ref_id_c ->
            simplex_common: (!ref && !ref_path && !kit)
                return [ meta, ref_id_c, ref_c, ref_index ]
            simplex_diff: (!kit && ref_path && ref)
                return [ meta, ref, file(ref_path)]
            multi_common: (kit && !ref_path && !ref)
                return tuple(id:id, ref_id_c, ref_c, ref_index)
            multi_diff: (kit && ref_path && ref)
                return tuple(id:id, ref, file(ref_path))
        }
        .set { ch_ref_branched }

    ch_ref_branched.simplex_diff
        .mix(ch_ref_branched.multi_diff)
        .set { ch_ref_diff }

    ch_ref_diff
        .map { meta, ref, ref_path ->
                def sorted_ref_path = ref_path.flatten().sort { path ->
                def filename = file(path).name
                // Give priority to reference files over index files
                if (filename.endsWith('.fai')) {
                    return filename.replaceAll('\\.fai$', '') + '_index'
                } else {
                    return filename
                }
            }
            tuple(meta, ref, sorted_ref_path).flatten()}
        .mix(ch_ref_branched.simplex_common)
        .mix(ch_ref_branched.multi_common)
        .set { ch_ref }


    emit:
    bed_sheet   = ch_adaptive
    tumor_type  = ch_tumor
    demux_sheet = ch_demux
    samplesheet = ch_in_samplesheet
    ref_ch      = ch_ref
    cfdna_ch    = ch_cfdna
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Function to modify metadata id field in a flexible way
// This is a generalized function that can handle various metadata transformations
//
def modifyMetaId(Map meta, String operation, String search_string = '', String replace_string = '', def suffix = '') {
    // Clone metadata and normalize to strings
    def new_meta = meta.collectEntries { k, v -> [k, v?.toString()] }

    switch(operation) {
        case 'remove_suffix':
            if (new_meta.id && suffix) {
                new_meta.id = new_meta.id.replaceFirst(suffix.toString(), '')
            }
            break

        case 'add_suffix':
            if (new_meta.id && suffix) {
                new_meta.id = new_meta.id + suffix
            }
            break

        case 'replace':
            if (new_meta.id && search_string) {
                new_meta.id = new_meta.id.replace(search_string, replace_string ?: '')
            }
            break

        case 'prefix':
            if (new_meta.id && suffix) {
                new_meta.id = suffix + new_meta.id
            }
            break
    }

    return new_meta
}

def getMinQC(model) {
    model_string = model.toString()
        if (model_string.contains('sup')) {
            return channel.of(10)
        } else if (model_string.contains('hac')) {
            return channel.of(9)
        } else {
            return channel.of(8)
        }
    }
