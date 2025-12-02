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

    ch_ubam         = channel.fromPath(params.ubam, checkIfExists: true)
    ch_bed          = channel.fromPath(params.bed, checkIfExists: true)
    ch_low_fidelity = channel.fromPath(params.low_fidelity, checkIfExists: true)
    ch_padding      = channel.of(params.padding)
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
            meta, project, input, ubam, _ref, _ref_path ,kit, barcode, _tumor_type, _bed, _padding, _low_fidelity, _purity, _filter, ubam_ph ->
            reg: (!ubam && !kit && !barcode)
                return [ meta, file(input), ubam_ph ]
            resume: (ubam && !kit && !barcode)
                return [ meta, file(input), file(ubam) ]
            demux_reg: (project && kit && barcode && !ubam)
                return tuple(id:project, file(input), ubam_ph )
            demux_resume: (project && kit && barcode)
                return tuple(id:project, file(input), file(ubam))

        }
        .set { ch_samplesheet_branched }

    ch_samplesheet_branched.reg
        .mix(ch_samplesheet_branched.resume)
        .mix(ch_samplesheet_branched.demux_reg)
        .mix(ch_samplesheet_branched.demux_resume)
        .unique()
        .set { ch_in_samplesheet }

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _tumor_type, _bed, _padding, _low_fidelity, _purity, _filter ->
            if (project && kit && barcode) {
                tuple(id:project, meta.id, barcode, kit)
            } else {
                tuple()
            }
        }
        .set { ch_demux }

    if (params.cfdna) {

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _tumor_type, _bed, _padding, _low_fidelity, purity, filter ->
                if(!purity){
                    throw new IllegalArgumentException("Please provide sample purity estimation for sample: ${meta.id ?: meta}")
                }
                if(!filter){
                    log.warn("Warning: No filter value provided for sample ${meta.id ?: meta} — defaulting to true (samples will be filtered for fragment length <= 700bp)")
                    filter = true
                }

                return tuple(meta, purity, filter)
            }
            .set { ch_cfdna }

        // Make it possible to merge reports together for cfdna
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .branch { meta, project, _input, _ubam, _ref, _ref_path ,_kit, _barcode, tumor_type, _bed, _padding, _low_fidelity, _purity, _filter ->
                individual: (!project)
                    return tuple(meta, meta.id)
                common: (project)
                    return tuple(meta, project) }
            .set { ch_id_branched }

        ch_id_branched.individual
            .mix(ch_id_branched.common)
            .set { ch_id }

        // Make channel with empty bed file for clair3 calling
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .combine(ch_bed)
            .combine(ch_padding)
            .combine(ch_low_fidelity)
            .map {
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _tumor_type, bed, padding, lf, _purity, _filter, bed_c, padding_c, lf_c ->
                tuple(meta, bed_c, padding_c, lf_c)
            }
            .set { ch_adaptive }

    } else if (params.adaptive) {

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .combine(ch_bed)
            .combine(ch_padding)
            .combine(ch_low_fidelity)
            .branch {
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _tumor_type, bed, padding, lf, _purity, _filter, bed_c, padding_c, lf_c ->
                common: (!bed && !padding && !lf)
                    return [ meta, bed_c, padding_c, lf_c ]
                bed: (bed && !padding && !lf)
                    return [ meta, file(bed), padding_c, lf_c  ]
                padding: (bed && padding && !lf)
                    return [ meta, file(bed), padding, lf_c  ]
                bed_diff: (bed && padding && lf)
                    return [ meta, file(bed), padding, file(lf) ]
            }
            .set { ch_adaptive_branched }

        ch_adaptive_branched.common
            .mix(ch_adaptive_branched.bed)
            .mix(ch_adaptive_branched.padding)
            .mix(ch_adaptive_branched.bed_diff)
            .set { ch_adaptive }

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { meta, project, _input, _ubam, _ref, _ref_path ,_kit, _barcode, tumor_type, _bed, _padding, _low_fidelity, _purity, _filter ->
                tuple(meta, meta.id)
            }
            .set { ch_id }

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ERROR MESSAGES
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .combine(ch_bed)
            .combine(ch_padding)
            .combine(ch_low_fidelity)
            .map { meta, project, input, ubam, ref, ref_path ,kit, barcode, tumor_type, bed, padding, low_fidelity, purity, filter, bed_c, padding_c, lf_c ->
                if (params.adaptive && !padding && padding_c == 20000){
                    log.warn("Using default padding of 20kb padding around ROI")
                }
                if (!bed && bed_c.name == "NO_BED") {
                    throw new IllegalArgumentException("Missing bed file for adaptive: please provide a bed file containing ROIs")
                }
            }

    } else {

        // Make channel with empty bed file for clair3 calling
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .combine(ch_bed)
            .combine(ch_padding)
            .combine(ch_low_fidelity)
            .map {
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _tumor_type, bed, padding, lf, _purity, _filter, bed_c, padding_c, lf_c ->
                tuple(meta, bed_c, padding_c, lf_c)
            }
            .set { ch_adaptive }

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { meta, project, _input, _ubam, _ref, _ref_path ,_kit, _barcode, tumor_type, _bed, _padding, _low_fidelity, _purity, _filter ->
                tuple(meta, meta.id)
            }
            .set { ch_id }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Tumor Type
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .branch {
            meta, project, _input, _ubam, _ref, _ref_path ,_kit, _barcode, tumor_type, _bed, _padding, _low_fidelity, _purity, _filter ->
            leukemia: (!tumor_type)
                return tuple(meta, "leukemia")
            other: (tumor_type)
                return tuple(meta, tumor_type)
        }
        .set { ch_tumor_branched }

   ch_tumor_branched.leukemia
        .mix(ch_tumor_branched.other)
        .set { ch_tumor }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    References files
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.ref == null) {
        ch_ref = channel.of(params.ref)
        ch_ref_index = channel.of(params.ref)
        ch_ref = ch_ref.combine(ch_ref_index)
    } else {
        ch_ref = channel.fromPath(params.ref, checkIfExists:true)
            .map { ref ->
                def index = file("${ref}.fai")
                if( !index.exists() )
                    throw new IllegalArgumentException("Missing index: ${index}")
                tuple(file(ref), index)
            }
    }

    // Throw error if no reference is provided

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .combine(ch_ref)
        .combine(ch_ref_id)
        .map {
            meta, project, _input, _ubam, ref, ref_path ,kit, _barcode, _tumor_type, _bed, _padding, _low_fidelity, _purity, _filter, ref_c, ref_index, ref_id_c ->
            if(!ref_path && ref_c == null) {
                throw new IllegalArgumentException("No reference file provided, please provide a reference throuh --ref or samplesheet")
            }
        }

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .combine(ch_ref)
        .combine(ch_ref_id)
        .branch {
            meta, project, _input, _ubam, ref, ref_path ,kit, _barcode, _tumor_type, _bed, _padding, _low_fidelity, _purity, _filter, ref_c, ref_index, ref_id_c ->
            common: (!ref && !ref_path)
                return [ meta, ref_id_c, ref_c, ref_index ]
            diff: (ref_path && ref)
                def fai = file("${ref_path}.fai")
                if(!fai.exists())
                    throw new IllegalArgumentException("Missing index for reference: ${ref} (expected: ${fai})")

                return [ meta, ref, file(ref_path), fai ]
        }
        .set { ch_ref_branched }

    ch_ref_branched.diff
        .mix(ch_ref_branched.common)
        .set { ch_ref }

    // Error message for demultiplexing
    if (params.demux) {
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _tumor_type, _bed, _padding, _low_fidelity, _purity, _filter ->
                if (project && kit && !barcode ) {
                    throw new IllegalArgumentException("Please provide barcode for sample ${meta.id ?: meta}")
                }
                if (project && !kit && barcode) {
                    throw new IllegalArgumentException("Please provide multplexing kit used for sample ${meta.id ?: meta}")
                }
                if (project && !kit && !barcode) {
                    throw new IllegalArgumentException("Please provide multplexing kit used and barcode for sample ${meta.id ?: meta}")
                }
            }
    } else {
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _tumor_type, _bed, _padding, _low_fidelity, _purity, _filter ->
                if (kit || barcode) {
                    throw new IllegalArgumentException("Please use --demux option for demultiplexing samples")
                }
            }
    }

    emit:
    bed_sheet   = ch_adaptive
    tumor_type  = ch_tumor
    demux_sheet = ch_demux
    samplesheet = ch_in_samplesheet
    ref_ch      = ch_ref
    cfdna_ch    = ch_cfdna
    id_ch       = ch_id
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
            return 10
        } else if (model_string.contains('hac')) {
            return 9
        } else {
            return 8
        }
}
