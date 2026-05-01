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
include { SAMTOOLS_FAIDX            } from '../../../modules/local/samtools/main.nf'
include { BGZIP_RECOMPRESS          } from '../../../modules/local/bcftools/main.nf'

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

    // Initialize empty channels to fall back on:

    ch_adaptive = channel.empty()
    ch_cfdna    = channel.empty()

    // Main input channel:

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .combine(ch_ubam)
        .branch {
            meta, project, input, ubam, _ref, _ref_path ,kit, barcode, _bed, _padding, _low_fidelity, _purity, _filter, ubam_ph ->
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
            meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _bed, _padding, _low_fidelity, _purity, _filter ->
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
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _bed, _padding, _low_fidelity, purity, filter ->
                if(!purity){
                    throw new IllegalArgumentException("Please provide sample purity estimation for sample: ${meta.id ?: meta}")
                }
                if(!filter){
                    log.warn("Warning: No filter value provided for sample ${meta.id ?: meta} — defaulting to yes (samples will be filtered for fragment length <= 700bp)")
                    filter = "yes"
                }

                return tuple(meta, purity, filter)
            }
            .set { ch_cfdna }

        // Make channel with empty bed file for clair3 calling
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .combine(ch_bed)
            .map {
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, bed, padding, lf, _purity, _filter, bed_c ->
                tuple(meta, bed_c)
            }
            .set { ch_adaptive }

    } else if (params.adaptive) {

        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .combine(ch_bed)
            .combine(ch_padding)
            .combine(ch_low_fidelity)
            .branch {
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, bed, padding, lf, _purity, _filter, bed_c, padding_c, lf_c ->
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
            .map { meta, project, input, ubam, ref, ref_path ,kit, barcode, bed, padding, low_fidelity, purity, filter, bed_c, padding_c, lf_c ->
                if (params.adaptive && !padding && padding_c == 20000){
                    log.warn("Using default padding of 20kb padding around ROI")
                }
                if (!bed && bed_c.name.contains("CRA-CHUSJ")) {
                    log.warn("Using default bed file ${bed_c.name}")
                }
            }

    } else {

        // Make channel with empty bed file for clair3 calling
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .combine(ch_bed)
            .map {
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, bed, padding, lf, _purity, _filter, bed_c->
                if (bed) {
                    tuple(meta,bed)
                } else {
                tuple(meta, bed_c)
                }
            }
            .set { ch_adaptive }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    References files
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def genomeMap = [
        'hg38'   : 'hg38',
        'GRCh38' : 'hg38',
        'hg19'   : 'hg19',
        'GRCh37' : 'hg19',
        'hs1'    : 'hs1',
        'CHM13'  : 'hs1'
    ]

    ch_genome    = params.genome ? channel.of(params.genome) : channel.empty()
    ch_ref_cache = params.ref_cache
        ? channel.fromPath(params.ref_cache, checkIfExists: true)
        : channel.empty()

    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .combine(ch_ref_cache.ifEmpty([null]))
        .combine(ch_genome.ifEmpty([null]))
        .map { meta, project, _input, _ubam, sample_genome, ref_path, kit, _barcode, _bed, _padding, _low_fidelity, _purity, _filter, ref_c, genome_c ->

            def faExts  = ['.fa', '.fasta', '.fa.gz', '.fasta.gz']
            def faiExts = ['.fa.fai', '.fasta.fai', '.fa.gz.fai', '.fasta.gz.fai']

            // Resolve genome: sample-level takes priority over global param
            def rawGenome = sample_genome ?: genome_c
            if (!rawGenome) error "ERROR: No genome provided. Set params.genome or provide genome in samplesheet."
            def genome  = genomeMap[rawGenome] ?: { error "ERROR: Unrecognized genome '${rawGenome}'. Valid values: ${genomeMap.keySet().join(', ')}" }()
            def aliases = genomeMap.findAll { k, v -> v == genome }.keySet()

            // Case-insensitive FASTA/FAI finders
            def findFa  = { files, aliasList -> files.find { f -> aliasList.any { a -> f.name.toLowerCase().contains(a.toLowerCase()) } && faExts.any  { ext -> f.name.toLowerCase().endsWith(ext) } } }
            def findFai = { files, aliasList -> files.find { f -> aliasList.any { a -> f.name.toLowerCase().contains(a.toLowerCase()) } && faiExts.any { ext -> f.name.toLowerCase().endsWith(ext) } } }

            def faFinal  = null
            def faiFinal = null

            // --- Case 1: sample-level ref_path provided ---
            if (ref_path) {
                def refFile = ref_path instanceof Path ? ref_path : file(ref_path)

                // --- Case 1a: ref_path is a direct FASTA file ---
                if (refFile.isFile()) {
                    if (faExts.any { ext -> refFile.name.endsWith(ext) }) {
                        faFinal  = refFile
                        def faiFile = file("${refFile.parent}/${refFile.name}.fai")
                        faiFinal = faiFile.exists() ? faiFile : null
                        if (!faiFinal) log.warn("No FAI found next to '${refFile}' — will index")
                    } else {
                        error "ERROR: ref_path '${refFile}' does not look like a FASTA file. Expected extensions: ${faExts}"
                    }

                // --- Case 1b: ref_path is a directory ---
                } else if (refFile.isDirectory()) {
                    def allFiles = refFile.listFiles() as List
                    def faFile   = findFa(allFiles, aliases)
                    def faiFile  = findFai(allFiles, aliases)

                    faFinal  = faFile  ?: file("ftp://hgdownload.cse.ucsc.edu/goldenPath/${genome}/bigZips/${genome}.fa.gz")
                    faiFinal = faiFile ?: null

                    if (!faFile)  log.warn("No FASTA found in ref_path '${refFile}' for '${genome}' — falling back to UCSC FTP")
                    if (!faiFile) log.warn("No FAI found in ref_path '${refFile}' for '${genome}' — will index")

                } else {
                    error "ERROR: ref_path '${refFile}' does not exist"
                }

            // --- Case 2: global ref_cache provided ---
            } else if (ref_c) {
                def refCacheFile = ref_c instanceof Path ? ref_c : file(ref_c)

                // --- Case 2a: ref_cache is a direct FASTA file ---
                if (refCacheFile.isFile()) {
                    if (faExts.any { ext -> refCacheFile.name.endsWith(ext) }) {
                        faFinal  = refCacheFile
                        def faiFile = file("${refCacheFile.parent}/${refCacheFile.name}.fai")
                        faiFinal = faiFile.exists() ? faiFile : null
                        if (!faiFinal) log.warn("No FAI found next to '${refCacheFile}' — will index")
                    } else {
                        error "ERROR: ref_cache '${refCacheFile}' does not look like a FASTA file. Expected extensions: ${faExts}"
                    }

                // --- Case 2b: ref_cache is a directory ---
                } else if (refCacheFile.isDirectory()) {
                    def allFiles = refCacheFile.listFiles() as List
                    def faFile   = findFa(allFiles, aliases)
                    def faiFile  = findFai(allFiles, aliases)

                    faFinal  = faFile  ?: file("ftp://hgdownload.cse.ucsc.edu/goldenPath/${genome}/bigZips/${genome}.fa.gz")
                    faiFinal = faiFile ?: null

                    if (!faFile)  log.warn("No FASTA found in ref_cache '${refCacheFile}' for '${genome}' — falling back to UCSC FTP")
                    if (!faiFile) log.warn("No FAI found in ref_cache '${refCacheFile}' for '${genome}' — will index")

                } else {
                    error "ERROR: ref_cache '${refCacheFile}' does not exist"
                }

            // --- Case 3: no ref_cache, no ref_path — FTP on demand ---
            } else {
                log.warn("No ref_cache or ref_path provided — staging FTP genome '${genome}' from UCSC")
                faFinal  = file("ftp://hgdownload.cse.ucsc.edu/goldenPath/${genome}/bigZips/${genome}.fa.gz")
                faiFinal = null
            }

            tuple(meta, genome, faFinal, faiFinal)
        }
        .set { ch_ref }

    ch_ref
        .branch { meta, ref_id_c, faFinal, faiFinal ->        // <-- missing -> fixed
            indexed  : faiFinal != null
            to_index : faiFinal == null
        }
        .set { ch_ref_branched }                              // <-- branch result must be .set{}

    // If genome is compressed, always recompress with bgzip otherwise samtools faidx will fail

    ch_ref_branched.to_index
        .branch { meta, ref_id_c, faFinal, faiFinal ->
            gzip: faFinal.extension == "gz"
            unzip: (faFinal.extension == "fa" || faFinal.extension == "fasta")
        }
        .set { ch_ref_to_index_branched }

    BGZIP_RECOMPRESS(ch_ref_to_index_branched.gzip
        .map { meta, ref_id_c, faFinal, _faiFinal -> tuple(meta, faFinal) }
        .groupTuple(by: 1)
        )

    // Reunite bgzipped FTP files with local files, then index all
    ch_to_index = BGZIP_RECOMPRESS.out.file
        .mix(ch_ref_to_index_branched.unzip.map { meta, ref_id_c, faFinal, _faiFinal -> tuple(meta, faFinal) })
        .transpose()

    SAMTOOLS_FAIDX(ch_to_index
        .groupTuple(by: 1)
    )

    // Rejoin with ref IDs
    ch_ref_ids = ch_ref
        .map { meta, ref_id_c, faFinal, faiFinal ->
        tuple(meta, ref_id_c) }

    ch_ref_indexed = ch_ref_ids
        .join(SAMTOOLS_FAIDX.out.fasta_index.transpose())
        .mix(ch_ref_branched.indexed)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Error messages
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.demux) {
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _bed, _padding, _low_fidelity, _purity, _filter ->
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
                meta, project, _input, _ubam, _ref, _ref_path ,kit, barcode, _bed, _padding, _low_fidelity, _purity, _filter ->
                if (kit || barcode) {
                    throw new IllegalArgumentException("Please use --demux option for demultiplexing samples")
                }
            }
    }

    emit:
    bed_sheet   = ch_adaptive
    demux_sheet = ch_demux
    samplesheet = ch_in_samplesheet
    ref_ch      = ch_ref_indexed
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
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Joudnal, DOI</li>" : "",
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

def selectLatestModel(model,model_path) {
    def model_name = model.toString()
    def requested_version = null

    def version_match = (model_name =~ /@v[0-9]+(\.[0-9]+)*/)
    if (version_match.find()) {
        requested_version = version_match.group(0).substring(2)
        model_name = model_name.split('@v')[0]
    }

    def version_pattern = requested_version
        ? java.util.regex.Pattern.quote(requested_version)
        : '[0-9]+(\\.[0-9]+)*'

    // 1) Find clean (unmodified) models
    def clean_models = model_path.listFiles().findAll { f ->
        f.name ==~ /.*${model_name}@v${version_pattern}$/
    }

    def clean_model_dna = clean_models.findAll { f ->
        f.name.contains('dna')
    }

    if (clean_model_dna.size() == 0) {
        def version_note = requested_version ? " (requested version ${requested_version})" : ""
        log.warn(
            "No model file found in the provided ${model_path} for model ${model_name}${version_note}: model will be downloaded"
        )
        return null
    } else {
        if (requested_version) {
            return clean_model_dna.first()
        }

        // 2) Select latest clean model
        return clean_model_dna
            .sort { a, b ->
                def va = a.name.split('@v')[-1].tokenize('.').collect { it as int }
                def vb = b.name.split('@v')[-1].tokenize('.').collect { it as int }

                int len = Math.max(va.size(), vb.size())
                for (int i = 0; i < len; i++) {
                    int ai = i < va.size() ? va[i] : 0
                    int bi = i < vb.size() ? vb[i] : 0
                    if (ai != bi) return ai <=> bi
                }
                return 0
            }
            .last()
    }
}

def selectLatestModif(model_path, model, modif) {
    def base_model = model.name
    def modif_full = modif.toString()
    def modif_name = modif_full
    def requested_version = null

    if (modif_full.startsWith(base_model + '_')) {
        modif_name = modif_full.substring(base_model.length() + 1)
    }

    def lastVIndex = modif_name.lastIndexOf('@v')
    if (lastVIndex != -1) {
        requested_version = modif_name.substring(lastVIndex + 2)
        modif_name = modif_name.substring(0, lastVIndex)
    }

    def version_pattern = requested_version
        ? java.util.regex.Pattern.quote(requested_version)
        : '[0-9]+(\\.[0-9]+)*'

    def clean_model_mod = model_path.listFiles().findAll { f ->
        f.name ==~ /${base_model}_${modif_name}@v${version_pattern}$/
    }

    if (!clean_model_mod) {
        def version_note = requested_version ? " (requested version ${requested_version})" : ""
        log.warn(
            "No modified model found for ${base_model} with ${modif_name}${version_note} in ${model_path}: model will be downloaded"
        )
        return null
    } else {
        if (requested_version) {
            return clean_model_mod.first()
        }
        return clean_model_mod
            .sort { a, b ->
                def va = a.name.split('@v')[-1].tokenize('.')*.toInteger().join('').toInteger()
                def vb = b.name.split('@v')[-1].tokenize('.')*.toInteger().join('').toInteger()
                va <=> vb
            }
            .last()
    }
}

def normalizeVersion(versionString) {
    versionString
        .tokenize('.')
        .collect { String.format('%03d', it as int) }
        .join()
}

def selectModelDownload(chModelsList, modelParam, subfield = null, chParentModel = null) {

    def chParsed = chModelsList
        .splitJson(path: 'dna_r10.4.1_e8.2_400bps_5khz.simplex_models')

    /*
     * CASE 1: Base model selection
     */
    if (!subfield) {
        return chParsed
            .filter { it.key.contains(modelParam) }
            .map { entry ->
                def v = (entry.key =~ /@v([0-9.]+)/)[0][1]
                tuple(entry.key, normalizeVersion(v))
            }
            .ifEmpty { error "No models found for: ${modelParam}" }
            .toSortedList { a, b -> a[1] <=> b[1] }
            .map { list ->
                def exact = list.find { it[0] == modelParam }
                exact ? exact[0] : list[-1][0]
            }
    }

    /*
     * CASE 2: Sub-model selection (modified, polish, etc.)
     */
    return chParsed
        .combine(chParentModel)
        .filter { entry, parent ->
            entry.key == parent
        }
        .map { entry, parent ->
            (entry.value[subfield] ?: [:]).entrySet()
        }
        .flatMap { it }
        .filter { entry ->
            entry.key.contains(modelParam) || entry.key == modelParam
        }
        .map { entry ->
            def v = (entry.key =~ /@v([0-9.]+)/)[0][1]
            tuple(entry.key, normalizeVersion(v))
        }
        .ifEmpty { error "No ${subfield} models found for: ${modelParam}" }
        .toSortedList { a, b -> a[1] <=> b[1] }
        .map { list ->
            def exact = list.find { it[0] == modelParam }
            exact ? exact[0] : list[-1][0]
        }
}
