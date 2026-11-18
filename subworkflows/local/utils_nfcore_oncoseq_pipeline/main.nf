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
include { SAMTOOLS_FAIDX            } from '../../../modules/nf-core/samtools/faidx/main.nf'
include { PIGZ_UNCOMPRESS as PIGZ_UNGZ_BED } from '../../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_UNCOMPRESS as PIGZ_UNGZ_FA  } from '../../../modules/nf-core/pigz/uncompress/main'
include { BEDTOOLS_SLOP             } from '../../../modules/nf-core/bedtools/slop/main'
include { DORADO_DOWNLOAD_LIST    } from '../../../modules/local/dorado/main.nf'
include { DORADO_DOWNLOAD_MODEL   } from '../../../modules/local/dorado/main.nf'
include { ENSEMBLVEP_DOWNLOAD     } from '../../../modules/nf-core/ensemblvep/download/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input_sheet             //  string: Path to input samplesheet
    basecall_model          //  string: basecalling model
    ref_cache            //  string: path to reference cache directory
    vep_cache            //  string: path to Ensembl VEP cache directory

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
        "${projectDir}/nextflow_schema.json",
        false,
        false,
        false,
        params.before_text ?: "",
        params.after_text ?: "",
        "nextflow run . --input ${input_sheet} --basecall_model ${basecall_model} --ref_cache ${ref_cache} --vep_cache ${vep_cache} --outdir ${outdir}",
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )


    def cols = ['meta', 'project', 'input', 'ubam', 'kit', 'barcode',
                'bed', 'padding', 'lf', 'purity', 'filter', 'tumor_type']

    ch_rows = channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { row -> [cols, row].transpose().collectEntries() }

    // Global defaults
    def ubam_ph     = file(params.ubam, checkIfExists: true)
    def bed_default = file(params.bed, checkIfExists: true)
    def lf_default  = file(params.low_fidelity, checkIfExists: true)

    // --- Main input: [ meta, input, ubam ] ---
    ch_in_samplesheet = ch_rows
        .map { row ->
            def demux = row.project && row.kit && row.barcode
            if (!demux && (row.kit || row.barcode)) {
                error "Sample '${row.meta.id}': kit and barcode require a project (demultiplexing)"
            }
            tuple(
                demux ? [id: row.project] : row.meta,
                file(row.input),
                row.ubam ? file(row.ubam) : ubam_ph
            )
        }
        .unique()

    // --- Demultiplexing info: [ [id: project], sample_id, barcode, kit ] ---
    ch_demux = ch_rows
        .filter { row -> row.project && row.kit && row.barcode }
        .map    { row -> tuple([id: row.project], row.meta.id, row.barcode, row.kit) }

    // --- cfDNA: [ meta, purity, filter ] ---
    ch_cfdna = params.cfdna
        ? ch_rows.map { row ->
            if (!row.purity) {
                error "Please provide sample purity estimation for sample: ${row.meta.id}"
            }
            if (!row.filter) {
                log.warn("No filter value provided for sample ${row.meta.id} — defaulting to yes (fragment length <= 700bp)")
            }
            tuple(row.meta, row.purity, row.filter ?: 'yes')
        }
        : channel.empty()

    // --- Tumor type: [ meta, tumor_type ] ---
    ch_tumor_type = ch_rows
        .map { row -> tuple(row.meta, row.tumor_type ?: 'other') }

    // --- Adaptive: [ meta, bed, padding, lf ], per-sample value or global default ---
    ch_adaptive_branched = ch_rows
        .map { row ->
            if (!row.padding) {
                def padding_kb = params.padding / 1000
                log.warn("Using padding of ${padding_kb}kb around ROI for all samples")
            }
            if (!row.bed && params.bed.contains("v2.0.1-pre-merge-panel-20kb-pad.bed")) {
                log.warn("Using default bed file ${params.bed} for all samples")
            }
            tuple(
                row.meta,
                row.bed ? file(row.bed) : bed_default,
                row.padding instanceof Number ? row.padding : params.padding,
                row.lf  ? file(row.lf)  : lf_default
            )
        }
        .branch { _meta, bed, _padding, _lf ->
            gz   : bed.name.endsWith('.gz') || bed.name.endsWith('.bgz')
            plain: true
        }

    // Decompress each unique gzipped BED once. The source path is kept in meta
    // so the result can be joined back without a hash-based key.
    ch_bed_to_decompress = ch_adaptive_branched.gz
        .map    { _meta, bed, _padding, _lf -> bed }
        .unique()
        .map    { bed -> tuple([id: bed.simpleName, src: bed], bed) }

    PIGZ_UNGZ_BED(
        ch_bed_to_decompress,
    )

    ch_bed_decompressed = PIGZ_UNGZ_BED.out.file
        .map { meta, bed -> tuple(meta.src, bed) }

    // combine(by: 0) because many samples can share one BED
    ch_adaptive = ch_adaptive_branched.plain.mix(
        ch_adaptive_branched.gz
            .map     { meta, bed, padding, lf -> tuple(bed, meta, padding, lf) }
            .combine(ch_bed_decompressed, by: 0)
            .map     { _src, meta, padding, lf, bed -> tuple(meta, bed, padding, lf) }
    )

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
        'CHM13'  : 'hs1',
        't2t'    : 'hs1'
    ]

    ch_genome    = params.genome ? channel.of(params.genome) : channel.empty()
    ch_ref_cache = params.ref_cache
        ? channel.fromPath(params.ref_cache, checkIfExists: true)
        : channel.empty()

    ch_genome
        .combine(ch_ref_cache.ifEmpty([null]))
        .map { genome, ref_dir ->

            def faExts  = ['.fa', '.fasta', '.fa.gz', '.fasta.gz']
            def faiExts = ['.fa.fai', '.fasta.fai', '.fa.gz.fai', '.fasta.gz.fai']

            // Resolve genome
            if (!genome) error "ERROR: No genome provided. Set params.genome or provide genome in samplesheet."
            def genomeRes = genomeMap[genome] ?: { error "ERROR: Unrecognized genome '${genome}'. Valid values: ${genomeMap.keySet().join(', ')}" }()
            def aliases   = genomeMap.findAll { k, v -> v == genomeRes }.keySet()

            def faFinal  = null
            def faiFinal = null

            // --- Case 1: ref_path provided ---
            if ( ref_dir ) {
                def refFile = file(ref_dir)

                // --- Case 1a: ref_path is a direct FASTA file ---
                if (refFile.isFile()) {
                    if (faExts.any { ext -> refFile.name.endsWith(ext) }) {
                        faFinal = refFile
                        def faiFile = file("${refFile.parent}/${refFile.name}.fai")
                        faiFinal = faiFile.exists() ? faiFile : null
                        if (!faiFinal) log.warn("No FAI found next to '${refFile}' — will index")
                    } else {
                        error "ERROR: ref_path '${refFile}' does not look like a FASTA file. Expected extensions: ${faExts}"
                    }

                // --- Case 1b: ref_path is a directory ---
                } else if (refFile.isDirectory()) {
                    def allFiles = refFile.listFiles() as List
                    def faFile   = findFaFile(allFiles, aliases, faExts)
                    def faiFile  = findFaFile(allFiles, aliases, faiExts)
                    def url = genomeRes == "hg38" ? "ftp://hgdownload.cse.ucsc.edu/goldenPath/${genomeRes}/bigZips/p14/${genomeRes}.p14.fa.gz"
                        : "ftp://hgdownload.cse.ucsc.edu/goldenPath/${genomeRes}/bigZips/${genomeRes}.fa.gz"

                    faFinal  = faFile  ?: file(url)
                    faiFinal = faiFile ?: null

                    if (!faFile)  log.warn("No FASTA found in ref_path '${refFile}' for '${genome}' — falling back to UCSC FTP")
                    if (!faiFile) log.warn("No FAI found in ref_path '${refFile}' for '${genome}' — will index")

                } else {
                    error "ERROR: ref_path '${refFile}' does not exist"
                }

            // --- Case 2: no ref_cache, no ref_path — FTP on demand ---
            } else {
                log.warn("No ref_cache or ref_path provided — staging FTP genome '${genome}' from UCSC")
                faFinal  = file("ftp://hgdownload.cse.ucsc.edu/goldenPath/${genomeRes}/bigZips/${genomeRes}.fa.gz")
                faiFinal = null
            }

            tuple([id: genomeRes], faFinal, faiFinal)
        }
        .set { ch_ref }

    ch_ref
        .branch { meta, faFinal, faiFinal ->
            // gzipped FASTA: decompress first (FAI may or may not exist)
            gz : faFinal.name.endsWith('gz')
            // uncompressed FASTA
            fa_toindex: (faFinal.name.endsWith('.fa') || faFinal.name.endsWith('.fasta')) && !faiFinal
                return tuple(meta, faFinal, [])
            fa_indexed: (faFinal.name.endsWith('.fa') || faFinal.name.endsWith('.fasta')) && faiFinal
        }
        .set { ch_ref_branched }

    // --- Decompress gzipped FASTA files ---
    // Input: tuple(meta, infile, infile_tbi, regions)
    ch_to_decompress = ch_ref_branched.gz
        .map { meta, fasta, _fai -> tuple(meta, fasta) }

    PIGZ_UNGZ_FA(
        ch_to_decompress
    )

    // --- Index if fai file is missing ---
    ch_to_index = ch_ref_branched.fa_toindex
        .mix(PIGZ_UNGZ_FA.out.file.map { meta, fasta -> tuple(meta, fasta, []) })

    SAMTOOLS_FAIDX(
        ch_to_index,
        false
    )

    // FASTA channel: tuple(meta, fasta) — always uncompressed
    ch_ref_fasta = ch_to_index
        .mix(ch_ref_branched.fa_indexed)
        .map { meta, fasta, _fai -> tuple(meta, fasta) }

    // FAI channel: tuple(meta, fai)
    ch_ref_fai = SAMTOOLS_FAIDX.out.fai
        .mix(ch_ref_branched.fa_indexed.map { meta, _fasta, fai -> tuple(meta, fai) })

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Create bed file wihtout padding around ROI
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_slop_in = ch_adaptive
        .groupTuple(by: [1,2])
        .map { metas, bed, padding, _lf ->
            def padding_minus = padding > 0 ? -padding : 0
            tuple([id: bed.baseName, padding: padding_minus, samples: metas], bed)
        }
        .filter { meta, _bed -> meta.padding != 0 }

    BEDTOOLS_SLOP(
        ch_slop_in,
        ch_ref_fai.map { _meta, fai -> fai }.first(),
    )

    ch_bednopad = BEDTOOLS_SLOP.out.bed
        .flatMap { meta, bed ->
            meta.samples.collect { sample_meta -> tuple(sample_meta, bed) }
        }
        .mix(
            ch_adaptive
                .filter { _meta, _bed, padding, _lf -> padding <= 0 }
                .map { meta, bed, _padding, _lf -> tuple(meta, bed) }
        )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Basecalling model checkup and download
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.skip_basecalling || params.skip_mapping ) {

        ch_model = channel.empty()
        ch_modif = channel.empty()

    } else {
        if (!params.m_bases) {
            log.warn("Tumor classification will be skipped — parameter --m_bases not set")
        }

        def ref_cache_dir      = params.ref_cache ? file(params.ref_cache) : null
        def ref_cache_resolved = ref_cache_dir ? (ref_cache_dir.isDirectory() ? ref_cache_dir.toString() : ref_cache_dir.getParent()) : null

        ch_model_dir = ref_cache_resolved
            ? channel.fromPath(ref_cache_resolved, checkIfExists: true)
            : channel.fromPath("${launchDir}")

        def model_resolved = ref_cache_resolved
            ? selectLatestModel(params.basecall_model, file("${ref_cache_resolved}/dorado_models"))
            : null
        def modif_resolved = params.m_bases && model_resolved ?
            selectLatestModif(file("${params.ref_cache}/dorado_models"), model_resolved, params.m_bases) :
            null
        ch_modif = channel.fromPath("${projectDir}/assets/NOMOD")

        if (model_resolved) {
            ch_model = channel.fromPath(model_resolved)
        }

        if (modif_resolved) {
            ch_modif = channel.fromPath(modif_resolved)
        }

        // download both
        if (!model_resolved && !modif_resolved && params.m_bases) {

            DORADO_DOWNLOAD_LIST()

            ch_model_to_download_resolved = selectBaseModelDownload(
                DORADO_DOWNLOAD_LIST.out.list,
                params.basecall_model
            )

            ch_model_to_download = ch_model_to_download_resolved
                .combine(ch_model_dir)
                .map { type, model -> tuple("base", type, model)}

            ch_modif_to_download = selectModifiedModelDownload(
                DORADO_DOWNLOAD_LIST.out.list,
                ch_model_to_download_resolved,
                params.m_bases
            )
                .combine(ch_model_dir)
                .map { type, model -> tuple("modif", type, model)}

            ch_in_download = ch_model_to_download
                .mix(ch_modif_to_download)

            DORADO_DOWNLOAD_MODEL(ch_in_download)

            ch_model = DORADO_DOWNLOAD_MODEL.out.model
                .filter { type, _model -> type == "base" }
                .map { _type, model -> model }

            ch_modif = DORADO_DOWNLOAD_MODEL.out.model
                .filter { type, _model -> type == "modif" }
                .map { _type, model -> model }

        // download only model
        } else if (!model_resolved && !params.m_bases) {

            DORADO_DOWNLOAD_LIST()

            ch_model_to_download = selectBaseModelDownload(
                DORADO_DOWNLOAD_LIST.out.list,
                params.basecall_model
            )
            .combine(ch_model_dir)
            .map { type, model -> tuple("base", type, model)}

            ch_in_download = ch_model_to_download

            DORADO_DOWNLOAD_MODEL(ch_model_to_download)

            ch_model = DORADO_DOWNLOAD_MODEL.out.model
                .map { _type, model -> model }

        // download only modif
        } else if (model_resolved && !modif_resolved && params.m_bases) {

            DORADO_DOWNLOAD_LIST()

            ch_model_match = ch_model.map { it -> it.name }

            ch_modif_to_download = selectModifiedModelDownload(
                DORADO_DOWNLOAD_LIST.out.list,
                ch_model_match,
                params.m_bases
            )
                .combine(ch_model_dir)
                .map { type, model -> tuple("modif", type, model)}

            DORADO_DOWNLOAD_MODEL(ch_modif_to_download)

            ch_modif = DORADO_DOWNLOAD_MODEL.out.model
                .map { _type, model -> model }
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENSEMBLVEP CACHE CHECKUP AND DOWNLOAD
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Use a supplied VEP cache when it contains the selected assembly. Otherwise,
    // download the exact cache required by the selected reference genome.
    def vepCacheDetails = getVepCacheDetails(params.genome)
    def refCacheFile = params.ref_cache ? file(params.ref_cache) : null
    def refCacheDir  = refCacheFile ? (refCacheFile.isDirectory() ? refCacheFile.toString() : refCacheFile.getParent()) : null

    def vepCacheCandidates = params.vep_cache
        ? [params.vep_cache]
        : (refCacheDir ? ['vep', 'vep_cache'].collect { cacheName -> "${refCacheDir}/${cacheName}" } : [])

    def vepCacheDirectory = vepCacheCandidates.find { candidate ->
        hasVepCache(candidate, vepCacheDetails, params.vep_version)
    }

    if (vepCacheDirectory) {
        ch_vep_cache = channel.fromPath(vepCacheDirectory, checkIfExists: true).collect()
        println "Using existing Ensembl VEP ${params.vep_version} cache for ${vepCacheDetails.assembly} from ${vepCacheDirectory}"
    } else {
        log.info("Downloading Ensembl VEP ${params.vep_version} cache for ${vepCacheDetails.assembly}")

        ch_vep_cache_to_download = channel.of([
            [id: "${params.vep_version}_${vepCacheDetails.assembly}"],
            vepCacheDetails.assembly,
            vepCacheDetails.species,
            params.vep_version
        ])

        ENSEMBLVEP_DOWNLOAD(ch_vep_cache_to_download, true)

        ch_vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache
            .map { _meta, cache -> cache }
            .collect()
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Error messages
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.demux) {
        channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map {
                meta, project, _input, _ubam, kit, barcode, _bed, _padding, _low_fidelity, _purity, _filter, _tumor_type ->
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
                _meta, _project, _input, _ubam, kit, barcode, _bed, _padding, _low_fidelity, _purity, _filter, _tumor_type ->
                if (kit || barcode) {
                    throw new IllegalArgumentException("Please use --demux option for demultiplexing samples")
                }
            }
    }

    emit:
    bed_sheet   = ch_adaptive
    bed_nopad   = ch_bednopad
    demux_sheet = ch_demux
    samplesheet = ch_in_samplesheet
    tumor_type  = ch_tumor_type
    ref_fa      = ch_ref_fasta
    ref_fai     = ch_ref_fai
    cfdna_ch    = ch_cfdna
    vep_cache   = ch_vep_cache
    model_ch    = ch_model
    modif_ch    = ch_modif
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

    if (operation == 'remove_suffix') {
        if (new_meta.id && suffix) {
            new_meta.id = new_meta.id.replaceFirst(suffix.toString(), '')
        }
    }
    else if (operation == 'add_suffix') {
        if (new_meta.id && suffix) {
            new_meta.id = new_meta.id + suffix
        }
    }
    else if (operation == 'replace') {
        if (new_meta.id && search_string) {
            new_meta.id = new_meta.id.replace(search_string, replace_string ?: '')
        }
    }
    else if (operation == 'prefix') {
        if (new_meta.id && suffix) {
            new_meta.id = suffix + new_meta.id
        }
    }

    return new_meta
}

def getMinQC(model) {
    def model_string = model.toString()
    if (model_string.contains('sup')) {
        return 10
    }
    else if (model_string.contains('hac')) {
        return 9
    }
    else {
        return 8
    }
}

def selectLatestModel(model, model_path) {
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
        f.name ==~ ".*${model_name}@v${version_pattern}\$"
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
    }
    else {
        if (requested_version) {
            return clean_model_dna.first()
        }

        // 2) Select latest clean model
        return clean_model_dna
            .sort { a, b ->
                def va = a.name.split('@v')[-1].tokenize('.').collect { it -> it as int }
                def vb = b.name.split('@v')[-1].tokenize('.').collect { it -> it as int }

                int len = Math.max(va.size(), vb.size())

                def diff_index = (0..<len).find { i ->
                    def ai = i < va.size() ? va[i] : 0
                    def bi = i < vb.size() ? vb[i] : 0
                    ai != bi
                }

                if (diff_index == null) return 0

                def ai = diff_index < va.size() ? va[diff_index] : 0
                def bi = diff_index < vb.size() ? vb[diff_index] : 0
                ai <=> bi
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
        f.name ==~ "${base_model}_${modif_name}@v${version_pattern}\$"
    }

    if (!clean_model_mod) {
        def version_note = requested_version ? " (requested version ${requested_version})" : ""
        log.warn(
            "No modified model found for ${base_model} with ${modif_name}${version_note} in ${model_path}: model will be downloaded"
        )
        return null
    }
    else {
        if (requested_version) {
            return clean_model_mod.first()
        }
        return clean_model_mod
            .sort { a, b ->
                def va = a.name.split('@v')[-1].tokenize('.').collect { it -> it.toInteger() }.join('').toInteger()
                def vb = b.name.split('@v')[-1].tokenize('.').collect { it -> it.toInteger() }.join('').toInteger()
                va <=> vb
            }
            .last()
    }
}

def normalizeVersion(versionString) {
    versionString
        .tokenize('.')
        .collect { item -> String.format('%03d', item as int) }
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
            .filter { it -> it.key.contains(modelParam) }
            .map { entry ->
                def v = (entry.key =~ /@v([0-9.]+)/)[0][1]
                tuple(entry.key, normalizeVersion(v))
            }
            .ifEmpty { error "No models found for: ${modelParam}" }
            .toSortedList { a, b -> a[1] <=> b[1] }
            .map { list ->
                def exact = list.find { item -> item[0] == modelParam }
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
        .map { entry, _parent ->
            (entry.value[subfield] ?: [:]).entrySet()
        }
        .flatMap { entrySet -> entrySet }
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
            def exact = list.find { item -> item[0] == modelParam }
            exact ? exact[0] : list[-1][0]
        }
}

def selectBaseModelDownload(chModelsList, modelParam) {
    selectModelDownload(chModelsList, modelParam)
}

def selectModifiedModelDownload(chModelsList, chBaseModel, modifParam) {
    selectModelDownload(chModelsList, modifParam, 'modified_models', chBaseModel)
}

def getVepCacheDetails(genome) {
    def assemblies = [
        hg38: [assembly: 'GRCh38', species: 'homo_sapiens'],
        hg19: [assembly: 'GRCh37', species: 'homo_sapiens'],
        hs1 : [assembly: 'T2T-CHM13v2.0', species: 'homo_sapiens_gca009914755v4']
    ]

    def aliases = [
        hg38  : 'hg38',
        GRCh38: 'hg38',
        hg19  : 'hg19',
        GRCh37: 'hg19',
        hs1   : 'hs1',
        CHM13 : 'hs1',
        t2t   : 'hs1'
    ]

    def canonicalGenome = aliases[genome]
    if (!canonicalGenome) {
        throw new IllegalArgumentException("Unsupported genome for VEP cache: ${genome}")
    }

    assemblies[canonicalGenome]
}

def hasVepCache(cacheDirectory, cacheDetails, cacheVersion) {
    cacheDirectory && file("${cacheDirectory}/${cacheDetails.species}/${cacheVersion}_${cacheDetails.assembly}").isDirectory()
}

// Case-insensitive FASTA/FAI finders
def findFaFile (files, aliasList, exts) {
    files.find { f -> aliasList.any { a -> f.name.toLowerCase().contains(a.toLowerCase()) } && exts.any { ext -> f.name.toLowerCase().endsWith(ext) } }
}

def findCached(String url) {
    if( !params.ref_cache ) return null
    def cache = file(params.ref_cache)
    if( !cache.isDirectory() ) return null

    def cached = cache.resolve(url.tokenize('/')[-1])
    return cached.exists() ? cached : null
}

def hasVal(v) {
    return v != null && v != '' && v != []
}

// Needed by the BED decompression step
def isGz(bed) {
    return bed.name ==~ /.*\.b?gz$/
}

def bedMeta(bed) {
    return [ id: "${bed.simpleName}_${Math.abs(bed.toString().hashCode())}" ]
}
