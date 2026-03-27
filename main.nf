#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oncoseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/oncoseq
    Website: https://nf-co.re/oncoseq
    Slack  : https://nfcore.slack.com/channels/oncoseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ADAPTIVE                } from './workflows/adaptive'
include { CFDNA                   } from './workflows/cfdna'
include { WGS                     } from './workflows/wgs'
include { LOCAL_REALTIME          } from './workflows/local_realtime'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { getMinQC                } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { select_latest_model     } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { select_latest_modif     } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { DORADO_DOWNLOAD_LIST    } from './modules/local/dorado/main.nf'
include { DORADO_DOWNLOAD_MODEL   } from './modules/local/dorado/main.nf'
include { STAGE_REFERENCE_FILES as STAGE_CLINICAL_REFERENCE_FILES } from './modules/local/reference_cache/main.nf'

def resolveIndexedVcfFiles(vcfSpec) {
    def resolved = file(vcfSpec, checkIfExists: true)
    def resolvedFiles = resolved instanceof List ? resolved : [resolved]
    def vcf = resolvedFiles.find { resolvedFile ->
        resolvedFile.name.endsWith('.vcf.gz') || resolvedFile.name.endsWith('.vcf')
    } ?: resolvedFiles.find { resolvedFile ->
        !resolvedFile.name.endsWith('.tbi')
    }

    if (!vcf) {
        throw new IllegalArgumentException(
            "Could not resolve a VCF file from clinical database path: ${vcfSpec}"
        )
    }

    def index = resolvedFiles.find { resolvedFile ->
        resolvedFile.name.endsWith('.tbi')
    } ?: file("${vcfSpec}.tbi", checkIfExists: true)

    tuple(vcf, index)
}

def normalizeStagedFiles(stagedFiles) {
    stagedFiles instanceof List ? stagedFiles : [stagedFiles]
}

def normalizeVersion(versionString) {
    versionString
        .tokenize('.')
        .collect { String.format('%03d', it as int) }
        .join()
}

def selectModel(chModelsList, modelParam, subfield = null, chParentModel = null) {

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

def selectBaseModel(chModelsList, modelParam) {
    selectModel(chModelsList, modelParam)
}

def selectModifiedModel(chModelsList, chBaseModel, modifParam) {
    selectModel(chModelsList, modifParam, 'modified_models', chBaseModel)
}

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_ONCOSEQ_ADAPTIVE {



    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    tumor_type  // channel: samplesheet read in from --input, contains only tumor type
    ch_ref_t2t  // channel : from params.ref_t2t
    bed         // channel: from path read from params.bed, bed file used for adaptive sampling
    clairs_model  // channel: model for calling snp with ClairS-TO
    basecall_model  // channel : basecalling model used with dorado
    ch_clin_database            // channel : from path, vcf containing the ClinVar database for annotating vcf
    sv_targets              // channel : list of genes with their position to represent in Figeno
    vep_cache

    main:

    //
    // WORKFLOW: Run pipeline
    //

    if (params.realtime == null) {
        ADAPTIVE (
            samplesheet,
            demux,
            ref,
            bed,
            clairs_model,
            basecall_model,
            ch_clin_database,
            sv_targets,
            tumor_type,
            ch_ref_t2t,
            vep_cache
        )
    } else {
        LOCAL_REALTIME (
            samplesheet,
            demux,
            ref,
            tumor_type,
            ch_ref_t2t,
            bed,
            basecall_model,
            ch_clin_database,
            sv_targets,
            vep_cache
        )
    }
}

workflow NFCORE_ONCOSEQ_CFDNA {

    take:
    samplesheet         // channel: samplesheet read in from --input
    demux               // channel: demux_samplesheet read in from --demux_samplesheet
    cfdna_samplesheet   // channel : from demux or samplesheeet
    tumor_type
    ref                 // channel : reference for mapping, either empty if skipping mapping, or a path
    ref_t2t
    max_len
    minqs
    ichor_bin
    mapq_wig
    basecall_model          // channel: model for basecalling
    ch_clin_database        // channel: clinical database for variant annotation
    clairs_model
    bed                     // channel: bed file used for adaptive sampling regions
    targets
    vep_cache

    main:

    //
    // WORKFLOW: Run pipeline
    //
    CFDNA (
        samplesheet,
        demux,
        cfdna_samplesheet,
        tumor_type,
        ref,
        ref_t2t,
        max_len,
        minqs,
        ichor_bin,
        mapq_wig,
        basecall_model,
        ch_clin_database,
        clairs_model,
        bed,
        targets,
        vep_cache
    )
}

workflow NFCORE_ONCOSEQ_WGS {

    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    clairs_model
    basecall_model
    ch_clin_database
    bed_empty       // channel with empty bed file to trigger subchrom
    sv_targets
    ch_minqs
    tumor_type
    ref_t2t
    vep_cache

    main:

    //
    // WORKFLOW: Run pipeline
    //
    WGS (
        samplesheet,
        demux,
        ref,
        clairs_model,
        basecall_model,
        ch_clin_database,
        bed_empty,
        sv_targets,
        ch_minqs,
        tumor_type,
        ref_t2t,
        vep_cache
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    if (!params.reference_cache_dir) {
        throw new IllegalArgumentException(
            'Please provide --reference_cache_dir to stage reference assets ' +
                'and Dorado models.'
        )
    }

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.ubam_samplesheet,
        params.demux_samplesheet,
        params.adaptive_samplesheet,
        params.bed,
        params.padding,
        params.low_fidelity
    )

    // Load model channels from parameters:
    if (!params.skip_basecalling && params.basecall_offline) {
        def model_resolved = select_latest_model(params.basecall_model,file("${params.reference_cache_dir}/dorado_models"))
        def modif_resolved = params.m_bases && model_resolved ?
            select_latest_modif(file("${params.reference_cache_dir}/dorado_models"), model_resolved, params.m_bases) :
            null
        ch_modif = channel.of('')

        if (model_resolved) {
            ch_model = channel.of(model_resolved)
        }

        if (modif_resolved) {
            ch_modif = channel.of(modif_resolved)
        }

        // download both
        if (!model_resolved && !modif_resolved && params.m_bases) {

            DORADO_DOWNLOAD_LIST()

            ch_model_to_download = selectBaseModel(
                DORADO_DOWNLOAD_LIST.out.list,
                params.basecall_model
            )

            ch_modif_to_download = selectModifiedModel(
                DORADO_DOWNLOAD_LIST.out.list,
                ch_model_to_download,
                params.m_bases
            )

            ch_in_download = ch_model_to_download
                .mix(ch_modif_to_download)

            DORADO_DOWNLOAD_MODEL(ch_in_download)

            ch_model = DORADO_DOWNLOAD_MODEL.out.model
                .combine(ch_model_to_download)
                .filter { model_out, model_in -> model_out.name == model_in }
                .map { it -> it[0].name }

            ch_modif = DORADO_DOWNLOAD_MODEL.out.model
                .combine(ch_modif_to_download)
                .filter { model_out, model_in -> model_out.name == model_in }
                .map { it -> it[0].name }

        // download only model
        } else if (!model_resolved && !params.m_bases) {

            DORADO_DOWNLOAD_LIST()

            ch_model_to_download = selectBaseModel(
                DORADO_DOWNLOAD_LIST.out.list,
                params.basecall_model
            )

            ch_in_download = ch_model_to_download

            DORADO_DOWNLOAD_MODEL(ch_in_download)

            ch_model = DORADO_DOWNLOAD_MODEL.out.model
                .combine(ch_in_download)
                .filter { model_out, model_in -> model_out.name == model_in }
                .map { it -> it[0].name }

        // download only modif
        } else if (model_resolved && !modif_resolved && params.m_bases) {

            DORADO_DOWNLOAD_LIST()

            ch_modif_to_download = selectModifiedModel(
                DORADO_DOWNLOAD_LIST.out.list,
                ch_model,
                params.m_bases
            )

            ch_in_download = ch_modif_to_download

            DORADO_DOWNLOAD_MODEL(ch_in_download)

            ch_modif = DORADO_DOWNLOAD_MODEL.out.model
                .combine(ch_in_download)
                .filter { model_out, model_in -> model_out.name == model_in }
                .map { it -> it[0].name }
        }

    ch_model_dir = channel.fromPath("${params.reference_cache_dir}/dorado_models")
        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .combine(ch_model_dir)
            .combine(ch_model)
            .combine(ch_modif)

    } else if (!params.skip_basecalling && !params.basecall_offline) {
        ch_model = channel.of(params.basecall_model)
        ch_modif = channel.of(params.m_bases)
        def model_dir = '.'

        if (params.reference_cache_dir) {
            def base = file(params.reference_cache_dir)
            def sub  = file("${params.reference_cache_dir}/dorado_models")

            if (sub.exists()) {
                model_dir = sub
            }
            else if (base.exists()) {
                model_dir = base
            }
        }

        ch_model_dir = channel.fromPath(model_dir)

        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .combine(ch_model_dir)
            .combine(ch_model)
            .combine(ch_modif)
    } else {
        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .map { meta, input, _ubam ->
            tuple(meta, input) }
    }

   // channels for SNP calling
    ch_clairs_model = channel.of(params.clairsto_model)

    // vep cache
    ch_vep_cache = Channel.fromPath(params.vep_cache, checkIfExists: true).collect()

    // channel for sv gene targets
    ch_sv_targets = channel.fromPath(params.sv_targets)

    // channels for chopper and IchorCNA (cfDNA)

    ch_max_len    = channel.of(params.max_length)
    ch_minqs      = ch_model.map { model -> getMinQC(model) }
    ch_ichor_bin  = channel.of(params.ichor_bin_size)
    ch_min_mapq   = channel.of(params.min_mapq_ichor)


   // WORKFLOW: Run main workflow


    if ( params.adaptive) {
        NFCORE_ONCOSEQ_ADAPTIVE (
            ch_input,
            PIPELINE_INITIALISATION.out.demux_sheet,
            PIPELINE_INITIALISATION.out.ref_ch,
            PIPELINE_INITIALISATION.out.tumor_type,
            PIPELINE_INITIALISATION.out.ref_t2t,
            ch_clairs_model,
            ch_model,
            ch_clin_database,
            PIPELINE_INITIALISATION.out.bed_sheet,
            ch_sv_targets,
            ch_vep_cache
        )
    } else if ( params.cfdna ) {
        NFCORE_ONCOSEQ_CFDNA (
            ch_input,
            PIPELINE_INITIALISATION.out.demux_sheet,
            PIPELINE_INITIALISATION.out.cfdna_ch,
            PIPELINE_INITIALISATION.out.tumor_type,
            PIPELINE_INITIALISATION.out.ref_ch,
            PIPELINE_INITIALISATION.out.ref_t2t,
            ch_max_len,
            ch_minqs,
            ch_ichor_bin,
            ch_min_mapq,
            ch_model,
            ch_clin_database,
            ch_clairs_model,
            PIPELINE_INITIALISATION.out.bed_sheet,
            ch_sv_targets,
            ch_vep_cache
        )
    } else if ( params.wgs ) {
        NFCORE_ONCOSEQ_WGS (
            ch_input,
            PIPELINE_INITIALISATION.out.demux_sheet,
            PIPELINE_INITIALISATION.out.ref_ch,
            ch_clairs_model,
            ch_model,
            ch_clin_database,
            PIPELINE_INITIALISATION.out.bed_sheet,
            ch_sv_targets,
            ch_minqs,
            PIPELINE_INITIALISATION.out.tumor_type,
            PIPELINE_INITIALISATION.out.ref_t2t,
            ch_vep_cache
        )
    }
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
