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

include { ADAPTIVE_WGS            } from './workflows/adaptive_wgs'
include { CFDNA                   } from './workflows/cfdna'
include { LOCAL_REALTIME          } from './workflows/local_realtime'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { getMinQC                } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { selectLatestModel       } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { selectLatestModif       } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { selectModelDownload     } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { DORADO_DOWNLOAD_LIST    } from './modules/local/dorado/main.nf'
include { DORADO_DOWNLOAD_MODEL   } from './modules/local/dorado/main.nf'
include { ENSEMBLVEP_DOWNLOAD     } from './modules/nf-core/ensemblvep/download/main.nf'

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
        hs1 : [assembly: 'CHM13', species: 'homo_sapiens']
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

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_ONCOSEQ_ADAPTIVE_WGS {



    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    bed         // channel: from path read from params.bed, bed file used for adaptive sampling
    clairs_model  // channel: model for calling snp with ClairS-TO
    basecall_model  // channel : basecalling model used with dorado
    sv_targets              // channel : list of genes with their position to represent in Figeno
    vep_cache

    main:

    //
    // WORKFLOW: Run pipeline
    //

    if (params.realtime == null) {
        ADAPTIVE_WGS (
            samplesheet,
            demux,
            ref,
            bed,
            clairs_model,
            basecall_model,
            sv_targets,
            vep_cache
        )
    } else {
        LOCAL_REALTIME (
            samplesheet,
            demux,
            ref,
            bed,
            basecall_model,
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
    ref                 // channel : reference for mapping, either empty if skipping mapping, or a path
    max_len
    minqs
    ichor_bin
    mapq_wig
    basecall_model          // channel: model for basecalling
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
        ref,
        max_len,
        minqs,
        ichor_bin,
        mapq_wig,
        basecall_model,
        clairs_model,
        bed,
        targets,
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
        params.bed,
        params.padding,
        params.low_fidelity
    )

    // Load model channels from parameters:
    if (params.skip_basecalling || params.skip_mapping ) {

        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .map { meta, input, _ubam ->
            tuple(meta, input) }

        ch_model = channel.of(params.basecall_model)

    } else {

        if (!params.m_bases) {
            log.warn("Tumor classification will be skipped — parameter --m_bases not set")
        }

        ch_model_dir = params.ref_cache ? channel.fromPath("${params.ref_cache}", checkIfExists:true) : channelfromPath("${launchDir}")
        def model_resolved = params.ref_cache  ? selectLatestModel(params.basecall_model,file("${params.ref_cache}/dorado_models")) : null
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
                .filter { type, model -> type == "base" }
                .map { type, model -> model }

            ch_modif = DORADO_DOWNLOAD_MODEL.out.model
                .filter { type, model -> type == "modif" }
                .map { type, model -> model }

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
                .map { type, model -> model }

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
                .map { type, model -> model }
        }

        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .combine(ch_model)
            .combine(ch_modif)
    }

   // channels for SNP calling
    ch_clairs_model = channel.of(params.clairsto_model)

    // Use a supplied VEP cache when it contains the selected assembly. Otherwise,
    // download the exact cache required by the selected reference genome.
    def vepCacheDetails = getVepCacheDetails(params.genome)
    def vepCacheDirectory = params.vep_cache ?: (
        params.ref_cache && file(params.ref_cache).isDirectory() ? "${params.ref_cache}/vep" : null
    )

    if (hasVepCache(vepCacheDirectory, vepCacheDetails, params.vep_version)) {
        ch_vep_cache = channel.fromPath(vepCacheDirectory, checkIfExists: true).collect()
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
    // channel for sv gene targets
    ch_sv_targets = channel.fromPath(params.sv_targets)

    // channels for chopper and IchorCNA (cfDNA)

    ch_max_len    = channel.of(params.max_length)
    ch_minqs      = ch_model.map { model -> getMinQC(model) }
    ch_ichor_bin  = channel.of(params.ichor_bin_size)
    ch_min_mapq   = channel.of(params.min_mapq_ichor)


   // WORKFLOW: Run main workflow

   if ( params.cfdna ) {
        NFCORE_ONCOSEQ_CFDNA (
            ch_input,
            PIPELINE_INITIALISATION.out.demux_sheet,
            PIPELINE_INITIALISATION.out.cfdna_ch,
            PIPELINE_INITIALISATION.out.ref_ch,
            ch_max_len,
            ch_minqs,
            ch_ichor_bin,
            ch_min_mapq,
            ch_model,
            ch_clairs_model,
            PIPELINE_INITIALISATION.out.bed_sheet,
            ch_sv_targets,
            ch_vep_cache
        )

    } else {
        NFCORE_ONCOSEQ_ADAPTIVE_WGS (
            ch_input,
            PIPELINE_INITIALISATION.out.demux_sheet,
            PIPELINE_INITIALISATION.out.ref_ch,
            ch_clairs_model,
            ch_model,
            PIPELINE_INITIALISATION.out.bed_sheet,
            ch_sv_targets,
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
