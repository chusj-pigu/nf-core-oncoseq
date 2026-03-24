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

def getDoradoModelsDirectory() {
    params.reference_cache_dir
        ? file(params.reference_cache_dir)
            .toAbsolutePath()
            .resolve('dorado_models')
            .toString()
        : null
}

def getCachedDoradoModelPath(modelName) {
    def modelsDirectory = getDoradoModelsDirectory()
    (modelsDirectory && modelName) ? "${modelsDirectory}/${modelName}" : null
}

def cachedDoradoModelExists(modelName) {
    def cachedModelPath = getCachedDoradoModelPath(modelName)
    cachedModelPath ? file(cachedModelPath).exists() : false
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

    def doradoModelsDirectory = getDoradoModelsDirectory()
    if (
        params.basecall_model &&
        !params.basecall_model_path &&
        !params.skip_basecalling &&
        doradoModelsDirectory
    ) {
        if (cachedDoradoModelExists(params.basecall_model)) {
            ch_model = channel.of(
                getCachedDoradoModelPath(params.basecall_model)
            )
        } else {
            ch_model = channel.of(
                tuple([id: 'dorado_model'], params.basecall_model,
                    doradoModelsDirectory)
            )

            DORADO_DOWNLOAD_MODEL(ch_model)

            ch_model = DORADO_DOWNLOAD_MODEL.out.model
                .map { _meta, cachedModelPath ->
                    cachedModelPath.toString()
                }
        }
    } else {
        ch_model = params.basecall_model
            ? channel.of(params.basecall_model)
            : channel.fromPath(params.basecall_model_path)
    }
    ch_modif = channel.of(params.m_bases)

    // Combine the samplesheet with the model :
    if (params.skip_basecalling || params.skip_mapping) {
        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .map { meta, input, _ubam ->
            tuple(meta, input) }
    } else {
        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .combine(ch_model)
            .combine(ch_modif)
    }

   // channels for SNP calling
    ch_clairs_model = channel.of(params.clairsto_model)
    def clinvarDefaultUrl = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz'
    def clinDatabasePath = (!params.clin_database || params.clin_database == 'ClinVar') ? clinvarDefaultUrl : params.clin_database
    ch_clin_database = channel.of(clinDatabasePath)
            .map { dbSpec ->
                def (vcf, index) = resolveIndexedVcfFiles(dbSpec)
                tuple([id: 'clinvar'], 'clinvar', [vcf, index])
            }

    STAGE_CLINICAL_REFERENCE_FILES(ch_clin_database)

    ch_clin_database = STAGE_CLINICAL_REFERENCE_FILES.out.staged
            .map { _meta, _label, stagedFiles ->
                def files = normalizeStagedFiles(stagedFiles)
                def vcf = files.find { stagedFile ->
                    stagedFile.name.endsWith('.vcf.gz') ||
                        stagedFile.name.endsWith('.vcf')
                } ?: files.find { stagedFile ->
                    !stagedFile.name.endsWith('.tbi')
                }
                def index = files.find { stagedFile ->
                    stagedFile.name.endsWith('.tbi')
                }
                tuple(vcf, index)
            }

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
