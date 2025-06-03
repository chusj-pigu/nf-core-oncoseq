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
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oncoseq_pipeline'

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_ONCOSEQ_ADAPTIVE {

    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    bed         // channel: from path read from params.bed, bed file used for adaptive sampling
    clairs_model  // channel: model for calling snp with ClairS-TO
    basecall_model  // channel : basecalling model used with dorado
    ch_clin_database            // channel : from path, vcf containing the ClinVar database for annotating vcf

    main:

    //
    // WORKFLOW: Run pipeline
    //
    ADAPTIVE (
        samplesheet,
        demux,
        ref,
        bed,
        clairs_model,
        basecall_model,
        ch_clin_database

    )
}

workflow NFCORE_ONCOSEQ_CFDNA {

    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path

    main:

    //
    // WORKFLOW: Run pipeline
    //
    CFDNA (
        samplesheet,
        demux,
        ref
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
        ch_clin_database
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
        params.ubam_samplesheet,
        params.demux_samplesheet,
        params.adaptive_samplesheet,
        params.bed,
        params.padding,
        params.low_fidelity
    )

    // Load Channels from parameters:

    ch_model = params.model ? Channel.of(params.model) : Channel.fromPath(params.model_path)

    // Combine the samplesheet with the model :
    if (params.skip_basecalling) {
        ch_input = PIPELINE_INITIALISATION.out.samplesheet
    } else if (params.ubam_samplesheet == null ) {
        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .combine(PIPELINE_INITIALISATION.out.ubam_ch)
            .combine(ch_model)
    } else {
        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .combine(ch_model)
    }

   // Channels for SNP calling
    ch_chr_list = Channel.of(params.chr_list)
    ch_clairs_model = Channel.of(params.clairsto_model)
    ch_clin_database = Channel.fromPath(params.clin_database)


   // WORKFLOW: Run main workflow


    if ( params.adaptive) {
        NFCORE_ONCOSEQ_ADAPTIVE (
        ch_input,
        PIPELINE_INITIALISATION.out.demux_sheet,
        PIPELINE_INITIALISATION.out.ref_ch,
        ch_clairs_model,
        ch_model,
        ch_clin_database,
        PIPELINE_INITIALISATION.out.bed_sheet,
        )
    } else if ( params.cfdna ) {
        NFCORE_ONCOSEQ_CFDNA (
            ch_input,
            PIPELINE_INITIALISATION.out.demux_sheet,
            PIPELINE_INITIALISATION.out.ref_ch
        )
    } else if ( params.wgs ) {
        NFCORE_ONCOSEQ_WGS (
        ch_input,
        PIPELINE_INITIALISATION.out.demux_sheet,
        PIPELINE_INITIALISATION.out.ref_ch,
        ch_clairs_model,
        ch_model,
        ch_clin_database
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
