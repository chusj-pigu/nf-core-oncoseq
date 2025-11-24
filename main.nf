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
    ch_id

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
            ch_id
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
            ch_id
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
    ch_id

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
        ch_id
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
    ch_id

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
        ch_id
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

    // Load model channels from parameters:

    ch_model = params.basecall_model ? channel.of(params.basecall_model) : channel.fromPath(params.basecall_model_path)
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
    ch_clin_database = channel.fromPath(params.clin_database)

    // channel for sv gene targets
    ch_sv_targets = channel.fromPath(params.sv_targets)

    // channels for chopper and IchorCNA (cfDNA)

    ch_max_len    = channel.of(params.max_length)
    ch_minqs      = ch_model.map { model -> getMinQC(model) }
        .view()
    ch_ichor_bin  = channel.of(params.ichor_bin_size)
    ch_min_mapq   = channel.of(params.min_mapq_ichor)

    // channel for T2T reference:

    ch_cns = PIPELINE_INITIALISATION.out.tumor_type
        .branch { meta, tumor ->
        tumor: tumor == "cns"
        }

    ch_ref_t2t_id = channel.of("t2t", "no_index")
        .toList()
    ch_ref_t2t = channel.fromPath(params.ref_t2t)
        .combine(ch_ref_t2t_id)
        .combine(ch_cns.tumor)
        .map { ref_path, ref_id, ref_index, meta, _tumor ->
            tuple(meta, ref_id, ref_path, ref_index)}
        .transpose()


   // WORKFLOW: Run main workflow


    if ( params.adaptive) {
        NFCORE_ONCOSEQ_ADAPTIVE (
            ch_input,
            PIPELINE_INITIALISATION.out.demux_sheet,
            PIPELINE_INITIALISATION.out.ref_ch,
            PIPELINE_INITIALISATION.out.tumor_type,
            ch_ref_t2t,
            ch_clairs_model,
            ch_model,
            ch_clin_database,
            PIPELINE_INITIALISATION.out.bed_sheet,
            ch_sv_targets,
            PIPELINE_INITIALISATION.out.id_ch
        )
    } else if ( params.cfdna ) {
        NFCORE_ONCOSEQ_CFDNA (
            ch_input,
            PIPELINE_INITIALISATION.out.demux_sheet,
            PIPELINE_INITIALISATION.out.cfdna_ch,
            PIPELINE_INITIALISATION.out.tumor_type,
            PIPELINE_INITIALISATION.out.ref_ch,
            ch_ref_t2t,
            ch_max_len,
            ch_minqs,
            ch_ichor_bin,
            ch_min_mapq,
            PIPELINE_INITIALISATION.out.id_ch
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
            PIPELINE_INITIALISATION.out.id_ch
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
