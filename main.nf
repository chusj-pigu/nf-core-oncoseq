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

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_ONCOSEQ_ADAPTIVE_WGS {



    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from samplesheet
    tumor_type  // channel: tumor type read in from samplesheet
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
            tumor_type,
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
            tumor_type,
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
    tumor_type          // channel: tumor type read in from samplesheet
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
        tumor_type,
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

params {

    // Bin size to use for CNV calling with QDNAseq (Kbp)
    qdnaseq_binsize: Integer = 500

    // Bed file used for adaptive mode (default is latest CRA-CHUSJ panel)
    bed: Path = "${projectDir}/assets/v2.0.1-pre-merge-panel-20kb-pad.bed"

    // turn on multiplex basecalling
    demux: Boolean

    // Modifications model to usefor basecalling
    m_bases: String?

    // Device to use for basecalling
    device: String?

    // Batch size for basecalling
    batch: String?

    // Base path for illumina genomes to pull from AWS, typically not changed
    igenomes_base: String = 's3://ngi-igenomes/igenomes/'

    // Option for finetuning time ressources allocation for minimap2. Set to false if a big genome is to be aligned (human whole genome 30X)
    mapping_small: Boolean = true

    // Run in cell-free DNA mode. Cannot be set with adaptive sampling or wgs mode.
    cfdna: Boolean

    // Scale resource allocations for very large datasets.
    huge: Boolean

    // ClairS_TO model to use for SNP calling
    clairsto_model: String = 'ont_r10_dorado_sup_5khz_ssrs'

    // Path to text file containing list of genes of low fidelity for adaptive sampling (default is empty list)
    low_fidelity: Path = "${projectDir}/assets/panel_low_fidelity.txt"

    // Padding around region of interest for adaptive sampling
    padding: Integer = 20000

    // Option to skip basecalling and start with fastq files in the samplesheet input column
    skip_basecalling: Boolean

    // Maximum memory any process can consume
    max_memory: MemoryUnit

    // Maximum cpus any process can use
    max_cpus: Integer

    // Maximum amount of time any process can take
    max_time: Duration

    time_series: Boolean

    time_points: String = '1,3,6,12,18,24,48,60,72'

    include_full: Boolean = true

    basecall_model: String

    // Use number of sequencing hour to activate real time mode
    realtime: Integer?

    subchrom_binsize: Integer = 500

    // List of genes to visualize with Figeno (columns GENE and pos)
    sv_targets: Path = "${projectDir}/assets/sv-list.csv"

    // List of important fusions to verify with stellerator
    fusion_targets: Path = "${projectDir}/assets/fusion-list.csv"

    // List of genes to exclude from SV analysis (possible artefacts)
    sv_exclude: Path = "${projectDir}/assets/sv_exclude.txt"

    // List of genes to exclude from snp analysis
    snp_exclude: Path = "${projectDir}/assets/snp_exclude.txt"

    skip_mapping: Boolean

    // -- ichorCNA parameters --------------------------------------------------

    // Maximum length (in bases) of reads to include in cfdna analysis (if filter true in samplesheet)
    max_length: Integer = 700

    // Minimum MAPQ of reads to include in ichorCNA analysis
    min_mapq_ichor: Integer = 20

    // Bin size to use for IchorCNA (in bases)
    ichor_bin_size: Integer = 500000

    // Comma-separated chromosome list used by readCounter to build the wig file
    chr_wig: String = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"

    // Candidate tumour ploidy values passed to run IchorCNA.R (R vector syntax)
    custom_ploidy: String = 'c(2,3)'

    // Maximum copy number to consider
    custom_maxCN: Integer = 5

    // Chromosome naming style ichorCNA should use internally (UCSC or Ensembl)
    genome_style: String = "UCSC"

    // Estimate subclonal prevalence within IchorCNA (default is false, which is faster)
    estimate_sc_prevalence: Boolean

    // ---------------------------------------------------------------------------

    report_template: Path = "${projectDir}/assets/templates/report_template.qmd"

    // Placeholder for incomplete ubam, to use basecalling resume, specify path to incomplete ubam in samplesheet
    ubam: Path = "${projectDir}/assets/NOFILE"

    // Bin size to use with Delly CNV calling
    delly_bin_size: Integer = 50000

    // Optional root directory for Ensembl VEP caches. The selected version and genome assembly are downloaded automatically on the Nextflow launch host when absent; that host must have internet access.
    vep_cache: Path?

    // Version of Emsembl to use
    vep_version: Integer = 115

    // Expression accepted by filtervep for filtering snp
    filtervep_expression: String = 'IMPACT != LOW and IMPACT != MODIFIER and (gnomADe_AF <= 0.01 or not gnomADe_AF) and not CLIN_SIG matches benign and MANE'

    // Path to directory containing reference genome
    ref_cache: Path?

    // Human genome assembly identifier
    genome: String = 'GRCh38'

    // Use gpu to run Clair3
    clair3_gpu: Boolean

    // Path to comma-separated file containing information about the samples in the experiment.
    input: Path

    // The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.
    outdir: Path

    // Email address for completion summary.
    email: String?

    // Git commit id for Institutional configs.
    custom_config_version: String = 'master'

    // Base directory for Institutional configs.
    custom_config_base: String = 'https://raw.githubusercontent.com/nf-core/configs/master'

    // Institutional config name.
    config_profile_name: String?

    // Institutional config description.
    config_profile_description: String?

    // Institutional config contact information.
    config_profile_contact: String?

    // Institutional config URL link.
    config_profile_url: String?

    // Display version and exit.
    version: Boolean

    // Method used to save pipeline results to output directory.
    publish_dir_mode: String = 'copy'

    // Email address for completion summary, only when pipeline fails.
    email_on_fail: String?

    // Send plain-text email instead of HTML.
    plaintext_email: Boolean

    // Do not use coloured log outputs.
    monochrome_logs: Boolean

    // Incoming hook URL for messaging service
    hook_url: String?

    // Boolean whether to validate parameters against the schema at runtime
    validate_params: Boolean = true

    // Base URL or local path to location of pipeline test dataset files
    pipelines_testdata_base_path: Path = "${projectDir}/test_data/"

    // Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss.
    trace_report_suffix: String?

    // Display full help information
    help_full: Boolean

    // DOI for the pipeline, used in the help message
    pipeline_doi: String = ""

    // Text to display before the help message
    before_text: String = "nf-core-oncoseq - Nextflow pipeline for realtime and non-realtime oncology analysis of adaptive sampling sequencing and cfDNA"

    // Text to display after the help message
    after_text: String = null

    // Version of the pipeline, used in the help message
    pipeline_version: String = 'pre-release'

    // Show hidden parameters in the help message
    show_hidden: Boolean

    wgs: Boolean
}

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
        params.basecall_model,
        params.ref_cache,
        params.vep_cache
    )

    // Load model channels from parameters:
    if (params.skip_basecalling || params.skip_mapping ) {

        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .map { meta, input, _ubam ->
            tuple(meta, input) }

        ch_model = channel.of(params.basecall_model)

    } else {

        ch_input = PIPELINE_INITIALISATION.out.samplesheet
            .combine(PIPELINE_INITIALISATION.out.model_ch)
            .combine(PIPELINE_INITIALISATION.out.modif_ch)

        ch_model = PIPELINE_INITIALISATION.out.model_ch
    }

   // channels for SNP calling
    ch_clairs_model = channel.of(params.clairsto_model)

    ch_vep_cache = PIPELINE_INITIALISATION.out.vep_cache
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
            PIPELINE_INITIALISATION.out.tumor_type,
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
            PIPELINE_INITIALISATION.out.tumor_type,
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
