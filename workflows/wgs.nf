// Basecalling subworkflows
include { BASECALL_SIMPLEX   } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'

// Core analysis subworkflows
include { MAPPING            } from '../subworkflows/local/mapping/mapping'

// Variant calling subworkflows
include { CLAIRS_TO_CALLING                    } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'
include { CLAIR3_CALLING                       } from '../subworkflows/local/variant_calling/clair3_calling.nf'
include { PHASING_VARIANTS as PHASING_SOMATIC  } from  '../subworkflows/local/variant_calling/phasing.nf'
include { PHASING_VARIANTS as PHASING_GERMLINE } from  '../subworkflows/local/variant_calling/phasing.nf'
include { SV_CALLING                           } from '../subworkflows/local/variant_calling/sv_calling.nf'
include { CNV_CALLING                          } from '../subworkflows/local/variant_calling/cnv_calling.nf'
include { SUBCHROM_CALL                        } from '../subworkflows/local/variant_calling/subchrom_call.nf'
include { modifyMetaId                         } from '../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'

// Variant processing and visualization subworkflow
include { VARIANT_PROCESS                       } from  '../subworkflows/local/variant_calling/variant_process.nf'

// Reporting
include { MIDNIGHT_REPORT } from '../subworkflows/local/report/final_report.nf'
include { FIGENO_REPORT   } from '../subworkflows/local/report/variants.nf'
include { WGS_REPORT      } from '../subworkflows/local/report/wgs.nf'

workflow WGS {

    take:
    samplesheet             // channel: samplesheet read in from --input
    demux_samplesheet       // channel : demux samplesheet read in from --demux_samplesheet
    ref                     // channel : reference for mapping, either empty if skipping mapping, or a path
    clairs_model
    basecall_model
    ch_clin_database
    bed_empty
    targets                 // channel : list of genes with their position to represent in Figeno
    minqs                   // channel obtained dynamically from params.basecall_model

    main:

    //
    // WORKFLOW: Run pipeline
    //

    ch_versions = channel.empty()

    if (params.skip_basecalling || params.skip_mapping) {

        MAPPING (
            samplesheet,
            ref
        )

        ch_seqkit = MAPPING.out.seqkit

    } else {

        if (params.demux) {

            BASECALL_MULTIPLEX (
                samplesheet,
                demux_samplesheet
            )

            MAPPING (
                BASECALL_MULTIPLEX.out.fastq,
                ref
            )

            ch_seqkit   = BASECALL_MULTIPLEX.out.stats_pass
            ch_versions = BASECALL_MULTIPLEX.out.versions

        } else {

            BASECALL_SIMPLEX (
                samplesheet
            )

            MAPPING (
                BASECALL_SIMPLEX.out.fastq,
                ref
            )

            ch_seqkit   = BASECALL_SIMPLEX.out.stats_pass
            ch_versions = BASECALL_SIMPLEX.out.versions
        }
    }

        CLAIRS_TO_CALLING (
            MAPPING.out.bam,
            ref,
            clairs_model,
            ch_clin_database
        )

        CLAIR3_CALLING (
            MAPPING.out.bam,
            ref,
            basecall_model,
            ch_clin_database,
            bed_empty
        )

        PHASING_SOMATIC (
            MAPPING.out.bam,
            ref,
            CLAIRS_TO_CALLING.out.vcf
        )

        PHASING_GERMLINE (
            MAPPING.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf
        )

        SV_CALLING (
            PHASING_GERMLINE.out.haptag_bam
                .map { meta, bamfile, bai ->
                    // Restore original sample ID for output naming
                    def meta_restore = modifyMetaId(meta, 'replace', '_somatic_snp_phased', '', '')
                    meta_restore = modifyMetaId(meta_restore, 'replace', '_germline_snp_phased', '', '')
                    tuple(meta_restore, bamfile, bai)
                },
            ref
        )

        CNV_CALLING (
            MAPPING.out.bam,
            ref
        )

        SUBCHROM_CALL (
            MAPPING.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf,
            bed_empty
        )
        // Filter variants to visualize :
        VARIANT_PROCESS (
            MAPPING.out.bam,
            SV_CALLING.out.vcf,
            CNV_CALLING.out.qdnaseq_bed,
            CNV_CALLING.out.qdnaseq_segs,
            targets,
            CNV_CALLING.out.delly_cov,
            CNV_CALLING.out.delly_segs
        )

    ch_versions = ch_versions
        .mix(MAPPING.out.versions)
        .mix(CLAIRS_TO_CALLING.out.versions)
        .mix(CLAIR3_CALLING.out.versions)
        .mix(PHASING_SOMATIC.out.versions)
        .mix(PHASING_GERMLINE.out.versions)
        .mix(SV_CALLING.out.versions)
        .mix(CNV_CALLING.out.versions)
        .mix(VARIANT_PROCESS.out.versions)
        .mix(SUBCHROM_CALL.out.versions)

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPILE SECTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

        WGS_REPORT(
            MAPPING.out.coverage,
            ch_seqkit,
            minqs
        )

        ch_subchrom_plot = SUBCHROM_CALL.out.subchrom_plot_wgs
            //.mix(SUBCHROM_CALL.out.subchrom_plot_panel)

        ch_subchrom_focal = SUBCHROM_CALL.out.subchrom_gene_plot_wgs
            //.mix(SUBCHROM_CALL.out.subchrom_gene_plot_panel)s

        ch_binsize_qdnaseq = channel.of(params.qdnaseq_binsize)
            .map { value ->
            def meta = "qDNAseq"
            tuple(meta, value) }
        ch_binsize_subchrom = channel.of(params.subchrom_binsize)
            .map { value ->
            def meta = "Subchrom"
            tuple(meta, value) }
        ch_binsize_delly = channel.of(params.delly_bin_size)
            .map { value ->
            def meta = "Delly"
            tuple(meta, value) }

        ch_binsizes = ch_binsize_qdnaseq
            .mix(ch_binsize_subchrom)
            .mix(ch_binsize_delly)

        FIGENO_REPORT(
            VARIANT_PROCESS.out.circos_plot,
            VARIANT_PROCESS.out.panchr_plot,
            ch_binsizes,
            VARIANT_PROCESS.out.sv_plot,
            VARIANT_PROCESS.out.fusion_plot,
            VARIANT_PROCESS.out.targets_plot,
            VARIANT_PROCESS.out.sv_table,
            VARIANT_PROCESS.out.fusion_table,
            ch_subchrom_plot,
            ch_subchrom_focal
        )
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SECTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Collect sections from all analysis steps
    ch_sections = WGS_REPORT.out.sections
        .mix(FIGENO_REPORT.out.sections)

    ch_mode = channel.of("WGS")

     // channel id containing only meta
    ch_id = MAPPING.out.bam
        .map { meta, _bam, _bai ->
        meta }

    MIDNIGHT_REPORT(
        ch_id,
        ch_sections,
        ch_versions,
        ch_mode
    )
}
