// Basecalling subworkflows
include { BASECALL_SIMPLEX   } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'

// Core analysis subworkflows
// Core analysis subworkflows
include { MAPPING as MAPPING_HG  } from '../subworkflows/local/mapping/mapping'
include { MAPPING as MAPPING_T2T } from '../subworkflows/local/mapping/mapping'


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

// Tumor classification
include { CLASSY                               } from '../subworkflows/local/methylation_analysis/marlin.nf'
include { STURGEON                             } from '../subworkflows/local/methylation_analysis/sturgeon.nf'
include { SUBSAMPLE_TIME as SUBSAMPLE_TIME_BAM } from '../subworkflows/local/read_processing/subsample_time.nf'
include { SUBSAMPLE_TIME as SUBSAMPLE_TIME_FQ  } from '../subworkflows/local/read_processing/subsample_time.nf'

// Reporting
include { MIDNIGHT_REPORT   } from '../subworkflows/local/report/final_report.nf'
include { FIGENO_REPORT     } from '../subworkflows/local/report/variants.nf'
include { WGS_REPORT        } from '../subworkflows/local/report/wgs.nf'
include { CLASSIFIER_REPORT } from '../subworkflows/local/report/methylation.nf'

workflow WGS {

    take:
    samplesheet             // channel: samplesheet read in from --input
    demux_samplesheet       // channel : demux samplesheet read in from --demux_samplesheet
    ref                     // channel : reference for mapping, either empty if skipping mapping, or a path
    clairs_model
    basecall_model
    ch_clin_database
    bed
    targets                 // channel : list of genes with their position to represent in Figeno
    minqs                   // channel obtained dynamically from params.basecall_model
    tumor_type
    ref_t2t

    main:

    //
    // WORKFLOW: Run pipeline
    //

    ch_versions = channel.empty()

    ch_tumor_type = tumor_type
            .branch { meta, tumor ->
                leukemia: tumor == "leukemia"
                cns: tumor == "cns"
                other: tumor == "other"
            }

    if (params.skip_mapping) {

        MAPPING_HG (
            samplesheet,
            ref
        )

        ch_seqkit = MAPPING_HG.out.seqkit

        ch_in_mapt2t = samplesheet
            .join(ch_tumor_type.cns)
            .map { meta, bam, _tumor ->
            tuple(meta, bam)
            }

        MAPPING_T2T(
            ch_in_mapt2t,
            ref_t2t
        )

        STURGEON(
            MAPPING_T2T.out.bam_t2t,
            ref_t2t
        )

    } else if (params.skip_basecalling) {

        MAPPING_HG (
            samplesheet,
            ref
        )

        ch_seqkit = MAPPING_HG.out.seqkit

        ch_fastq = samplesheet

        // Downsample for running Sturgeon

        ch_in_subsample_sturgeon = ch_fastq
            .join(ch_tumor_type.cns)
            .map { meta, fastq, _tumor_type ->
                tuple(meta, fastq, 0, 1)
            }

        SUBSAMPLE_TIME_FQ(
            ch_in_subsample_sturgeon,
            'fq'
        )

        MAPPING_T2T(
            SUBSAMPLE_TIME_FQ.out.fq,
            ref_t2t
        )

        STURGEON(
            MAPPING_T2T.out.bam,
            ref_t2t
        )

    } else {

        if (params.demux) {

            BASECALL_MULTIPLEX (
                samplesheet,
                demux_samplesheet
            )

            MAPPING_HG (
                BASECALL_MULTIPLEX.out.fastq,
                ref
            )

            ch_fastq = BASECALL_MULTIPLEX.out.fastq

            ch_seqkit   = BASECALL_MULTIPLEX.out.stats_pass
            ch_versions = BASECALL_MULTIPLEX.out.versions

        } else {

            BASECALL_SIMPLEX (
                samplesheet
            )

            MAPPING_HG (
                BASECALL_SIMPLEX.out.fastq,
                ref
            )

            ch_fastq = BASECALL_SIMPLEX.out.fastq

            ch_seqkit   = BASECALL_SIMPLEX.out.stats_pass
            ch_versions = BASECALL_SIMPLEX.out.versions
        }

        // Downsample for running Sturgeon

        ch_in_subsample_sturgeon = ch_fastq
            .join(ch_tumor_type.cns)
            .map { meta, fastq, _tumor_type ->
                tuple(meta, fastq, 0, 1)
            }

        SUBSAMPLE_TIME_FQ(
            ch_in_subsample_sturgeon,
            'fq'
        )

        MAPPING_T2T(
            SUBSAMPLE_TIME_FQ.out.fq,
            ref_t2t
        )

        STURGEON(
            MAPPING_T2T.out.bam,
            ref_t2t
        )
    }

        // Downsample to 1h to run methylation classification

        ch_in_subsample = MAPPING_HG.out.bam
            .join(ch_tumor_type.leukemia)
            .map { meta, bam, index, _tumor_type ->
            tuple(meta, bam, index, 0, 1)
            }

        SUBSAMPLE_TIME_BAM(
            ch_in_subsample,
            'bam'
        )

        CLASSY(
            SUBSAMPLE_TIME_BAM.out.bam,
            ref
        )

        CLAIRS_TO_CALLING (
            MAPPING_HG.out.bam,
            ref,
            clairs_model,
            ch_clin_database,
            bed,
            vep_cache
        )

        CLAIR3_CALLING (
            MAPPING_HG.out.bam,
            ref,
            basecall_model,
            ch_clin_database,
            bed,
            vep_cache
        )

        PHASING_SOMATIC (
            MAPPING_HG.out.bam,
            ref,
            CLAIRS_TO_CALLING.out.vcf_snpeff
        )

        PHASING_GERMLINE (
            MAPPING_HG.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf_snpeff
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
            MAPPING_HG.out.bam,
            ref
        )

        SUBCHROM_CALL (
            MAPPING_HG.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf_snpeff,
            bed
        )

        ch_snp_to_process = CLAIR3_CALLING.out.vcf_vep
            .mix(CLAIRS_TO_CALLING.out.vcf_vep)

        // Filter variants to visualize :
        VARIANT_PROCESS (
            MAPPING_HG.out.bam,
            SV_CALLING.out.vcf,
            CNV_CALLING.out.qdnaseq_bed,
            CNV_CALLING.out.qdnaseq_segs,
            targets,
            CNV_CALLING.out.delly_cov,
            CNV_CALLING.out.delly_segs,
            ch_snp_to_process
        )

    ch_versions = ch_versions
        .mix(MAPPING_HG.out.versions)
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
            MAPPING_HG.out.coverage,
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
            ch_subchrom_focal,
            VARIANT_PROCESS.out.snp_table
        )

        ch_classifiers_plots = STURGEON.out.plot
            .mix(CLASSY.out.plot)

        ch_classifiers_pred = STURGEON.out.pred
            .mix(CLASSY.out.pred)

        CLASSIFIER_REPORT(
            ch_classifiers_plots,
            ch_classifiers_pred
        )
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SECTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Collect sections from all analysis steps
    ch_sections = WGS_REPORT.out.sections
        .mix(FIGENO_REPORT.out.sections)
        .mix(CLASSIFIER_REPORT.out.sections)

    ch_mode = channel.of("WGS")

     // channel id containing only meta
    ch_id = MAPPING_HG.out.bam
        .map { meta, _bam, _bai ->
        meta }

    MIDNIGHT_REPORT(
        ch_id,
        ch_sections,
        ch_versions,
        ch_mode
    )
}
