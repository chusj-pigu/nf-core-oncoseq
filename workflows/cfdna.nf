include { BASECALL_MULTIPLEX     } from '../subworkflows/local/basecalling/basecall_multiplex'
include { BASECALL_SIMPLEX       } from '../subworkflows/local/basecalling/basecall_simplex'
include { READS_FILTER           } from '../subworkflows/local/read_processing/reads_filter.nf'
include { MAPPING as MAPPING_HG  } from '../subworkflows/local/mapping/mapping'
include { MAPPING as MAPPING_T2T } from '../subworkflows/local/mapping/mapping'
include { TIDEHUNTER_CONCENSUS   } from '../subworkflows/local/read_processing/tidehunter'
include { CNV_CALLING            } from '../subworkflows/local/variant_calling/cnv_calling.nf'
include { SV_CALLING             } from  '../subworkflows/local/variant_calling/sv_calling.nf'
include { CLAIR3_CALLING         } from '../subworkflows/local/variant_calling/clair3_calling.nf'
include { VARIANT_PROCESS        } from  '../subworkflows/local/variant_calling/variant_process.nf'
include { ICHORCNA_CALLING       } from '../subworkflows/local/variant_calling/ichor_calling.nf'
include { CLASSY                 } from '../subworkflows/local/methylation_analysis/marlin.nf'
include { STURGEON               } from '../subworkflows/local/methylation_analysis/sturgeon.nf'
include { SUBCHROM_CALL          } from '../subworkflows/local/variant_calling/subchrom_call.nf'

// Reporting
include { MIDNIGHT_REPORT      } from '../subworkflows/local/report/final_report.nf'
include { CFNDA_REPORT         } from '../subworkflows/local/report/cfdna.nf'
include { CLASSIFIER_REPORT    } from '../subworkflows/local/report/methylation.nf'
include { FIGENO_REPORT        } from '../subworkflows/local/report/variants.nf'

workflow CFDNA {

    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    cfdna_samplesheet   // channel : from demux or samplesheeet
    tumor_type      // channel: to decide if marlin or sturgeon is run
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
    ref_t2t
    max_len
    minqs
    ichor_bin
    mapq_wig
    basecall_model          // channel: model for basecalling
    ch_clin_database        // channel: clinical database for variant annotation
    bed                     // channel: bed file used for adaptive sampling regions
    targets

    main:

    //
    // WORKFLOW: Run pipeline
    //

    ch_versions = channel.empty()

    if (params.demux) {

        BASECALL_MULTIPLEX (
            samplesheet,
            demux
        )

        ch_fastq = BASECALL_MULTIPLEX.out.fastq

        ch_versions = BASECALL_MULTIPLEX.out.versions

    } else if (params.skip_basecalling || params.skip_mapping) {

        ch_fastq = samplesheet

    } else {

         BASECALL_SIMPLEX (
            samplesheet
        )

        ch_fastq = BASECALL_SIMPLEX.out.fastq

        ch_versions = BASECALL_SIMPLEX.out.versions

    }

    if (params.rca) {
        TIDEHUNTER_CONCENSUS(ch_fastq)

        ch_fastq_processed = TIDEHUNTER_CONCENSUS.out.fastq

    } else {
        READS_FILTER (
            cfdna_samplesheet,
            ch_fastq,
            max_len,
            minqs
        )

        ch_fastq_processed = READS_FILTER.out.reads

        ch_tumor_type = tumor_type
            .branch { meta, tumor ->
                leukemia: tumor == "leukemia"
                cns: tumor == "cns"
                other: tumor == "other"
            }

        ch_mapping_t2t = ch_tumor_type.cns
            .join(ch_fastq_processed)
            .map { meta, _tumor, input ->
                tuple(meta, input)}

        MAPPING_HG(
            ch_fastq_processed,
            ref
        )

        MAPPING_T2T(
            ch_mapping_t2t,
            ref_t2t
        )

        // Run snp calling only for higher coverage samples
        ch_coverage = MAPPING_HG.out.coverage
            .map { meta, table ->
                def lines = table.readLines()
                def cov = lines[5].tokenize('\t')[1].toDouble()
                tuple(meta, cov)
            }
            .branch {
                meta, cov ->
                high: cov >= 10
                 return meta
            }

        ch_high_cov_bam = ch_coverage.high
            .join(MAPPING_HG.out.bam)

        CLAIR3_CALLING (
            ch_high_cov_bam,
            ref,
            basecall_model,
            ch_clin_database,
            bed
        )

        // Run Subchrom only for 15X samples

        ch_subchrom_meta = MAPPING_HG.out.coverage
            .map { meta, table ->
                def lines = table.readLines()
                def cov = lines[5].tokenize('\t')[1].toDouble()
                tuple(meta, cov)
            }
            .branch { meta, cov ->
                higher: cov >= 15
                    return meta
            }

        ch_subcrhom_bam = ch_subchrom_meta.higher
            .join(MAPPING_HG.out.bam)

        SUBCHROM_CALL (
            ch_subcrhom_bam,
            ref,
            CLAIR3_CALLING.out.vcf,
            bed
        )

        SV_CALLING(
            MAPPING_HG.out.bam,
            ref
        )

        // Run Marlin on leukemia samples only:
        ch_marlin_bam = MAPPING_HG.out.bam
            .join(ch_tumor_type.leukemia)
            .map { meta, bam, bai, _tumor ->
                tuple(meta, bam, bai)}

        CLASSY(
            ch_marlin_bam,
            ref
        )

        STURGEON(
            MAPPING_T2T.out.bam
        )

        CNV_CALLING (
            MAPPING_HG.out.bam,
            ref
        )

        ICHORCNA_CALLING (
            MAPPING_HG.out.bam,
            cfdna_samplesheet,
            ichor_bin,
            mapq_wig
        )

        // Sub reports

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

        VARIANT_PROCESS (
            MAPPING_HG.out.bam,
            SV_CALLING.out.vcf,
            CNV_CALLING.out.qdnaseq_bed,
            CNV_CALLING.out.qdnaseq_segs,
            targets,
            CNV_CALLING.out.delly_cov,
            CNV_CALLING.out.delly_segs
        )

        ch_subchrom_plot = SUBCHROM_CALL.out.subchrom_plot_wgs
        ch_subchrom_focal = SUBCHROM_CALL.out.subchrom_gene_plot_wgs

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
            COLLECT VERSIONS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */


        ch_versions = ch_versions
            .mix(MAPPING_HG.out.versions)
            .mix(MAPPING_T2T.out.versions)
            .mix(CLASSY.out.versions)
            .mix(STURGEON.out.versions)
            .mix(CNV_CALLING.out.versions)
            .mix(ICHORCNA_CALLING.out.versions)

        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            COMPILE SECTIONS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        ch_classifiers_plots = STURGEON.out.plot
            .mix(CLASSY.out.plot)

        ch_classifiers_pred = STURGEON.out.pred
            .mix(CLASSY.out.pred)

        CFNDA_REPORT(
            READS_FILTER.out.stats,
            MAPPING_HG.out.coverage,
            ICHORCNA_CALLING.out.ichorcna_plot
        )

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
        ch_sections = CFNDA_REPORT.out.sections
            .mix(CLASSIFIER_REPORT.out.sections)
            .mix(FIGENO_REPORT.out.sections)


        ch_mode = channel.of("cfDNA")

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
}
