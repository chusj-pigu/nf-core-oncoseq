include { BASECALL_MULTIPLEX     } from '../subworkflows/local/basecalling/basecall_multiplex'
include { BASECALL_SIMPLEX       } from '../subworkflows/local/basecalling/basecall_simplex'
include { READS_FILTER           } from '../subworkflows/local/read_processing/reads_filter.nf'
include { MAPPING                } from '../subworkflows/local/mapping/mapping'
include { TIDEHUNTER_CONCENSUS   } from '../subworkflows/local/read_processing/tidehunter'
include { CNV_CALLING            } from '../subworkflows/local/variant_calling/cnv_calling.nf'
include { SV_CALLING             } from  '../subworkflows/local/variant_calling/sv_calling.nf'
include { CLAIRS_TO_CALLING      } from '../subworkflows/local/variant_calling/clairs_to_calling.nf'
include { CLAIR3_CALLING         } from '../subworkflows/local/variant_calling/clair3_calling.nf'
include { VARIANT_PROCESS        } from  '../subworkflows/local/variant_calling/variant_process.nf'
include { ICHORCNA_CALLING       } from '../subworkflows/local/variant_calling/ichor_calling.nf'
include { CLASSY                 } from '../subworkflows/local/methylation_analysis/classy.nf'
include { SUBSAMPLE_TIME as SUBSAMPLE_TIME_BAM } from '../subworkflows/local/read_processing/subsample_time.nf'
include { SAMTOOLS_COUNT_READS                 } from '../modules/local/samtools/main.nf'

// Reporting
include { MIDNIGHT_REPORT      } from '../subworkflows/local/report/final_report.nf'
include { CFDNA_REPORT         } from '../subworkflows/local/report/cfdna.nf'
include { CLASSIFIER_REPORT    } from '../subworkflows/local/report/methylation.nf'
include { FIGENO_REPORT        } from '../subworkflows/local/report/variants.nf'

workflow CFDNA {

    take:
    samplesheet // channel: samplesheet read in from --input
    demux       // channel: demux_samplesheet read in from --demux_samplesheet
    cfdna_samplesheet   // channel : from demux or samplesheeet
    ref         // channel : reference for mapping, either empty if skipping mapping, or a path
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

    ch_versions = channel.empty()

    if (params.demux) {

        BASECALL_MULTIPLEX (
            samplesheet,
            demux
        )

        ch_versions = BASECALL_MULTIPLEX.out.versions

        READS_FILTER (
            cfdna_samplesheet,
            BASECALL_MULTIPLEX.out.fastq,
            max_len,
            minqs
        )

        MAPPING(
            READS_FILTER.out.reads,
            ref
        )

        ch_to_classify = MAPPING.out.bam
            .map { meta, _bam, _bai -> meta }

    } else if (params.skip_basecalling || params.skip_mapping) {

        READS_FILTER (
            cfdna_samplesheet,
            samplesheet,
            max_len,
            minqs
        )

        ch_fastq_processed = READS_FILTER.out.reads

        MAPPING(
            ch_fastq_processed,
            ref
        )

        SAMTOOLS_COUNT_READS(MAPPING.out.bam.map { meta, bam, _bai -> tuple(meta, bam)})

        ch_bam_methylation_counts = SAMTOOLS_COUNT_READS.out.txt
            .map { meta, txt ->
                def count = txt.text.trim().toInteger()
                tuple(meta, count)
            }
            .branch { meta, count ->
                pos:  count > 0
                    return meta
                none: true
            }

        ch_bam_methylation_counts.none
            .subscribe { meta, count ->
                log.warn "Tumor classification will be skipped for ${meta.id} -- no methylation tags found in bam"
            }
        ch_to_classify = ch_bam_methylation_counts.pos

    } else {

         BASECALL_SIMPLEX (
            samplesheet
        )

        ch_fastq = BASECALL_SIMPLEX.out.fastq

        ch_versions = BASECALL_SIMPLEX.out.versions

        READS_FILTER (
            cfdna_samplesheet,
            ch_fastq,
            max_len,
            minqs
        )

        ch_fastq_processed = READS_FILTER.out.reads

        MAPPING(
            ch_fastq_processed,
            ref
        )

        ch_to_classify = MAPPING.out.bam
            .map { meta, _bam, _bai -> meta }

    }

    ch_vcf_subchrom = channel.empty()

    if (params.rca) {
        TIDEHUNTER_CONCENSUS(ch_fastq)

        ch_fastq_processed = TIDEHUNTER_CONCENSUS.out.fastq

    } else {

        // Run snp calling only for higher coverage samples
        ch_coverage = MAPPING.out.coverage
            .map { meta, table ->
                def lines = table.readLines()
                def cov = lines[5].tokenize('\t')[1].toDouble()
                tuple(meta, cov)
            }
            .branch {
                meta, cov ->
                high: cov >= 10
                 return meta
                low: cov < 10
                 return meta
            }

        ch_high_cov_bam = ch_coverage.high
            .join(MAPPING.out.bam)

        CLAIR3_CALLING (
            ch_high_cov_bam,
            ref,
            basecall_model,
            bed,
            vep_cache
        )

        CLAIRS_TO_CALLING (
            ch_high_cov_bam,
            ref,
            clairs_model,
            bed,
            vep_cache
        )

        ch_vcf_subchrom = ch_vcf_subchrom
            .mix(CLAIR3_CALLING.out.vcf_snpeff)

        SV_CALLING(
            MAPPING.out.bam,
            ref
        )

        if (params.m_bases || params.skip_basecalling || params.skip_mapping) {
            ch_in_subsample = ch_high_cov_bam
                .join(ch_to_classify)
                .map { meta, bam, index ->
                tuple(meta, bam, index, 0, 8)}          // Subsample to first 8h of sequencing to avoid slowing down classification

            SUBSAMPLE_TIME_BAM(
                ch_in_subsample
            )

            ch_in_classy = SUBSAMPLE_TIME_BAM.out.bam
                .mix(ch_coverage.low.join(MAPPING.out.bam))
                .join(ch_to_classify)

            CLASSY(
                ch_in_classy,
                ref
            )

            CLASSIFIER_REPORT(
                CLASSY.out.plot,
                CLASSY.out.pred
            )

            ch_classy_section = CLASSIFIER_REPORT.out.sections

            ch_versions = ch_versions
                .mix(SUBSAMPLE_TIME_BAM.out.versions)
                .mix(CLASSY.out.versions)

        } else {
            ch_classy_section = channel.empty()

            ch_versions = ch_versions
        }

        CNV_CALLING (
            MAPPING.out.bam,
            ch_vcf_subchrom,
            ref
        )

        ICHORCNA_CALLING (
            MAPPING.out.bam,
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

        ch_snp_to_process = CLAIR3_CALLING.out.vcf_vep
            .mix(CLAIRS_TO_CALLING.out.vcf_vep)

        VARIANT_PROCESS (
            MAPPING.out.bam,
            SV_CALLING.out.vcf,
            CNV_CALLING.out.qdnaseq_bed,
            CNV_CALLING.out.qdnaseq_segs,
            targets,
            CNV_CALLING.out.delly_cov,
            CNV_CALLING.out.delly_segs,
            ch_snp_to_process
        )

        ch_subchrom_plot = CNV_CALLING.out.subchrom_plot_wgs
        ch_subchrom_focal = CNV_CALLING.out.subchrom_gene_plot_wgs

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

        CFDNA_REPORT(
            READS_FILTER.out.stats,
            MAPPING.out.coverage,
            ICHORCNA_CALLING.out.ichorcna_plot
        )

         /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            COLLECT VERSIONS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */


        ch_versions = ch_versions
            .mix(MAPPING.out.versions)
            .mix(CNV_CALLING.out.versions)
            .mix(ICHORCNA_CALLING.out.versions)
            .mix(CLAIR3_CALLING.out.versions)
            .mix(CLAIRS_TO_CALLING.out.versions)
            .mix(READS_FILTER.out.versions)
            .mix(SV_CALLING.out.versions)
            .mix(VARIANT_PROCESS.out.versions)
            .mix(FIGENO_REPORT.out.versions)
            .mix(CFDNA_REPORT.out.versions)

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            COLLECT SECTIONS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // Collect sections from all analysis steps
        ch_sections = CFDNA_REPORT.out.sections
            .mix(ch_classy_section)
            .mix(FIGENO_REPORT.out.sections)

         // channel id containing only meta
        ch_id = MAPPING.out.bam
            .map { meta, _bam, _bai ->
            meta }

        ch_title = ch_id
            .map { meta ->
            tuple(meta, "OncoSeq cfDNA Report — ${meta.id}")}

        MIDNIGHT_REPORT(
            ch_id,
            ch_sections,
            ch_versions,
            ch_title
        )
    }
}
