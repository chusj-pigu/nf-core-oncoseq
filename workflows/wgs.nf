// Basecalling subworkflows
include { BASECALL_SIMPLEX   } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'

// Core analysis subworkflows
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

// Tumor classification
include { SAMTOOLS_COUNT_READS                 } from '../modules/local/samtools/main.nf'
include { CLASSY                               } from '../subworkflows/local/methylation_analysis/classy.nf'
include { SUBSAMPLE_TIME as SUBSAMPLE_TIME_BAM } from '../subworkflows/local/read_processing/subsample_time.nf'

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
    bed
    targets                 // channel : list of genes with their position to represent in Figeno
    minqs                   // channel obtained dynamically from params.basecall_model
    vep_cache

    main:

    //
    // WORKFLOW: Run pipeline
    //

    ch_versions = channel.empty()

    if (params.skip_mapping || params.skip_basecalling) {

        MAPPING (
            samplesheet,
            ref
        )

        ch_seqkit = MAPPING.out.seqkit

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
            .join(MAPPING.out.bam)

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

            ch_fastq = BASECALL_MULTIPLEX.out.fastq

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

            ch_fastq = BASECALL_SIMPLEX.out.fastq

            ch_seqkit   = BASECALL_SIMPLEX.out.stats_pass
            ch_versions = BASECALL_SIMPLEX.out.versions
        }
        ch_to_classify = MAPPING.out.bam
    }

    if (params.m_bases || params.skip_basecalling || params.skip_mapping) {
        // Downsample to 1h to run methylation classification

        ch_in_subsample = ch_to_classify
            .map { meta, bam, index ->
            tuple(meta, bam, index, 0, 1)
            }

        SUBSAMPLE_TIME_BAM(
            ch_in_subsample
        )

        ch_in_classy = SUBSAMPLE_TIME_BAM.out.bam

        CLASSY(
            ch_in_classy,
            ref
        )

        CLASSIFIER_REPORT(
            CLASSY.out.plot,
            CLASSY.out.pred
        )

        ch_classy_section = CLASSIFIER_REPORT.out.sections
    } else {
        ch_classy_section = channel.empty()
    }

        CLAIRS_TO_CALLING (
            MAPPING.out.bam,
            ref,
            clairs_model,
            bed,
            vep_cache
        )

        CLAIR3_CALLING (
            MAPPING.out.bam,
            ref,
            basecall_model,
            bed,
            vep_cache
        )

        PHASING_SOMATIC (
            MAPPING.out.bam,
            ref,
            CLAIRS_TO_CALLING.out.vcf_snpeff
        )

        PHASING_GERMLINE (
            MAPPING.out.bam,
            ref,
            CLAIR3_CALLING.out.vcf_snpeff
        )

        SV_CALLING (
            PHASING_GERMLINE.out.haptag_bam
                .map { meta, bamfile, bai ->
                    // Restore original sample ID for output naming
                    def meta_restore = modifyMetaId(meta, 'replace', '_somatic_snp_snpeff_phased', '', '')
                    meta_restore = modifyMetaId(meta_restore, 'replace', '_germline_snp_snpeff_phased', '', '')
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
            CLAIR3_CALLING.out.vcf_snpeff,
            bed
        )

        ch_snp_to_process = CLAIR3_CALLING.out.vcf_vep
            .mix(CLAIRS_TO_CALLING.out.vcf_vep)

        // Filter variants to visualize :
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
            ch_subchrom_focal,
            VARIANT_PROCESS.out.snp_table
        )
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SECTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Collect sections from all analysis steps
    ch_sections = WGS_REPORT.out.sections
        .mix(FIGENO_REPORT.out.sections)
        .mix(ch_classy_section)

    ch_id = MAPPING.out.bam
        .map { meta, _bam, _bai ->
        meta }

    ch_title = ch_id
        .map { meta ->
        tuple(meta, "OncoSeq WGS Report — ${meta.id}")}

    MIDNIGHT_REPORT(
        ch_id,
        ch_sections,
        ch_versions,
        ch_title
    )
}
