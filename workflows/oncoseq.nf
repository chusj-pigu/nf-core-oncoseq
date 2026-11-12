//
// OncoSeq analysis: one workflow, three modes
//
//   --cfdna     -> cfDNA mode
//   --realtime  -> local realtime mode   (ignored if --cfdna; cfDNA handles realtime titles itself)
//   (otherwise) -> adaptive WGS mode
//

// Basecalling
include { BASECALL_SIMPLEX   } from '../subworkflows/local/basecalling/basecall_simplex'
include { BASECALL_MULTIPLEX } from '../subworkflows/local/basecalling/basecall_multiplex'

// Read processing
include { READS_FILTER   } from '../subworkflows/local/read_processing/reads_filter.nf'
include { SUBSAMPLE_TIME } from '../subworkflows/local/read_processing/subsample_time.nf'

// Mapping
include { MAPPING } from '../subworkflows/local/mapping/mapping'

// Variant calling
include { SNP_CALLING                                     } from '../subworkflows/local/variant_calling/snp_calling.nf'
include { VCF_ANNOTATION as SNP_ANNOTATION                } from '../subworkflows/local/variant_calling/annotation.nf'
include { VCF_ANNOTATION as VAR_ANNOTATION                } from '../subworkflows/local/variant_calling/annotation.nf'
include { BCFTOOLS_CONCAT as BCFTOOLS_SNP_PHASED_CONCAT   } from '../modules/nf-core/bcftools/concat'
include { BCFTOOLS_CONCAT as BCFTOOLS_SNP_UNPHASED_CONCAT } from '../modules/nf-core/bcftools/concat'
include { PHASING_VARIANTS                                } from '../subworkflows/local/variant_calling/phasing.nf'
include { PHARMCAT                                        } from '../subworkflows/local/variant_calling/pharmcat.nf'
include { SV_CALLING                                      } from '../subworkflows/local/variant_calling/sv_calling.nf'
include { CNV_CALLING                                     } from '../subworkflows/local/variant_calling/cnv_calling.nf'
include { NASVAR                                          } from '../subworkflows/local/variant_calling/nasvar_pipeline.nf'
include { ICHORCNA_CALLING                                } from '../subworkflows/local/variant_calling/ichor_calling.nf'
include { VARIANT_PROCESS                                 } from '../subworkflows/local/variant_calling/variant_process.nf'

// Adaptive-specific
include { COVERAGE_SEPARATE } from '../subworkflows/local/adaptive_specific/coverage_separate'

// Tumor classification
include { CLASSY               } from '../subworkflows/local/methylation_analysis/classy.nf'
include { SAMTOOLS_COUNT_READS } from '../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TIME_FASTQ } from '../modules/local/samtools/main.nf'

// Time series
include { SPLIT_BAMS_TIME } from '../subworkflows/local/time_series_evaluation/split_bams.nf'
include { CRAMINO         } from '../modules/nf-core/cramino/main'

// Realtime
include { ONTIME_TIME_RANGE } from '../modules/local/ontime/main.nf'

// Reporting
include { MIDNIGHT_REPORT   } from '../subworkflows/local/report/final_report.nf'
include { VARIANT_REPORT    } from '../subworkflows/local/report/variants.nf'
include { ADAPTIVE_REPORT   } from '../subworkflows/local/report/adaptive.nf'
include { CLASSIFIER_REPORT } from '../subworkflows/local/report/methylation.nf'
include { CFDNA_REPORT      } from '../subworkflows/local/report/cfdna.nf'

// Utilities
include { modifyMetaId           } from '../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

// Useful functions to handle time parsing
def parseToInstant(str) {
    try {
        return java.time.Instant.parse(str)
    } catch (e) {}
    try {
        return java.time.OffsetDateTime.parse(str).toInstant()
    } catch (e) {}
    try {
        def fmt = java.time.format.DateTimeFormatter
            .ofPattern("yyyy-MM-dd'T'HH:mm:ss.SSSz")
            .withZone(java.time.ZoneId.of("UTC"))
        return java.time.ZonedDateTime.parse(str, fmt).toInstant()
    } catch (e) {}
    throw new RuntimeException("Could not parse timestamp: ${str}")
}


/**
 * Generate a time-point sample ID.
 */
def timed_id(meta, t) {
    return "${meta.id}_0h_${t}h"
}

/**
 * Add a FULL version of each channel item.
 * Preserves all metadata fields and the remaining tuple elements.
 */
def add_full(ch) {
    return ch.map { item ->
        def l = item.toList()
        l[0] = l[0] + [id: "${l[0].id}_FULL"]
        l
    }
}

/**
 * Include FULL samples when params.include_full is enabled.
 * Otherwise, return only the time-point samples.
 */
def with_full(ch, timed) {
    return params.include_full ? add_full(ch).mix(timed) : timed
}


workflow ONCOSEQ {

    take:
    samplesheet             // channel: samplesheet read in from --input
    demux_samplesheet       // channel: demux samplesheet read in from --demux_samplesheet
    tumor_type              // channel: tumor type read in from samplesheet
    cfdna_samplesheet       // channel: cfDNA samplesheet (cfDNA mode only)
    ref                     // channel: reference fasta
    ref_fai                 // channel : reference fasta file index
    genome                  // variable: standardized genome ID
    ichor_bin               // cfDNA mode only
    mapq_wig                // cfDNA mode only
    basecall_model          // channel: model for basecalling
    bed                     // channel: bed file used for adaptive sampling regions
    bed_nopad                // channel: bed file without padding
    targets                 // channel: list of genes with their position to represent in Figeno
    vep_cache

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MODE FLAGS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    def is_cfdna       = params.cfdna ? true : false
    def is_realtime    = (!is_cfdna && params.realtime) ? true : false
    def is_full_wgs    = (!is_cfdna && !is_realtime)
    def realtime_hours = is_realtime ? params.realtime?.toInteger() : 0
    def time_series    = (params.time_series && !is_realtime) ? true : false
    def check_tags     = (params.skip_basecalling || params.skip_mapping) ? true : false
    def do_classify    = ((params.m_bases || check_tags) && genome != "hs1") ? true : false
    def run_snp        = (!is_realtime || realtime_hours >= 6)

    def title_prefix = is_cfdna     ? 'OncoSeq cfDNA Report'
                     : is_realtime  ? 'OncoSeq Adaptive Sampling Report'
                     :                'OncoSeq Workflow Report'

    ch_versions = channel.empty()

    // Transform ref channel into value channel so it doesn't get consumed only once
    ch_ref = ref.first()
    ch_ref_fai = ref_fai.first()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        BASECALLING -> (READS FILTER) -> MAPPING -> (READS FILTER)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.skip_basecalling || params.skip_mapping) {

        ch_reads   = samplesheet

    } else if (params.demux) {

        BASECALL_MULTIPLEX (
            samplesheet,
            demux_samplesheet
        )

        ch_reads   = BASECALL_MULTIPLEX.out.fastq

        ch_versions = ch_versions
            .mix(BASECALL_MULTIPLEX.out.versions)

    } else {

        BASECALL_SIMPLEX (
            samplesheet
        )

        ch_reads   = BASECALL_SIMPLEX.out.fastq

        ch_versions = ch_versions
            .mix(BASECALL_SIMPLEX.out.versions)
    }

    // cfDNA: filter reads before mapping (unless the input is already mapped)
    if (is_cfdna && !params.skip_mapping) {

        READS_FILTER (
            cfdna_samplesheet,
            ch_reads
        )

        ch_to_map = READS_FILTER.out.reads

    } else {
        ch_to_map = ch_reads
    }

    MAPPING(
        ch_to_map,
        ch_ref,
        ch_ref_fai
    )

    // cfDNA with mapped input: filter the BAM after mapping
    if (is_cfdna && params.skip_mapping) {

        READS_FILTER (
            cfdna_samplesheet,
            MAPPING.out.bam
        )

        ch_bam_res = READS_FILTER.out.reads

    } else {
        ch_bam_res = MAPPING.out.bam
    }

    ch_versions = ch_versions
        .mix(MAPPING.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        REALTIME: sequencing time range (used for report title)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.realtime) {

        if (params.skip_mapping) {
            // We have to convert bam to fastq when bam input is from MinKnow :

            ch_bam = MAPPING.out.bam
            .map { meta, bam, _bai ->
                tuple(meta,bam) }

            SAMTOOLS_TIME_FASTQ(ch_bam)

            ONTIME_TIME_RANGE(SAMTOOLS_TIME_FASTQ.out.fq)
        } else {
            ONTIME_TIME_RANGE(ch_reads)
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        METHYLATION TAG CHECK: which samples can be classified
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (do_classify) {

        if (check_tags) {

            SAMTOOLS_COUNT_READS(ch_bam_res.map { meta, bam, _bai -> tuple(meta, bam)})

            ch_bam_methylation_counts = SAMTOOLS_COUNT_READS.out.txt
                .map { meta, txt ->
                    def count = txt.text.trim().toInteger()
                    tuple(meta, count)
                }
                .branch { meta, count ->
                    pos:  count > 0
                        return meta
                    none: true
                        return meta
                }

            ch_bam_methylation_counts.none
                .subscribe { meta ->
                    log.warn "Tumor classification will be skipped for ${meta.id} -- no methylation tags found in bam"
                }

            ch_to_classify = ch_bam_methylation_counts.pos

        } else {
            // m_bases: all samples carry methylation tags
            ch_to_classify = ch_bam_res
                .map { meta, _bam, _bai -> meta }
        }
    }

    /*
     ~ *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     TIME SERIES SPLIT (or full BAM)
     Samples are renamed <id>_0h_<N>h per time point, plus <id>_FULL with --include_full.
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (time_series) {

        SPLIT_BAMS_TIME(
            ch_bam_res,
            bed
        )

        // Use time series outputs for downstream variant calling
        ch_bam_for_calling = SPLIT_BAMS_TIME.out.bam
        .map { meta, bam, bai ->
            tuple(id:meta.id, bam, bai)}

            ch_bed_for_calling = SPLIT_BAMS_TIME.out.bed
            .map { meta,bedfile,padding,low_fidelity ->
                tuple(id:meta.id,bedfile,padding,low_fidelity)}

                // One item per time point
                ch_time_points = channel.of(params.time_points)
                .flatMap { times -> times.tokenize(',').collect { t -> t.replace('"', '') } }

                ch_tumor_type = with_full(
                    tumor_type,
                    tumor_type
                    .combine(ch_time_points)
                    .map { meta, tumor, t -> tuple(id:timed_id(meta, t), tumor) }
                )

                bed_nopad_final = with_full(
                    bed_nopad,
                    bed_nopad
                    .combine(ch_time_points)
                    .map { meta, bedfile, t -> tuple(id:timed_id(meta, t), bedfile) }
                )

                ch_versions = ch_versions
                .mix(SPLIT_BAMS_TIME.out.versions)

    } else {
        // Standard mode: Use the full BAM directly for variant calling
        ch_bam_for_calling = ch_bam_res
        ch_bed_for_calling = bed
        ch_tumor_type      = tumor_type
        bed_nopad_final = bed_nopad
    }

    /*
     ~ *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     cfDNA-SPECIFIC SAMPLE CHANNELS
     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */

    if (is_cfdna) {

        if (time_series) {

            ch_cfdna_ss = cfdna_samplesheet
                .combine(ch_time_points)
                .map { meta, purity, filter, t -> tuple(id:timed_id(meta, t), purity, filter) }

                CRAMINO(ch_bam_for_calling.filter{ meta, _bam, _bai -> !meta.id.contains('FULL') })

                ch_cfdna_full     = with_full(cfdna_samplesheet, ch_cfdna_ss)
                ch_coverage_final = with_full(MAPPING.out.coverage, CRAMINO.out.stats)
                ch_seqkit_stats   = with_full(READS_FILTER.out.stats, channel.empty())

                ch_bam_to_process = params.include_full
                ? ch_bam_for_calling.filter{ meta, _bam, _bai -> meta.id.contains('FULL') }
                : ch_bam_for_calling

        } else {
            ch_cfdna_full     = cfdna_samplesheet
            ch_coverage_final = MAPPING.out.coverage
            ch_seqkit_stats   = READS_FILTER.out.stats
            ch_bam_to_process = ch_bam_res
        }

        ch_cov_final = ch_coverage_final
            .map { meta, table ->
                def lines = table.readLines()
                def cov = lines[5].tokenize('\t')[1].toDouble()
                tuple(meta, cov)
            }

    } else {

        ch_bam_to_process = ch_bam_for_calling
        ch_cfdna_full     = channel.empty()
        ch_coverage_final = channel.empty()
        ch_seqkit_stats   = channel.empty()
        ch_cov_final      = channel.empty()
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COVERAGE SEPARATION (target vs background)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    COVERAGE_SEPARATE(
        ch_bam_for_calling,
        ch_bed_for_calling,
        bed_nopad_final,
        ch_ref_fai
    )

    ADAPTIVE_REPORT(
        COVERAGE_SEPARATE.out.coverage_tbl,
        COVERAGE_SEPARATE.out.coverage_plot
    )

    ch_versions = ch_versions
        .mix(COVERAGE_SEPARATE.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TUMOR CLASSIFICATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (do_classify) {

        if (time_series) {

            // Classify the split BAMs of samples that carry methylation tags.
            // Split sample ids look like <id>_0h_<N>h or <id>_FULL; adaptive mode only uses the first hour.
            ch_classy_in = ch_bam_for_calling
                .filter { meta, _bam, _bai ->
                    is_cfdna || meta.id.endsWith('_0h_1h') }
                .map { meta, bam, bai ->
                    tuple(meta.id.replaceFirst(/_(FULL|0h_[^_]+h)$/, ''), meta, bam, bai) }
                .join(ch_to_classify.map { meta -> tuple(meta.id, meta) })
                .map { _key, meta, bam, bai, _orig_meta ->
                    tuple(meta, bam, bai) }

        } else if (is_full_wgs) {

            // Downsample to 1h (when coverage allows) to run methylation classification
            ch_classify_bams = ch_to_classify
                .join(MAPPING.out.bam)

            ch_adaptive_cov = COVERAGE_SEPARATE.out.coverage_tbl
                .map { meta, table ->
                    def lines = table.readLines()
                    def cov = lines[1].tokenize(',')[1].toDouble()
                    tuple(meta, cov)
                }
                .branch { meta, cov ->
                    high: cov >= 10
                        return meta
                    low: true
                        return meta
                }

            ch_in_subsample_classy = ch_classify_bams
                .map { meta, bam, index ->
                    tuple(meta, bam, index, 0, 1)
                }
                .join(ch_adaptive_cov.high)

            SUBSAMPLE_TIME(
                ch_in_subsample_classy
            )

            ch_versions = ch_versions
                .mix(SUBSAMPLE_TIME.out.versions)

            ch_classy_in = SUBSAMPLE_TIME.out.bam
                .mix(ch_classify_bams.join(ch_adaptive_cov.low))

        } else {
            ch_classy_in = ch_to_classify
                .join(MAPPING.out.bam)
        }

        CLASSY(
            ch_classy_in,
            ch_ref,
            ch_tumor_type
        )

        CLASSIFIER_REPORT(
            CLASSY.out.plot,
            CLASSY.out.pred,
            ch_tumor_type
        )

        ch_classy_section = CLASSIFIER_REPORT.out.sections

        ch_versions = ch_versions
            .mix(CLASSY.out.versions)

    } else {
        ch_classy_section = channel.empty()
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SNP CALLING
        (cfDNA: high-coverage samples only; realtime: only from 6h on)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (run_snp) {

        if (is_cfdna) {
            ch_coverage_branched = ch_cov_final
                .branch { meta, cov ->
                    high: cov >= 10
                        return meta
                    low: true
                        return meta
                    }

            ch_snp_in = ch_coverage_branched.high
                .join(ch_bam_to_process)
        } else {
            ch_snp_in = ch_bam_to_process
        }

        SNP_CALLING(
            ch_snp_in,
            ch_ref,
            ch_ref_fai,
            basecall_model,
            bed_nopad_final
        )

        ch_versions = ch_versions
            .mix(SNP_CALLING.out.versions)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SV CALLING
          full workflow : annotate SNPs -> phase -> SV calling on phased BAM
          cfDNA and/or relatime : unphased SV calling (no SNP vcf)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (is_full_wgs) {

        SNP_ANNOTATION(
            SNP_CALLING.out.vcf,
            ch_ref,
            vep_cache,
            genome
        )

        ch_vcf_to_concat = SNP_ANNOTATION.out.vcf_vep
            .map { meta, vcf_file, tbi ->
                def merged_id = meta.id.replaceAll(/_(clairsto|deepsomatic|clair3)_vep$/, '')
                tuple(id: merged_id, vcf_file, tbi)
            }
            .groupTuple()

        BCFTOOLS_SNP_PHASED_CONCAT(
            ch_vcf_to_concat
        )

        ch_snp_vcf_concat = BCFTOOLS_SNP_PHASED_CONCAT.out.vcf
            .join(BCFTOOLS_SNP_PHASED_CONCAT.out.index)

        // Phase somatic variants
        PHASING_VARIANTS (
            ch_bam_to_process,
            ch_ref,
            ch_ref_fai,
            ch_snp_vcf_concat
        )

        ch_versions = ch_versions
            .mix(PHASING_VARIANTS.out.versions)

        // Restore original sample ID for output naming
        ch_sv_bam = PHASING_VARIANTS.out.haptag_bam
            .map { meta, bamfile, bai ->
                def meta_restore = modifyMetaId(meta, 'replace', '_snp_phased', '', '')
                tuple(meta_restore, bamfile, bai)
            }

        ch_snp_vcf = PHASING_VARIANTS.out.phased_vcf
            .map { meta, vcf ->
                def meta_restore = modifyMetaId(meta, 'replace', '_snp_phased', '', '')
                tuple(meta_restore, vcf)
            }

    } else {
        ch_sv_bam = ch_bam_to_process
        ch_snp_vcf = channel.empty()    // placeholder: no phased SNP vcf
    }

    SV_CALLING(
        ch_sv_bam,
        ch_ref,
        genome,
        ch_snp_vcf
    )

    ch_versions = ch_versions
        .mix(SV_CALLING.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ANNOTATION OF VARIANTS TO VISUALIZE
        (full workflow annotated SNPs above, before phasing; the other modes annotate SV + SNP vcfs here without phasing)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_sv_to_annotate = SV_CALLING.out.vcf_to_annotate
            .filter{ meta, _vcf, _index -> !meta.id.contains("stellerator") }

    ch_vcf_to_annotate = run_snp && !is_full_wgs ? ch_sv_to_annotate.mix(SNP_CALLING.out.vcf) : ch_sv_to_annotate

    VAR_ANNOTATION(
        ch_vcf_to_annotate,
        ch_ref,
        vep_cache,
        genome
    )

    ch_vcf_annotated = VAR_ANNOTATION.out.vcf_filt
            .mix(SV_CALLING.out.vcf_to_annotate.filter{ meta, _vcf, _index -> meta.id.contains("stellerator") })

    ch_vcf_to_process =  is_full_wgs ? ch_vcf_annotated.mix(SNP_ANNOTATION.out.vcf_filt) : ch_vcf_annotated

    ch_vcf_to_concat = VAR_ANNOTATION.out.vcf_vep
        .filter { meta, _vcf_file, _tbi -> meta.id.contains('snp')}
        .map { meta, vcf_file, tbi ->
            def merged_id = meta.id.replaceAll(/_(clairsto|deepsomatic|clair3)/, '')
            tuple(id: merged_id, vcf_file, tbi)
        }
        .groupTuple()

    BCFTOOLS_SNP_UNPHASED_CONCAT(
        ch_vcf_to_concat
    )

    if (is_full_wgs) {
        // Run pharmcat on phased vcf in full workflow
        ch_snp_to_pharmcat = PHASING_VARIANTS.out.phased_vcf
            .join(PHASING_VARIANTS.out.phased_tbi)
    } else {
        ch_snp_to_pharmcat = BCFTOOLS_SNP_UNPHASED_CONCAT.out.vcf
            .join(BCFTOOLS_SNP_UNPHASED_CONCAT.out.index)
    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PHARMCAT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (genome == "hg38") {
        PHARMCAT(
            ch_snp_to_pharmcat,
            ref,
            ch_ref_fai
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CNV / NASVAR / ICHORCNA
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    CNV_CALLING(
        ch_bam_to_process,
        ch_ref,
        genome
    )

    if (genome != "hg19") {
        NASVAR(
            ch_bam_to_process,
            ch_bed_for_calling,
            bed_nopad_final,
            ch_ref,
            ch_ref_fai
        )

        ch_kar_baf_plot = NASVAR.out.baf_plot
        ch_kar_out = NASVAR.out.json
    } else {
        ch_kar_baf_plot = channel.empty()
        ch_kar_out = channel.empty()
    }

    ch_versions = ch_versions
        .mix(CNV_CALLING.out.versions)

    if (is_cfdna) {

        ICHORCNA_CALLING (
            ch_bam_for_calling,
            ch_ref,
            ch_cfdna_full,
            ichor_bin,
            mapq_wig
        )

        ch_versions = ch_versions
            .mix(ICHORCNA_CALLING.out.versions)
            .mix(READS_FILTER.out.versions)
    }

    // Filter variants to visualize :
    VARIANT_PROCESS (
        ch_bam_to_process,
        bed_nopad_final,
        ch_vcf_to_process,
        CNV_CALLING.out.qdnaseq_bed,
        CNV_CALLING.out.qdnaseq_segs,
        targets,
        CNV_CALLING.out.delly_cov,
        CNV_CALLING.out.delly_segs
    )

    ch_versions = ch_versions
        .mix(VARIANT_PROCESS.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SUB REPORTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_binsize_qdnaseq = channel.of(params.qdnaseq_binsize)
        .map { value ->
        def meta = "qDNAseq"
        tuple(meta, value) }
    ch_binsize_delly = channel.of(params.delly_bin_size)
        .map { value ->
        def meta = "Delly"
        tuple(meta, value) }

    ch_binsizes = ch_binsize_qdnaseq.mix(ch_binsize_delly)

    VARIANT_REPORT(
        VARIANT_PROCESS.out.circos_plot,
        VARIANT_PROCESS.out.panchr_plot,
        ch_kar_baf_plot,
        ch_kar_out,
        ch_binsizes,
        VARIANT_PROCESS.out.sv_plot,
        VARIANT_PROCESS.out.fusion_plot,
        VARIANT_PROCESS.out.targets_plot,
        VARIANT_PROCESS.out.sv_table,
        VARIANT_PROCESS.out.fusion_table,
        VARIANT_PROCESS.out.snp_table
    )

    if (is_cfdna) {

        CFDNA_REPORT(
            ch_cfdna_full,
            ch_seqkit_stats,
            ch_coverage_final,
            ICHORCNA_CALLING.out.ichorcna_plot
        )

        ch_versions = ch_versions
            .mix(VARIANT_REPORT.out.versions)
            .mix(CFDNA_REPORT.out.versions)

        ch_cfdna_section = CFDNA_REPORT.out.sections

    } else {
        ch_cfdna_section = channel.empty()
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT SECTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_sections = ch_classy_section
        .mix(ADAPTIVE_REPORT.out.sections)
        .mix(ch_cfdna_section)
        .mix(VARIANT_REPORT.out.sections)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Realtime mode also collects versions emitted through the "versions" topic.
    */

    if (is_realtime) {

        def topic_versions = channel.topic("versions")
            .distinct()
            .branch { entry ->
                versions_file: entry instanceof Path
                versions_tuple: true
            }

        def topic_versions_string = topic_versions.versions_tuple
            .map { process, tool, version ->
                [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
            }
            .groupTuple(by:0)
            .map { process, tool_versions ->
                tool_versions.unique().sort()
                "${process}:\n${tool_versions.join('\n')}"
            }

        ch_final_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
            .mix(topic_versions_string)

    } else {
        ch_final_versions = ch_versions
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FINAL REPORT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.realtime) {

        // Automatically detect the time stamp
        ch_title = ONTIME_TIME_RANGE.out.txt
            .map { meta, file ->
                def lines     = file.readLines()
                def start_str = lines[0].split(':\\s+')[1].trim()
                def end_str   = lines[1].split(':\\s+')[1].trim()

                def start    = parseToInstant(start_str)
                def end      = parseToInstant(end_str)
                def duration = java.time.Duration.between(start, end)
                def hours    = duration.toHours()
                def minutes  = duration.toMinutesPart()

                def time_str = "${hours}h${minutes > 0 ? String.format('%02dm', minutes) : ''}"

                tuple(meta, "${title_prefix} in realtime — ${meta.id} (atfer ${time_str} of sequencing)")
            }

    } else {

        ch_title_meta = is_cfdna
            ? ch_bed_for_calling.map { params_tuple -> params_tuple[0] }
            : MAPPING.out.bam.map { meta, _bam, _bai -> meta }

        ch_title = ch_title_meta
            .map { meta ->
                tuple(meta, "${title_prefix} — ${meta.id}")}
    }

    MIDNIGHT_REPORT(
        ch_bed_for_calling,
        ch_cfdna_full,
        ch_sections,
        ch_final_versions,
        ch_title
    )
}
