/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_INDEX    } from '../../../modules/local/bcftools/main.nf'
include { WHATSHAP_PHASE    } from '../../../modules/local/whatshap/main.nf'
include { WHATSHAP_HAPLOTAG } from '../../../modules/local/whatshap/main.nf'
include { WHATSHAP_STATS    } from '../../../modules/local/whatshap/main.nf'
include { SAMTOOLS_INDEX    } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_FAIDX } from '../../../modules/local/samtools/main.nf'
include { modifyMetaId } from '../utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHASING_VARIANTS {

    //TODO Add reports for coverage stats figure ?

    take:
    bam  // channel: from mapping workflow (tuple include bai)
    ref     // channel : from input samplesheet
    vcf_ch
    main:

    ch_ref = ref
        .map { meta, _ref, ref_fasta, ref_fai ->
            tuple(meta, ref_fasta, ref_fai) }


    // Remove snv and indel from meta of vcf and separate them in different channels:
    ch_snv_vcf = vcf_ch
        .filter { meta, _vcf, _vcf_tbi -> !meta.id.endsWith('_clinvar') }

    ch_snv_clinvar_vcf = vcf_ch
        .filter { meta, _vcf, _vcf_tbi -> meta.id.contains('_clinvar') }
        .map {meta, vcf, vcf_tbi ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_clinvar')
            tuple(new_meta, vcf, vcf_tbi) }


    ch_phase_snv_in = bam
        .join(ch_snv_vcf)
        .join(ch_ref)
        .map { _meta, bamfile, bai, meta_type, vcf, vcf_tbi, ref_fasta, ref_idx ->
                tuple(id:meta_type, bamfile, bai, vcf, vcf_tbi, ref_fasta, ref_idx) }

bam.view()

    ch_phase_clin_snv_in = bam
        .join(ch_snv_clinvar_vcf)
        .join(ch_ref)
        .map { _meta, bamfile, bai, meta_type, vcf, vcf_tbi, ref_fasta, ref_idx ->
                tuple(id:meta_type, bamfile, bai, vcf, vcf_tbi, ref_fasta, ref_idx) }

    ch_phase_in = ch_phase_snv_in
        .mix(ch_phase_clin_snv_in)

    ch_phase_in.view {
        "Phasing input channel: ${it}"
    }

    ch_phase_in = ch_phase_snv_in
        .mix(ch_phase_clin_snv_in)

    WHATSHAP_PHASE(ch_phase_in)

    // Add phased in the meta to name correctly output

    ch_phased_vcf = WHATSHAP_PHASE.out.vcf_phased
        .map { meta, vcf ->
            def meta_phased = meta.id + '_phased'
                tuple(id:meta_phased, vcf) }

    BCFTOOLS_INDEX(ch_phased_vcf)

    // Add bam, bai, ref and red_idx from ch_pased_in by replacing vcf to prepare channel input for WHATSHAP_HAPLOTAG:

    ch_phased_indexed_vcf = ch_phase_in
        .map { meta, bamfile, bai, _vcf, _vcf_tbi, ref_fasta, ref_idx ->
            def meta_phased = meta.id + '_phased'
                tuple(id:meta_phased,bamfile, bai,ref_fasta, ref_idx)
        }
        .join(BCFTOOLS_INDEX.out.vcf_tbi)
        .filter { meta, _bamfile, _bai, _ref_fasta, _ref_idx, _vcf, _vcf_tbi -> meta.id.endsWith('_snv_phased') }            // Only has to run once per sample (SNV calls are usually enough to produce haplotagged cram and gtf)
        .map { meta, bamfile, bai, ref_fasta, ref_idx, vcf, vcf_tbi ->
            def meta_restore = meta.id.replace('_snv_phased', '')               // We want to name outputs according to sample_id only at this point
                tuple(id:meta_restore,bamfile,bai,vcf,vcf_tbi,ref_fasta,ref_idx)
        }

    WHATSHAP_HAPLOTAG(ch_phased_indexed_vcf)
    SAMTOOLS_INDEX(WHATSHAP_HAPLOTAG.out.bam_hap)

    ch_whatshap_stats_in = ch_phased_indexed_vcf
        .map { meta,_bamfile,_bai,vcf,vcf_tbi,_ref,_ref_idx ->
            tuple(meta,vcf,vcf_tbi) }

    WHATSHAP_STATS(ch_whatshap_stats_in)

    //
    // Collate and save software versions
    //
    ch_versions = WHATSHAP_PHASE.out.versions
        .mix(WHATSHAP_HAPLOTAG.out.versions)
        .mix(WHATSHAP_STATS.out.versions)
        .mix(BCFTOOLS_INDEX.out.versions)

    emit:
    haptag_bam       = SAMTOOLS_INDEX.out.bamfile_index
    phased_vcf       = WHATSHAP_PHASE.out.vcf_phased
    versions         = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
