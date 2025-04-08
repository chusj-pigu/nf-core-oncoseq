/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_INDEX    } from '../../../modules/local/bcftools/main.nf'
include { WHATSHAP_PHASE    } from '../../../modules/local/whatshap/main.nf'
include { WHATSHAP_HAPLOTAG } from '../../../modules/local/whatshap/main.nf'
include { WHATSHAP_STATS    } from '../../../modules/local/whatshap/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHASING_VARIANTS {

    //TODO Add reports for coverage stats figure ?

    take:
    bam  // channel: from mapping workflow (tuple include bai)
    ref_ch   // reference channel with index, output from clairs_to_calling workflow
    vcf_ch
    main:

    // Remove snv and indel from meta of vcf and separate them in different channels:
    ch_snv_vcf = vcf_ch
        .filter { meta, _vcf, _vcf_tbi -> meta.id.endsWith('_snv') }
        .map { meta, vcf, vcf_tbi ->
            def meta_restore = meta.id.replace('_snv', '')
                tuple(id:meta_restore, vcf, vcf_tbi)}

    ch_indel_vcf = vcf_ch
        .filter { meta, _vcf, _vcf_tbi -> meta.id.endsWith('_indel') }
        .map { meta, vcf, vcf_tbi ->
            def meta_restore = meta.id.replace('_indel', '')
                tuple(id:meta_restore, vcf, vcf_tbi)}

    ch_snv_clinvar_vcf = vcf_ch
        .filter { meta, _vcf, _vcf_tbi -> meta.id.contains('_snv_clinvar') }
        .map { meta, vcf, vcf_tbi ->
            def meta_restore = meta.id.replace('_snv_clinvar', '')
                tuple(id:meta_restore, vcf, vcf_tbi)}

    ch_indel_clinvar_vcf = vcf_ch
        .filter { meta, _vcf, _vcf_tbi -> meta.id.contains('_indel_clinvar') }
        .map { meta, vcf, vcf_tbi ->
            def meta_restore = meta.id.replace('_indel_clinvar', '')
                tuple(id:meta_restore, vcf, vcf_tbi)}

    // Combined bam with vcf and ref, put clinvar, indel and snv back in meta to process all tuples together:
    ch_phase_snv_in = bam
        .join(ch_snv_vcf)
        .combine(ref_ch)
        .map { meta, bamfile, bai, vcf, vcf_tbi, ref, ref_idx ->
            def meta_type = meta.id + '_snv'
                tuple(id:meta_type, bamfile, bai, vcf, vcf_tbi, ref, ref_idx) }
    ch_phase_indel_in = bam
        .join(ch_indel_vcf)
        .combine(ref_ch)
        .map { meta, bamfile, bai, vcf, vcf_tbi, ref, ref_idx ->
            def meta_type = meta.id + '_indel'
                tuple(id:meta_type, bamfile, bai, vcf, vcf_tbi, ref, ref_idx) }
    ch_phase_clin_snv_in = bam
        .join(ch_snv_clinvar_vcf)
        .combine(ref_ch)
        .map { meta, bamfile, bai, vcf, vcf_tbi, ref, ref_idx ->
            def meta_type = meta.id + '_snv_clinvar'
                tuple(id:meta_type, bamfile, bai, vcf, vcf_tbi, ref, ref_idx) }
    ch_phase_clin_indel_in = bam
        .join(ch_indel_clinvar_vcf)
        .combine(ref_ch)
        .map { meta, bamfile, bai, vcf, vcf_tbi, ref, ref_idx ->
            def meta_type = meta.id + '_indel_clinvar'
                tuple(id:meta_type, bamfile, bai, vcf, vcf_tbi, ref, ref_idx) }

    ch_phase_in = ch_phase_snv_in
        .mix(ch_phase_indel_in,
            ch_phase_clin_snv_in,
            ch_phase_clin_indel_in)

    WHATSHAP_PHASE(ch_phase_in)

    // Add phased in the meta to name correctly output

    ch_phased_vcf = WHATSHAP_PHASE.out.vcf_phased
        .map { meta, vcf ->
            def meta_phased = meta.id + '_phased'
                tuple(id:meta_phased, vcf) }

    BCFTOOLS_INDEX(ch_phased_vcf)

    // Add bam, bai, ref and red_idx from ch_pased_in by replacing vcf to prepare channel input for WHATSHAP_HAPLOTAG:

    ch_phased_indexed_vcf = ch_phase_in
        .map { meta, bamfile, bai, _vcf, _vcf_tbi, ref, ref_idx ->
            def meta_phased = meta.id + '_phased'
                tuple(id:meta_phased,bamfile, bai,ref, ref_idx)
        }
        .join(BCFTOOLS_INDEX.out.vcf_tbi)
        .filter { meta, _bamfile, _bai, _ref, _ref_idx, _vcf, _vcf_tbi -> meta.id.endsWith('_snv_phased') }            // Only has to run once per sample (SNV calls are usually enough to produce haplotagged cram and gtf)
        .map { meta, bamfile, bai, ref, ref_idx, vcf, vcf_tbi ->
            def meta_restore = meta.id.replace('_snv_phased', '')               // We want to name outputs according to sample_id only at this point
                tuple(id:meta_restore,bamfile,bai,vcf,vcf_tbi,ref,ref_idx)
        }

    WHATSHAP_HAPLOTAG(ch_phased_indexed_vcf)

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

    emit:
    phased_vcf       = WHATSHAP_PHASE.out.vcf_phased
    versions         = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
