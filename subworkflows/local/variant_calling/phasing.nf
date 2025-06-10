/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// BCFtools for VCF indexing
include { BCFTOOLS_INDEX    } from '../../../modules/local/bcftools/main.nf'

// WhatsHap suite for phasing operations
include { WHATSHAP_PHASE    } from '../../../modules/local/whatshap/main.nf'
include { WHATSHAP_HAPLOTAG } from '../../../modules/local/whatshap/main.nf'
include { WHATSHAP_STATS    } from '../../../modules/local/whatshap/main.nf'

// SAMtools for BAM operations
include { SAMTOOLS_INDEX    } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_FAIDX } from '../../../modules/local/samtools/main.nf'

// Utility function for metadata manipulation
include { modifyMetaId } from '../utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Phase variants using WhatsHap
//
// This workflow performs variant phasing using WhatsHap to determine which variants
// are on the same chromosome copy (haplotype). It processes both regular variants
// and ClinVar-annotated variants, producing phased VCFs and haplotype-tagged BAMs.
//
workflow PHASING_VARIANTS {


    take:
    bam         // channel: BAM files with index from mapping workflow [meta, bam, bai]
    ref         // channel: reference genome files [meta, ref_dir, ref_fasta, ref_fai]
    vcf_ch      // channel: VCF files to be phased [meta, vcf, vcf_tbi]

    main:

    // Extract reference FASTA and index from reference channel
    // TODO: Validate that reference index exists and is current
    ch_ref = ref
        .map { meta, _ref, ref_fasta, ref_fai ->
            tuple(meta, ref_fasta, ref_fai) }

    // Separate regular variants from ClinVar-annotated variants
    // Regular variants (somatic_snp, germline_snp) for phasing
    ch_snv_vcf = vcf_ch
        .filter { meta, _vcf, _vcf_tbi -> !meta.id.endsWith('_clinvar') }
        .map { meta, vcf, vcf_tbi ->
            // Remove variant type suffixes to get base sample ID, but keep original for tracking
            def base_meta = modifyMetaId(meta, 'replace', '_somatic_snp', '', '')
            base_meta = modifyMetaId(base_meta, 'replace', '_germline_snp', '', '')
            tuple(base_meta, meta, vcf, vcf_tbi)
        }

    // ClinVar variants for separate processing
    ch_snv_clinvar_vcf = vcf_ch
        .filter { meta, _vcf, _vcf_tbi -> meta.id.contains('_clinvar') }
        .map {meta, vcf, vcf_tbi ->
            def base_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_clinvar')
            tuple(base_meta, meta, vcf, vcf_tbi)
        }
        .map { base_meta, meta, vcf, vcf_tbi ->
            // Remove variant type suffixes to get base sample ID, but keep original for tracking
            base_meta = modifyMetaId(base_meta, 'replace', '_somatic_snp', '', '')
            base_meta = modifyMetaId(base_meta, 'replace', '_germline_snp', '', '')
            tuple(base_meta, meta, vcf, vcf_tbi)
        }

    // Prepare input channels for phasing by joining BAM, VCF, and reference files
    // Regular variants channel
    ch_phase_snv_in = bam
        .join(ch_snv_vcf)
        .join(ch_ref)
        .map { _meta, bamfile, bai, base_meta, vcf, vcf_tbi, ref_fasta, ref_idx ->
                tuple(base_meta, bamfile, bai, vcf, vcf_tbi, ref_fasta, ref_idx)
            }

    // ClinVar variants channel
    ch_phase_clin_snv_in = bam
        .join(ch_snv_clinvar_vcf)
        .join(ch_ref)
        .map { _meta, bamfile, bai, base_meta, vcf, vcf_tbi, ref_fasta, ref_idx ->
                tuple(base_meta, bamfile, bai, vcf, vcf_tbi, ref_fasta, ref_idx) }

    // Combine both channels for unified phasing
    ch_phase_in = ch_phase_snv_in
        .mix(ch_phase_clin_snv_in)

    // STEP 1: Phase variants using WhatsHap
    // This determines which variants are on the same haplotype
    WHATSHAP_PHASE(ch_phase_in)

    // STEP 2: Index the phased VCF files
    // Add '_phased' suffix to metadata for proper file naming
    ch_phased_vcf = WHATSHAP_PHASE.out.vcf_phased
        .map { meta, vcf ->
            def meta_phased = modifyMetaId(meta, 'add_suffix', '', '', '_phased')
            tuple(meta_phased, vcf) }

    // Index phased VCF files for downstream analysis
    BCFTOOLS_INDEX(ch_phased_vcf)

    // STEP 3: Prepare input for haplotype tagging
    // Reconstruct channel with BAM, reference, and indexed phased VCF
    // TODO: Simplify this complex channel transformation
    ch_phased_indexed_vcf = ch_phase_in
        .map { meta, bamfile, bai, _vcf, _vcf_tbi, ref_fasta, ref_idx ->
            def meta_phased = modifyMetaId(meta, 'add_suffix', '', '', '_phased')
            tuple(meta_phased, bamfile, bai, ref_fasta, ref_idx)
        }
        .join(BCFTOOLS_INDEX.out.vcf_tbi)
        // Only process non-ClinVar variants for haplotagging (avoids duplicate processing)
        .filter { meta, _bamfile, _bai, _ref_fasta, _ref_idx, _vcf, _vcf_tbi -> !meta.id.endsWith('_clinvar_phased') }
        .map { meta, bamfile, bai, ref_fasta, ref_idx, vcf, vcf_tbi ->
            // Restore original sample ID for output naming
            // def meta_restore = modifyMetaId(meta, 'replace', '_somatic_snp_phased', '', '')
            // meta_restore = modifyMetaId(meta_restore, 'replace', '_germline_snp_phased', '', '')
            tuple(meta, bamfile, bai, vcf, vcf_tbi, ref_fasta, ref_idx)
        }

    // STEP 4: Create haplotype-tagged BAM files
    // Tag reads in BAM with haplotype information based on phased variants
    WHATSHAP_HAPLOTAG(ch_phased_indexed_vcf)
    
    // Index the haplotype-tagged BAM files
    SAMTOOLS_INDEX(WHATSHAP_HAPLOTAG.out.bam_hap)

    // STEP 5: Generate phasing statistics
    // Extract VCF files for statistics calculation
    // TODO: Add more comprehensive phasing quality metrics
    ch_whatshap_stats_in = ch_phased_indexed_vcf
        .map { meta,_bamfile,_bai,vcf,vcf_tbi,_ref,_ref_idx ->
            tuple(meta,vcf,vcf_tbi) }

    // Calculate phasing statistics (block lengths, switch errors, etc.)
    WHATSHAP_STATS(ch_whatshap_stats_in)

    
    ch_versions = WHATSHAP_PHASE.out.versions
        .mix(WHATSHAP_HAPLOTAG.out.versions)
        .mix(WHATSHAP_STATS.out.versions)
        .mix(BCFTOOLS_INDEX.out.versions)

    emit:
    haptag_bam       = SAMTOOLS_INDEX.out.bamfile_index    // Haplotype-tagged BAM files with index
    phased_vcf       = WHATSHAP_PHASE.out.vcf_phased       // Phased VCF files
    versions         = ch_versions                          // Software versions for reporting

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW SUMMARY & IDENTIFIED IMPROVEMENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WORKFLOW FLOW:
Input: BAM files, VCF files (regular + ClinVar), Reference genome
1. Separate regular variants from ClinVar variants
2. Phase variants using WhatsHap
3. Index phased VCF files
4. Create haplotype-tagged BAM files
5. Generate phasing statistics
Output: Haplotype-tagged BAMs, Phased VCFs, Statistics

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
