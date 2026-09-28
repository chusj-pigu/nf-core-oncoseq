/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NASVAR_COVERAGE        } from '../../../modules/local/nasvar/main.nf'
include { NASVAR_KARYOTYPE       } from '../../../modules/local/nasvar/main.nf'
include { NASVAR_MAF             } from '../../../modules/local/nasvar/main.nf'
include { BEDTOOLS_INTERSECT_MAF } from '../../../modules/local/nasvar/main.nf'
include { modifyMetaId           } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow KARYOTYPE {

    take:
    bam        // channel: from mapping workflow, includes index
    bed        // snp vcf from Clair3
    ref         // reference channel with index
    main:

    ch_ref_id = ref
        .filter { _meta, ref_id, _ref_fasta, _ref_fai -> ref_id == 'hs1' || ref_id == 'hg38' }
        .map { meta, ref_id, _ref_fasta, _ref_fai ->
            tuple(meta,ref_id) }

    ch_in_coverage = bam
        .join(ch_ref_id)

    NASVAR_COVERAGE(ch_in_coverage)

    ch_chr_names = bed
        .map { meta, bedfile, _padding, _low_fidelity ->
            def chr= bedfile.readLines().first().split('\t')[0].contains('chr') ? true : false
            tuple(meta, chr)
        }
        .view()

    ch_bed = bed
        .map { meta, bedfile, _padding, _low_fidelity ->
            tuple(meta, bedfile)
        }

    // Intersect the MAF sites with the target bed file
    ch_in_maf_prepare = ch_bed
        .join(ch_chr_names)
        .join(ch_ref_id)

    BEDTOOLS_INTERSECT_MAF(ch_in_maf_prepare)

    ch_maf_in = bam
        .join(ch_bed)
        .join(BEDTOOLS_INTERSECT_MAF.out.maf_sites)
        .join(ch_ref_id)

    NASVAR_MAF(ch_maf_in)

    ch_in_karyotype = NASVAR_COVERAGE.out.cov
        .join(ch_ref_id)
        .join(NASVAR_MAF.out.maf)
        .join(BEDTOOLS_INTERSECT_MAF.out.maf_sites)

    NASVAR_KARYOTYPE(ch_in_karyotype)

    ch_versions = NASVAR_COVERAGE.out.versions
        .mix(BEDTOOLS_INTERSECT_MAF.out.versions)
        .mix(NASVAR_MAF.out.versions)
        .mix(NASVAR_KARYOTYPE.out.versions)


    emit:
    karyotype              = NASVAR_KARYOTYPE.out.karyo_json
    baf_plot               = NASVAR_KARYOTYPE.out.karyo_baf
    versions               = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
