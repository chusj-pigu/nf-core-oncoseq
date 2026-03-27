/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLAIR3_CALL                             } from '../../../modules/local/clair3/main.nf'
include { SNPEFF_ANNOTATE                         } from '../../../modules/local/snpeff/main.nf'
include { BCFTOOLS_FILTER_REGION                  } from '../../../modules/local/bcftools/main.nf'
include { BGZIP_VCF as BGZIP_VCF_FINAL            } from '../../../modules/local/bcftools/main.nf'
include { BGZIP_VCF as BGZIP_VCF_INTER            } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FINAL  } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_RAW    } from '../../../modules/local/bcftools/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HG38       } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HG19       } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_FILTERVEP            } from '../../../modules/nf-core/ensemblvep/filtervep/main.nf'
include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { modifyMetaId                    } from '../utils_nfcore_oncoseq_pipeline' // Function to modify meta IDs

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLAIR3_CALLING {

    //TODO Add reports for coverage stats figure ?

    take:
    bam  // channel: from mapping workflow (tuple include bai)
    ref     // channel: from input samplesheet
    basecall_model       // channel: basecalling model
    bed
    vep_cache
    main:

    ch_versions = channel.empty()

    ch_ref = ref
        .map { meta, _ref, ref_fasta, ref_fai ->
            tuple(meta, ref_fasta, ref_fai) }

    ch_input_clair3 = bam
        .join(ch_ref)
        .combine(basecall_model)
        .map { meta, bamfile, bai, ref_fasta, ref_fai, model ->
            def model_str = model instanceof Path ? model.getBaseName() : model.toString()
            def model_clair3 = model_str.contains('sup')
                ? 'r1041_e82_400bps_sup_v500'
                : (model_str.contains('hac') || model_str.contains('fast'))
                    ? 'r1041_e82_400bps_hac_v500'
                    : { throw new IllegalArgumentException("Unsupported model: ${model}") }()
            tuple(meta, bamfile, bai, ref_fasta, ref_fai, model_clair3)
        }
        .join(bed)

    ch_ref_type = ref
        .map { meta, refid, _ref_fasta, _ref_fai ->
            tuple(meta, refid)
            }

    CLAIR3_CALL(ch_input_clair3)


// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     VCF ANNOTATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    // Filter for regions inside adaptive bed file

    if (params.bed != "${projectDir}/assets/NO_BED") {

        BCFTOOLS_INDEX_RAW(CLAIR3_CALL.out.vcf)
        ch_in_filter_bcftools = BCFTOOLS_INDEX_RAW.out.vcf_tbi
            .join(bed)
        BCFTOOLS_FILTER_REGION(ch_in_filter_bcftools)

        ch_clair3_out = BCFTOOLS_FILTER_REGION.out.filt_vcf
    } else {
        ch_clair3_out = CLAIR3_CALL.out.vcf
    }

    // Branch ref channel to create database channel
    ch_databases = ch_ref_type.branch {
        hg38: { meta, refid -> refid.matches('hg38|GRCh38') }
        hg19: { meta, refid -> refid.matches('hg19|GRCh37') }
        other: true
            return 'Error'
    }

    ch_vep_hg38 = ch_clair3_out
        .join(ch_databases.hg38)
        .map { meta, vcf, _refid ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_germline_snp_vep')
            tuple(new_meta,vcf, [])}

    ch_vep_hg19 = ch_clair3_out
        .join(ch_databases.hg19)
        .map { meta, vcf, _refid ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_germline_snp_vep')
            tuple(new_meta,vcf, [])}

    ch_fasta_hg38 = ref
        .filter { _meta, refid, _ref_fasta, _ref_index ->
            refid.matches('hg38|GRCh38') }
        .map { meta, _refid, ref_fasta, _ref_index ->
            tuple(meta, ref_fasta)}

    ch_fasta_hg19 = ref
        .filter { _meta, refid, _ref_fasta, _ref_index ->
            refid.matches('hg19|GRCh37') }
        .map { meta, _refid, ref_fasta, _ref_index ->
            tuple(meta, ref_fasta)}

    ENSEMBLVEP_HG38(
        ch_vep_hg38,
        "GRCh38",
        "homo_sapiens",
        params.vep_version,
        vep_cache,
        ch_fasta_hg38,
        []
    )

    ENSEMBLVEP_HG19(
        ch_vep_hg19,
        "GRCh37",
        "homo_sapiens",
        params.vep_version,
        vep_cache,
        ch_fasta_hg19,
        []
    )

    ch_vep_to_filter = ENSEMBLVEP_HG38.out.vcf
        .mix(ENSEMBLVEP_HG19.out.vcf)

    ENSEMBLVEP_FILTERVEP(
        ch_vep_to_filter,
        []
    )

    ch_databases_hg38 = ch_databases.hg38
        .map { meta, _refid -> tuple(meta, 'GRCh38.p14') }
    ch_databases_hg19 = ch_databases.hg19
        .map { meta, _refid -> tuple(meta, 'GRCh37.p13') }

    ch_databases_ref = ch_databases_hg38
        .mix(ch_databases_hg19)

    ch_snp_annotate = CLAIR3_CALL.out.vcf
        .join(ch_databases_ref)
        .map { meta, output, database ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_germline_snp_snpeff')
            tuple(new_meta, output, database)
        }

    SNPEFF_ANNOTATE(ch_snp_annotate)

    ch_vcf_final = SNPEFF_ANNOTATE.out.vcf

    // Compress and index vcf, using different names to publish only when realtime or cfdna is used (no phasing)

    if (params.realtime == null || !params.cfdna ) {
        BGZIP_VCF_INTER(ch_vcf_final)
        ch_vcf_zip = BGZIP_VCF_INTER.out.vcf_gz
        ch_versions = BGZIP_VCF_INTER.out.versions
    } else {
        BGZIP_VCF_FINAL(ch_vcf_final)
        ch_vcf_zip = BGZIP_VCF_FINAL.out.vcf_gz
        ch_versions = BGZIP_VCF_FINAL.out.versions
    }

    BCFTOOLS_INDEX_FINAL(ch_vcf_zip)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_versions = ch_versions
            .mix(CLAIR3_CALL.out.versions)
            .mix(SNPEFF_ANNOTATE.out.versions)
            .mix(BCFTOOLS_INDEX_FINAL.out.versions)
            //.mix(ENSEMBLVEP_VEP.out.versions_ensemblvep)
            //.mix(ENSEMBLVEP_FILTERVEP.out.versions_ensemblvep)

    emit:
    vcf_snpeff       = BCFTOOLS_INDEX_FINAL.out.vcf_tbi
    vcf_vep          = ENSEMBLVEP_FILTERVEP.out.output
    versions         = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
