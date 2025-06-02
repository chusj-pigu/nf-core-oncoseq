/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLAIR3_CALL                  } from '../../../modules/local/clair3/main.nf'
include { SNPEFF_ANNOTATE              } from '../../../modules/local/snpeff/main.nf'
include { SNPSIFT_ANNOTATE             } from '../../../modules/local/snpeff/main.nf'
include { BGZIP_VCF                    } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX               } from '../../../modules/local/bcftools/main.nf'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

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
    clinic_database
    main:

    ch_versions = Channel.empty()

    ch_ref = ref
        .map { meta, _ref, ref_fasta, ref_fai ->
            tuple(meta, ref_fasta, ref_fai) }

    ch_input_clair3 = bam
        .join(ch_ref)
        .combine(basecall_model)
        .map { meta, bamfile, bai, ref_fasta, ref_fai, model ->
            def model_clai3 = model.contains('sup')
                ? 'r1041_e82_400bps_sup_v500'
                : (model.contains('hac') || model.contains('fast'))
                    ? 'r1041_e82_400bps_hac_v500'
                    : { throw new IllegalArgumentException("Unsupported model: ${model}") }()
            tuple(meta, bamfile, bai, ref_fasta, ref_fai, model_clai3)
        }

    CLAIR3_CALL(ch_input_clair3)

    ch_ref_type = ref
        .map { meta, refid, _ref_fasta, _ref_fai ->
            tuple(meta, refid) }

    // Branch ref channel to create database channel
    ch_databases = ch_ref_type.branch {
        hg38: { meta, refid -> refid.matches('hg38|GRCh38') }
        hg19: { meta, refid -> refid.matches('hg19|GRCh37') }
        other: true
            return 'Error'
    }

    // Generate error if reference is not hg38 nor hg19
    ch_error = ch_databases.other.map {
        throw new IllegalArgumentException("Unsupported reference genome: ${it.name}. Currently, only hg38/GRCh38 and hg19/GRCh37 are supported.")
    }

    ch_databases_hg38 = ch_databases.hg38
        .map { meta, _refid -> tuple(meta, 'GRCh38.p14') }
    ch_databases_hg19 = ch_databases.hg19
        .map { meta, _refid -> tuple(meta, 'GRCh37.p13') }

    ch_databases_ref = ch_databases_hg38
        .mix(ch_databases_hg19)

    ch_snp_annotate = CLAIR3_CALL.out.vcf
        .join(ch_databases_ref)
        .map { meta, output, database ->
            def meta_type = meta.id + '_germline'
                tuple(id:meta_type, output, database) }

    SNPEFF_ANNOTATE(ch_snp_annotate)

    ch_clin_db = clinic_database.toSortedList()

    ch_snpsift_annotate  = SNPEFF_ANNOTATE.out.vcf
        .combine(ch_clin_db)

    SNPSIFT_ANNOTATE(ch_snpsift_annotate)

    // Add clinvar into meta_id of SNPSIFT output

    ch_snipsift_out = SNPSIFT_ANNOTATE.out.vcf
        .map { meta, vcf ->
            def meta_type = meta.id + '_clinvar'
                tuple(id:meta_type, vcf) }

    ch_vcf_final = SNPEFF_ANNOTATE.out.vcf
        .mix(ch_snipsift_out)

    // Compress and index vcf :

    BGZIP_VCF(ch_vcf_final)

    BCFTOOLS_INDEX(BGZIP_VCF.out.vcf_gz)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_versions = CLAIR3_CALL.out.versions
            .mix(SNPEFF_ANNOTATE.out.versions)
            .mix(BGZIP_VCF.out.versions)
            .mix(BCFTOOLS_INDEX.out.versions)

    emit:
    vcf              = BCFTOOLS_INDEX.out.vcf_tbi
    versions         = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
