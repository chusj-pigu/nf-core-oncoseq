/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SNIFFLES_CALL   } from '../../../modules/local/sniffles/main.nf'
include { SNPEFF_ANNOTATE } from '../../../modules/local/snpeff/main.nf'
include { BCFTOOLS_SORT   } from '../../../modules/local/bcftools/main.nf'
include { BGZIP_VCF       } from '../../../modules/local/bcftools/main.nf'
include { modifyMetaId    } from '../utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SV_CALLING {

    //TODO Add reports for coverage stats figure ?

    take:
    cram        // channel: from phasing workflow, includes index
    ref         // reference channel with index
    main:


    ch_in_sniffles = cram
        .join(ref)

    SNIFFLES_CALL(ch_in_sniffles)

    ch_ref_type = ref
        .map { meta, refid, _ref_fasta, _ref_fai ->
            tuple(meta, refid) }

    // Branch ref channel to create database channel
    ch_databases = ch_ref_type.branch {
        hg38: { _meta, refid -> refid.matches('hg38|GRCh38') }
        hg19: { _meta, refid -> refid.matches('hg19|GRCh37') }
        other: true
            return 'Error'
    }

    ch_databases_hg38 = ch_databases.hg38
        .map { meta, _refid -> tuple(meta, 'GRCh38.p14') }
    ch_databases_hg19 = ch_databases.hg19
        .map { meta, _refid -> tuple(meta, 'GRCh37.p13') }

    ch_databases_ref = ch_databases_hg38
        .mix(ch_databases_hg19)

    ch_sv_annotate = SNIFFLES_CALL.out.vcf
        .join(ch_databases_ref)
        .map { meta, output, database ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_sv')
            tuple(new_meta, output, database) }

    SNPEFF_ANNOTATE(ch_sv_annotate)
    BCFTOOLS_SORT(SNPEFF_ANNOTATE.out.vcf)
    BGZIP_VCF(BCFTOOLS_SORT.out.vcf)

    ch_versions = SNIFFLES_CALL.out.versions


    emit:
    sv_vcf           = BGZIP_VCF.out.vcf_gz
    versions         = ch_versions

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
