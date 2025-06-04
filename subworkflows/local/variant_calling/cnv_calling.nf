/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { QDNASEQ_CALL   } from '../../../modules/local/qdnaseq/main.nf'
include { BGZIP_VCF      } from '../../../modules/local/bcftools/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CNV_CALLING {

    //TODO Add reports for coverage stats figure ?

    take:
    bam        // channel: from mapping workflow, includes index
    ref         // reference channel with index
    main:

    ch_ref = ref
        .map { meta, ref_id, _ref_fasta, _ref_fai ->
            tuple(meta,ref_id) }

    ch_in_qdnaseq = bam
        .join(ch_ref)

    QDNASEQ_CALL(ch_in_qdnaseq)

    ch_versions = QDNASEQ_CALL.out.versions


    emit:
    qdnaseq_vcf      = QDNASEQ_CALL.out.call_vcf
    qdnaseq_plot     = QDNASEQ_CALL.out.cov_png             // TODO: Quarto report
    versions         = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
