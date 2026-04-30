/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { QDNASEQ_CALL                         } from '../../../modules/local/qdnaseq/main.nf'
include { DELLY_CNV                            } from '../../../modules/local/delly/main.nf'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_CNV } from '../../../modules/local/bcftools/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CNV_CALLING {

    take:
    bam        // channel: from mapping workflow, includes index
    ref         // reference channel with index
    main:

    ch_ref_qdnaseq = ref
        .map { meta, ref_id, _ref_fasta, _ref_fai ->
            tuple(meta,ref_id) }

    ch_in_qdnaseq = bam
        .join(ch_ref_qdnaseq)

    QDNASEQ_CALL(ch_in_qdnaseq)

    ch_ref_delly = ref
        .map { meta, ref_id, ref_fasta, _ref_fai ->
            tuple(meta,ref_id,ref_fasta) }

    ch_in_delly = bam
        .join(ch_ref_delly)

    DELLY_CNV(ch_in_delly)
    BCFTOOLS_QUERY_CNV(DELLY_CNV.out.bcf)

    ch_versions = QDNASEQ_CALL.out.versions
        .mix(DELLY_CNV.out.versions)
        .mix(BCFTOOLS_QUERY_CNV.out.versions)

    emit:
    qdnaseq_plot        = QDNASEQ_CALL.out.cov_png             // TODO: Quarto report
    qdnaseq_bed         = QDNASEQ_CALL.out.calls_bed
    qdnaseq_segs        = QDNASEQ_CALL.out.segs_bed
    delly_segs          = BCFTOOLS_QUERY_CNV.out.bed
    delly_cov           = DELLY_CNV.out.cov
    versions            = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
