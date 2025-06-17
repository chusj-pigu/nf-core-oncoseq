/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { QDNASEQ_CALL       } from '../../../modules/local/qdnaseq/main.nf'
include { SUBCHROM_CALL_WGS  } from '../../../modules/local/subchrom/main.nf'

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
    vcf         // SNP vcf
    main:

    ch_ref_qdnaseq = ref
        .map { meta, ref_id, _ref_fasta, _ref_fai ->
            tuple(meta,ref_id) }

    ch_in_qdnaseq = bam
        .join(ch_ref_qdnaseq)

    QDNASEQ_CALL(ch_in_qdnaseq)

    ch_ref_subchrom = ref
        .map { meta, ref_id, ref_fasta, _ref_fai ->
            tuple(meta,ref_id, ref_fasta) }

    ch_in_subchrom = vcf
        .filter { meta, _vcf_file, _vcf_tbi -> meta.id.endsWith('germline_snp') }       // Only keep the snp file created by clair3 annotated with SnpEff
        .map { meta, vcf_file, _vcf_tbi ->
            def meta_restore = meta.id.replaceAll('_germline_snp', '')       // Restore meta to be sample id only to join with ref
                tuple(id:meta_restore, vcf_file) }
        .join(ch_ref_subchrom)

    SUBCHROM_CALL_WGS(ch_in_subchrom)

    ch_versions = QDNASEQ_CALL.out.versions
        .mix(ch_in_subchrom)


    emit:
    qdnaseq_vcf         = QDNASEQ_CALL.out.call_vcf
    qdnaseq_plot        = QDNASEQ_CALL.out.cov_png             // TODO: Quarto report
    subchrom_plot       = SUBCHROM_CALL_WGS.out.cnv_png
    subchrom_gene_plot  = SUBCHROM_CALL_WGS.out.focal_png
    versions            = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
