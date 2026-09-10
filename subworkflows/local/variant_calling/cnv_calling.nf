/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { QDNASEQ_CALL                         } from '../../../modules/local/qdnaseq/main.nf'
include { DELLY_CNV                            } from '../../../modules/local/delly/main.nf'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_CNV } from '../../../modules/local/bcftools/main.nf'
include { SUBCHROM_CALL_WGS     } from '../../../modules/local/subchrom/main.nf'
include { modifyMetaId          } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CNV_CALLING {

    take:
    bam        // channel: from mapping workflow, includes index
    vcf         // snp vcf from Clair3
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

    if (!params.realtime) {
        ch_ref_subchrom = ref
            .map { meta, ref_id, ref_fasta, _ref_fai ->
                tuple(meta, ref_id, ref_fasta) }

        ch_in_subchrom_wgs = vcf
            .filter { meta, _vcf_file, _vcf_tbi -> meta.id.endsWith('_germline_snp_snpeff') }       // Only keep the snp file created by clair3 annotated with SnpEff
            .map { meta, vcf_file, _vcf_tbi ->
                def meta_restore = modifyMetaId(meta, 'replace', '_germline_snp_snpeff', '', '')       // Restore meta to be sample id only to join with ref
                    tuple(meta_restore, vcf_file)
                    }
            .join(ch_ref_subchrom)

        SUBCHROM_CALL_WGS(ch_in_subchrom_wgs)

        ch_subchrom_png = SUBCHROM_CALL_WGS.out.cnv_png
        ch_subchrom_focal = SUBCHROM_CALL_WGS.out.focal_png
        ch_subchrom_versions = SUBCHROM_CALL_WGS.out.versions

    } else {
        ch_subchrom_png = channel.empty()
        ch_subchrom_focal = channel.empty()
        ch_subchrom_versions = channel.empty()
    }

    ch_versions = QDNASEQ_CALL.out.versions
        .mix(DELLY_CNV.out.versions)
        .mix(BCFTOOLS_QUERY_CNV.out.versions)
        .mix(ch_subchrom_versions)

    emit:
    qdnaseq_plot           = QDNASEQ_CALL.out.cov_png
    qdnaseq_bed            = QDNASEQ_CALL.out.calls_bed
    qdnaseq_segs           = QDNASEQ_CALL.out.segs_bed
    delly_segs             = BCFTOOLS_QUERY_CNV.out.bed
    delly_cov              = DELLY_CNV.out.cov
    subchrom_plot_wgs      = ch_subchrom_png
    subchrom_gene_plot_wgs = ch_subchrom_focal
    versions               = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
