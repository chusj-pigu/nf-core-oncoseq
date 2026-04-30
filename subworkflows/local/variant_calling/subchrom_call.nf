/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SUBCHROM_CALL_WGS     } from '../../../modules/local/subchrom/main.nf'
include { SUBCHROM_CALL_PANEL   } from '../../../modules/local/subchrom/main.nf'
include { modifyMetaId          } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'

workflow SUBCHROM_CALL {

    //TODO Add reports for coverage stats figure ?

    take:
    bam        // channel: from mapping workflow, includes index
    ref         // reference channel with index
    vcf         // SNP vcf
    ch_panel_bin // subchrom panel bin file

    main:

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

    if (params.adaptive) {
        ch_panel_bin_processed = ch_panel_bin
    } else {
        ch_panel_bin_processed = ch_panel_bin                   // make a bed channel with a different ID so that subchrom_call_panel is skipped with joined channels
            .map { _meta,bedfile ->
                tuple(id:'skip', bedfile) }
    }

    ch_in_subchrom_panel = vcf
        .filter { meta, _vcf_file, _vcf_tbi -> meta.id.endsWith('_germline_snp_snpeff') }       // Only keep the snp file created by clair3 annotated with SnpEff
        .map { meta, vcf_file, _vcf_tbi ->
            def meta_restore = modifyMetaId(meta, 'replace', '_germline_snp_snpeff', '', '')       // Restore meta to be sample id only to join with ref and bam
                tuple(meta_restore, vcf_file)
            }
        .join(bam)
        .join(ch_ref_subchrom)
        .join(ch_panel_bin_processed)
        .map { meta, vcf_file, bam_file, bai, ref_id, ref_fasta, panel ->
            tuple(meta, bam_file, bai, vcf_file, ref_id, ref_fasta, panel) }        // Reorder the files to match input declaration

    SUBCHROM_CALL_PANEL(ch_in_subchrom_panel)

    ch_versions = SUBCHROM_CALL_WGS.out.versions
        .mix(SUBCHROM_CALL_PANEL.out.versions)


    emit:
    subchrom_plot_panel       = SUBCHROM_CALL_PANEL.out.cnv_png             // TODO: Quarto report
    subchrom_gene_plot_panel  = SUBCHROM_CALL_PANEL.out.focal_png           // TODO: Quarto report
    subchrom_plot_wgs         = SUBCHROM_CALL_WGS.out.cnv_png               // TODO: Quarto report
    subchrom_gene_plot_wgs    = SUBCHROM_CALL_WGS.out.focal_png             // TODO: Quarto report
    versions                  = ch_versions

}
