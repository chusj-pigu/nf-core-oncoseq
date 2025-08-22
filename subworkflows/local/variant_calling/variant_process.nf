/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_INDEX                    } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_FILTER_SUPPORT           } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_FILTER_ID                } from '../../../modules/local/bcftools/main.nf'
include { SV_PROCESS                        } from '../../../modules/local/vcf_process/main.nf'
include { FIGENO_SV_FIGURE as FIGENO_FUSION } from '../../../modules/local/figeno/main.nf'
include { FIGENO_SV_FIGURE as FIGENO_OTHER  } from '../../../modules/local/figeno/main.nf'
include { FIGENO_CIRCOS                     } from '../../../modules/local/figeno/main.nf'
include { BGZIP_VCF                         } from '../../../modules/local/bcftools/main.nf'
include { QDNASEQ_PROCESS                   } from '../../../modules/local/vcf_process/main.nf'
include { modifyMetaId                      } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT_PROCESS {

    //TODO Add reports for coverage stats figure ?

    take:
    bam        // channel: from mapping workflow, includes index
    sv_vcf
    qdnaseq_calls
    qdnaseq_segm
    main:

    // Circos with cnv

    ch_cnv_process = qdnaseq_calls
        .join(qdnaseq_segm)

    // Process qdnaseq outputs to match expected ControlFREEC input of figeno:
    QDNASEQ_PROCESS(ch_cnv_process)

    // Filter all variants with support >4
    BCFTOOLS_FILTER_SUPPORT(sv_vcf)

    // Filter variants with HIGH and MODERATE impacts, create region file expected by figeno and list of selected IDS
    SV_PROCESS(BCFTOOLS_FILTER_SUPPORT.out.filt_vcf)

    // Create a filtered vcf file with only HIGH and MODERATE effects variants with support > 4
    ch_to_filter_id = sv_vcf
        .join(SV_PROCESS.out.filt_ids)

    BCFTOOLS_FILTER_ID(ch_to_filter_id)

    BGZIP_VCF(BCFTOOLS_FILTER_ID.out.filt_vcf)

    BCFTOOLS_INDEX(BGZIP_VCF.out.vcf_gz)

    ch_to_circos = BCFTOOLS_INDEX.out.vcf_tbi
        .map { meta, vcf, index ->
            def meta_restore = modifyMetaId(meta, 'replace', '_sv', '', '')
            tuple(meta_restore, vcf, index) }
        .join(QDNASEQ_PROCESS.out.cnv_file)
        .join(QDNASEQ_PROCESS.out.ratio_file)

    FIGENO_CIRCOS(ch_to_circos)

    // Figure for structural variants

    ch_bam = bam
        .map {meta, bamfile, bai ->
           def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_sv')
           tuple(new_meta, bamfile, bai)}             // Match meta_id with the vcf and region files

    ch_figeno_fusion = ch_bam
        .join(sv_vcf)
        .join(SV_PROCESS.out.fusion_txt)

    ch_figeno_sv = ch_bam
        .join(sv_vcf)
        .join(SV_PROCESS.out.indel_txt)

    FIGENO_OTHER(ch_figeno_sv)

    FIGENO_FUSION(ch_figeno_fusion)

    ch_versions = BCFTOOLS_FILTER_SUPPORT.out.versions
        .mix(BCFTOOLS_FILTER_ID.out.versions)
        .mix(SV_PROCESS.out.versions)
        .mix(BCFTOOLS_INDEX.out.versions)
        .mix(BGZIP_VCF.out.versions)
        .mix(FIGENO_OTHER.out.versions)

    emit:
    circos_plot         = FIGENO_CIRCOS.out.figure              // TODO: Quarto report
    sv_plot             = FIGENO_OTHER.out.figure               // TODO: Quarto report
    fusion_plot         = FIGENO_FUSION.out.figure              // TODO: Quarto report
    versions            = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
