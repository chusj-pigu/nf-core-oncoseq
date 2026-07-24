/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_FILTER                      } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_SNP } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_VIEW as BCFTOOLS_RM_HEADER  } from '../../../modules/local/bcftools/main.nf'
include { SV_PROCESS                           } from '../../../modules/local/vcf_process/main.nf'
include { STELLERATOR_PROCESS                  } from '../../../modules/local/vcf_process/main.nf'
include { FIGENO_SV_FIGURE as FIGENO_FUSION    } from '../../../modules/local/figeno/main.nf'
include { FIGENO_SV_FIGURE as FIGENO_OTHER     } from '../../../modules/local/figeno/main.nf'
include { FIGENO_SV_FIGURE as FIGENO_TARGETS   } from '../../../modules/local/figeno/main.nf'
include { FIGENO_CIRCOS                        } from '../../../modules/local/figeno/main.nf'
include { FIGENO_PAN_CHR                       } from '../../../modules/local/figeno/main.nf'
include { QDNASEQ_PROCESS                      } from '../../../modules/local/vcf_process/main.nf'
include { ENSEMBL_VEP_TABLE                    } from '../../../modules/local/vcf_process/main.nf'
include { modifyMetaId                         } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT_PROCESS {

    //TODO Add reports for coverage stats figure ?

    take:
    bam        // channel: from mapping workflow, includes index
    bed
    sv_vcf
    stellerator
    qdnaseq_calls
    qdnaseq_segm
    sv_targets
    delly_cov
    delly_bed
    snp_filt
    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     CNV CIRCOS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_cnv_process = qdnaseq_calls
        .join(qdnaseq_segm)

    // Process qdnaseq outputs to match expected ControlFREEC input of figeno:
    QDNASEQ_PROCESS(ch_cnv_process)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SV FILTERING AND PROCESSING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    BCFTOOLS_FILTER(sv_vcf)

    ch_severus_sniffles_filtered = BCFTOOLS_FILTER.out.filt_vcf

    ch_in_rm_header = ch_severus_sniffles_filtered
        .mix(stellerator)

    BCFTOOLS_RM_HEADER(ch_in_rm_header)

    ch_exclude = channel.fromPath(params.sv_exclude)

    // Filter variants with HIGH and MODERATE impacts, create region file expected by figeno and tables for reports
    ch_sv_vcf = BCFTOOLS_RM_HEADER.out.vcf
        .map { meta, vcf ->
            def new_meta = modifyMetaId(meta, 'replace', '_sv_severus', '', '')
            def meta_replace = modifyMetaId(new_meta, 'replace', '_sv_sniffles', '', '')
            def meta_final = modifyMetaId(meta_replace, 'replace', '_sv_stellerator', '', '')
            tuple(meta_final, vcf)}
        .groupTuple()
        .map { meta, list ->
            tuple(meta, list[0], list[1], list[2])}


    ch_to_process = ch_sv_vcf
        .join(bed)
        .combine(sv_targets)
        .combine(ch_exclude)

    SV_PROCESS(ch_to_process)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     VIZUALISATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_to_circos = SV_PROCESS.out.figeno_table
        .join(QDNASEQ_PROCESS.out.cnv_file)
        .join(QDNASEQ_PROCESS.out.ratio_file)

    FIGENO_CIRCOS(ch_to_circos)

    // Figure for structural variants

    ch_figeno_fusion = bam
        .join(SV_PROCESS.out.figeno_table)
        .join(SV_PROCESS.out.fusion_txt)

    ch_figeno_sv = bam
        .join(SV_PROCESS.out.figeno_table)
        .join(SV_PROCESS.out.other_txt)

    ch_figeno_targets = bam
        .join(SV_PROCESS.out.figeno_table)
        .join(SV_PROCESS.out.targets)

    FIGENO_OTHER(ch_figeno_sv)

    FIGENO_FUSION(ch_figeno_fusion)

    FIGENO_TARGETS(ch_figeno_targets)

    ch_delly_figeno = delly_cov
        .join(delly_bed)

    FIGENO_PAN_CHR(ch_delly_figeno)

    // Table for Filtered snps :
    BCFTOOLS_QUERY_SNP(snp_filt)

    ch_snp_exclude = channel.fromPath(params.snp_exclude)

    ch_snp_to_table = BCFTOOLS_QUERY_SNP.out.bed
        .combine(ch_snp_exclude)

    ENSEMBL_VEP_TABLE(ch_snp_to_table)

    ch_other_sv = SV_PROCESS.out.other_tsv
        .transpose()

    ch_fusion = SV_PROCESS.out.fusion_tsv
        .transpose()

    ch_versions = BCFTOOLS_FILTER.out.versions
        .mix(SV_PROCESS.out.versions)
        .mix(FIGENO_OTHER.out.versions)

    emit:
    circos_plot         = FIGENO_CIRCOS.out.figure
    sv_plot             = FIGENO_OTHER.out.figure
    fusion_plot         = FIGENO_FUSION.out.figure
    targets_plot        = FIGENO_TARGETS.out.figure
    sv_table            = ch_other_sv
    fusion_table        = ch_fusion
    panchr_plot         = FIGENO_PAN_CHR.out.figure
    snp_table           = ENSEMBL_VEP_TABLE.out.csv
    versions            = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
