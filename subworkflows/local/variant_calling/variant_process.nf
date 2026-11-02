/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_VIEW as BCFTOOLS_FILTER_SV  } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_SNP } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_VIEW as BCFTOOLS_RM_HEADER  } from '../../../modules/nf-core/bcftools/view/main'
include { SV_PROCESS                           } from '../../../modules/local/vcf_process/main.nf'
include { FIGENO_SV_FIGURE as FIGENO_FUSION    } from '../../../modules/local/figeno/main.nf'
include { FIGENO_SV_FIGURE as FIGENO_OTHER     } from '../../../modules/local/figeno/main.nf'
include { FIGENO_SV_FIGURE as FIGENO_TARGETS   } from '../../../modules/local/figeno/main.nf'
include { FIGENO_CIRCOS                        } from '../../../modules/local/figeno/main.nf'
include { FIGENO_PAN_CHR                       } from '../../../modules/local/figeno/main.nf'
include { QDNASEQ_PROCESS                      } from '../../../modules/local/vcf_process/main.nf'
include { ENSEMBL_VEP_TABLE                    } from '../../../modules/local/vcf_process/main.nf'
// functions
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { modifyMetaId                 } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'


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
    vep_filt_vcf              // channel containing vcf and its index
    qdnaseq_calls
    qdnaseq_segm
    sv_targets
    delly_cov
    delly_bed

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

    // only filter with >= 5 reads support for severus and sniffles

    ch_sv_to_filt = vep_filt_vcf
        .filter { meta, _vcf, _index ->
            meta.id.contains('sv') && !(meta.id =~ /stellerator|empty/)
        }
        .map { meta, vcf, index ->
            def new_meta = meta.id.replaceAll('filt', 'filt2')
            tuple(id:new_meta, vcf, index)}

    BCFTOOLS_FILTER_SV(
        ch_sv_to_filt,
        [],
        [],
        []
    )

    ch_severus_sniffles_filtered = BCFTOOLS_FILTER_SV.out.vcf
        .join(BCFTOOLS_FILTER_SV.out.index)

    ch_in_rm_header = ch_severus_sniffles_filtered
        .mix(vep_filt_vcf.filter { meta, _vcf, _index -> meta.id =~ /stellerator|empty/})

    BCFTOOLS_RM_HEADER(
        ch_in_rm_header,
        [],
        [],
        []
    )

    ch_exclude = channel.fromPath(params.sv_exclude)

    // Filter variants with HIGH and MODERATE impacts, create region file expected by figeno and tables for reports
    ch_sv_vcf = BCFTOOLS_RM_HEADER.out.vcf
        .map { meta, vcf ->
            def new_meta = meta.id.replaceFirst(/_sv_(severus|sniffles)_vep_(filt2|filt_empty)$|_sv_stellerator$/, '')
            tuple(id:new_meta, vcf)
        }
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

    ch_snp_filt = vep_filt_vcf
        .filter { meta, _vcf, _index -> meta.id.contains('snp')}

    // Table for Filtered snps :
    BCFTOOLS_QUERY_SNP(
        ch_snp_filt,
        [],
        [],
        []
    )

    ch_snp_exclude = channel.fromPath(params.snp_exclude)

    ch_snp_to_table = BCFTOOLS_QUERY_SNP.out.output
        .combine(ch_snp_exclude)

    ENSEMBL_VEP_TABLE(ch_snp_to_table)

    ch_other_sv = SV_PROCESS.out.other_tsv
        .transpose()

    ch_fusion = SV_PROCESS.out.fusion_tsv
        .transpose()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collect software version information from all processes for inclusion in reports
        and pipeline provenance tracking. Essential for reproducibility.
    */
    ch_versions = SV_PROCESS.out.versions
        .mix(FIGENO_OTHER.out.versions)

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)

    emit:
    circos_plot         = FIGENO_CIRCOS.out.figure
    sv_plot             = FIGENO_OTHER.out.figure
    fusion_plot         = FIGENO_FUSION.out.figure
    targets_plot        = FIGENO_TARGETS.out.figure
    sv_table            = ch_other_sv
    fusion_table        = ch_fusion
    panchr_plot         = FIGENO_PAN_CHR.out.figure
    snp_table           = ENSEMBL_VEP_TABLE.out.csv
    versions            = ch_collated_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
