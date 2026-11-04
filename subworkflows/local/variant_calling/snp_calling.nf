/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This section imports all necessary modules, subworkflows, and functions required
    for variant calling with ClairS-TO, annotation with SNPEff/SNPSift, and VCF
    processing with BCFtools.
*/
include { CLAIRS_TO_CALL                            } from '../../../modules/local/clairsto/main.nf'
include { CLAIR3_CALL                               } from '../../../modules/local/clair3/main.nf'
include { DEEPSOMATIC_TUMORONLY                     } from '../../../modules/local/deepsomatic/main'
include { BCFTOOLS_CONCAT                           } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_SORT                             } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_CLAIR3   } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_CLAIRSTO } from '../../../modules/nf-core/bcftools/index'

// functions

include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { modifyMetaId                 } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SNP_CALLING {

    //TODO Add reports for coverage stats figure and quality metrics visualization

    take:
    bam           // channel: from mapping workflow (tuple containing [meta, bam, bai])
    ref           // channel: reference genome information from input samplesheet (tuple containing [meta, refid, ref_fasta, ref_fai]
    ref_fai
    basecall_model      // channel: basecalling model name used for platform-specific optimizations in variant calling
    bed

    main:

    // Initialize empty channel for collecting software versions
    ch_versions = channel.empty()

    // Extract just the reference FASTA and its index from the reference channel
    ch_bam_reg= bam
        .join(bed)

    CLAIRS_TO_CALL(
        ch_bam_reg,
        ref,
        ref_fai,
        params.clairsto_model
    )

    // Add indel and snv to meta to name index correctly while preserving ts
    ch_snv_to_index = CLAIRS_TO_CALL.out.snv
        .map { meta, vcf ->
        def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_snv')
        tuple(new_meta,vcf)
        }

    ch_indel_to_index = CLAIRS_TO_CALL.out.indel
        .map { meta, vcf ->
        def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_indel')
        tuple(new_meta,vcf)
        }

    ch_vcf_to_index = ch_snv_to_index
        .mix(ch_indel_to_index)

    BCFTOOLS_INDEX_CLAIRSTO(ch_vcf_to_index)

    // Restore original meta id to join indel with snv together
    ch_snv_to_concat = ch_snv_to_index
        .join(BCFTOOLS_INDEX_CLAIRSTO.out.index)
        .map { meta, vcf, tbx ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_snv')
            tuple(new_meta,vcf,tbx)
        }

    ch_indel_to_concat = ch_indel_to_index
        .join(BCFTOOLS_INDEX_CLAIRSTO.out.index)
        .map { meta, vcf, tbx ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_indel')
            tuple(new_meta,vcf,tbx)
        }

    // Merge SNV and INDEL together
    ch_to_concat =  ch_snv_to_concat
        .mix(ch_indel_to_concat)
        .groupTuple()

    BCFTOOLS_CONCAT(ch_to_concat)

    ch_to_sort = BCFTOOLS_CONCAT.out.vcf
        .map { meta, vcf ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_snp_clairsto')
            tuple(new_meta,vcf)
        }
    BCFTOOLS_SORT(ch_to_sort)

    ch_clairsto_final = BCFTOOLS_SORT.out.vcf
        .join(BCFTOOLS_SORT.out.index)

    ch_versions = ch_versions.mix(CLAIRS_TO_CALL.out.versions)

    DEEPSOMATIC_TUMORONLY(
        ch_bam_reg,
        ref,
        ref_fai,
        [[], []]
    )

    ch_deeps_final = DEEPSOMATIC_TUMORONLY.out.vcf
        .join(DEEPSOMATIC_TUMORONLY.out.vcf_tbi)
        .map { meta, vcf, vcf_index ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_snp_deepsomatic')
            tuple(new_meta, vcf, vcf_index)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CLAIR3
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Germline variant calling
    */
    ch_model = basecall_model
        .map { model ->
            def model_str = model instanceof Path ? model.getBaseName() : model.toString()
            def model_clair3 = model_str.contains('sup')
                ? 'r1041_e82_400bps_sup_v500'
                : (model_str.contains('hac') || model_str.contains('fast'))
                    ? 'r1041_e82_400bps_hac_v500'
                    : { throw new IllegalArgumentException("Unsupported model: ${model}") }()
            model_clair3
        }
        .first()

    CLAIR3_CALL(
        ch_bam_reg,
        ref,
        ref_fai,
        ch_model
    )

    ch_to_index_clair3 = CLAIR3_CALL.out.vcf
        .map { meta, vcf ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_snp_clair3')
            tuple(new_meta,vcf)
        }

    BCFTOOLS_INDEX_CLAIR3(ch_to_index_clair3)

    ch_clair3_final = ch_to_index_clair3
        .join(BCFTOOLS_INDEX_CLAIR3.out.index)

    ch_snp_final = ch_clairsto_final
        .mix(ch_clair3_final)
        .mix(ch_deeps_final)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collect software version information from all processes for inclusion in reports
        and pipeline provenance tracking. Essential for reproducibility.
    */
    ch_versions = ch_versions
        .mix(CLAIR3_CALL.out.versions)

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
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Define the output channels that will be available to downstream processes.
        These include the annotated and indexed VCF files, version information,
        and any error channels for proper error handling.
    */

    emit:
    versions         = ch_collated_versions                       // Output collected version information for MultiQC and reports
    vcf              = ch_snp_final

}
