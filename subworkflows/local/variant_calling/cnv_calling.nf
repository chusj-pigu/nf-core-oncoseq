/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { QDNASEQ_CALL                         } from '../../../modules/local/qdnaseq/main.nf'
include { DELLY_CNV                            } from '../../../modules/local/delly/main.nf'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_CNV } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_INDEX                       } from '../../../modules/nf-core/bcftools/index/main'
// functions
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { modifyMetaId                 } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CNV_CALLING {

    take:
    bam        // channel: from mapping workflow, includes index
    ref         // reference channel with index
    genome

    main:

    if (genome != "hs1") {
        QDNASEQ_CALL(
            bam,
            genome
        )

        ch_qdnaseq_cov_png   = QDNASEQ_CALL.out.cov_png
        ch_qdnaseq_calls_bed = QDNASEQ_CALL.out.calls_bed
        ch_qdnaseq_segs_bed  = QDNASEQ_CALL.out.segs_bed

        ch_versions = QDNASEQ_CALL.out.versions

    } else {
        ch_qdnaseq_cov_png   = channel.empty()
        ch_qdnaseq_calls_bed = channel.empty()
        ch_qdnaseq_segs_bed  = channel.empty()

        ch_versions = channel.empty()
    }

    DELLY_CNV(
        bam,
        ref
    )

    BCFTOOLS_INDEX(
        DELLY_CNV.out.bcf
    )

    ch_delly_to_query = DELLY_CNV.out.bcf
        .join(BCFTOOLS_INDEX.out.index)

    BCFTOOLS_QUERY_CNV(
        ch_delly_to_query,
        [],
        [],
        []
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collect software version information from all processes for inclusion in reports
        and pipeline provenance tracking. Essential for reproducibility.
    */
    ch_versions = ch_versions
        .mix(DELLY_CNV.out.versions)

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
    qdnaseq_plot           = ch_qdnaseq_cov_png
    qdnaseq_bed            = ch_qdnaseq_calls_bed
    qdnaseq_segs           = ch_qdnaseq_segs_bed
    delly_segs             = BCFTOOLS_QUERY_CNV.out.output
    delly_cov              = DELLY_CNV.out.cov
    versions               = ch_collated_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
