/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DORADO_BASECALL                           } from '../../../modules/local/dorado/main.nf'
include { SAMTOOLS_QSFILTER                         } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ                          } from '../../../modules/local/samtools/main.nf'
include { HTSLIB_BGZIPTABIX                         } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { SEQKIT_STATS                              } from '../../../modules/nf-core/seqkit/stats/main'

include { softwareVersionsToYAML                    } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { modifyMetaId                              } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BASECALL_SIMPLEX {

    //TODO Add reports for read stats figure and tables

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    DORADO_BASECALL(ch_samplesheet)

    SAMTOOLS_QSFILTER(DORADO_BASECALL.out.ubam)

    // Add pass and fail to meta in tuples for output naming

    ch_ubam_pass = SAMTOOLS_QSFILTER.out.ubam_pass
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_${meta_suffix}")
            tuple(new_meta, ubam)
            }

    ch_ubam_fail = SAMTOOLS_QSFILTER.out.ubam_fail
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_${meta_suffix}")
            tuple(new_meta, ubam)
            }

    ch_to_fastq = ch_ubam_pass
        .mix(ch_ubam_fail)

    SAMTOOLS_TOFASTQ(ch_to_fastq)

    SEQKIT_STATS(SAMTOOLS_TOFASTQ.out.fq)

    ch_stats_pass = SEQKIT_STATS.out.stats
        .filter { meta, _tbl -> meta.id.contains('pass')}
        .map { meta, table ->
        def meta_restore = modifyMetaId(meta, 'remove_suffix', '', '', '_pass')
        tuple(meta_restore, table)}

    ch_fastq_pass = SAMTOOLS_TOFASTQ.out.fq
        .filter { meta, _reads -> meta.id.contains('pass')}
        .map { meta, reads ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_pass')
            tuple(new_meta, reads)}


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_versions = DORADO_BASECALL.out.versions
        .mix(SAMTOOLS_QSFILTER.out.versions)
        .mix(SAMTOOLS_TOFASTQ.out.versions)

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
    fastq          = ch_fastq_pass
    stats_pass     = ch_stats_pass
    versions       = ch_collated_versions

}
