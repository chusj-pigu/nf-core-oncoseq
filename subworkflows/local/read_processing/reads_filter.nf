
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ                            } from '../../../modules/nf-core/cat/fastq/main'
include { SAMTOOLS_VIEW as SAMTOOLS_FILTER_LEN } from '../../../modules/nf-core/samtools/view/main'
include { CHOPPER                              } from '../../../modules/nf-core/chopper/main'
include { SEQKIT_STATS                         } from '../../../modules/nf-core/seqkit/stats/main'

include { modifyMetaId    } from '../utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN MAPPING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow READS_FILTER {
    // Input channels:
    take:
    cfdna_samplesheet    // Channel: subset of input samplesheet or demux_samplesheet
    input               // Channel: output of basecalling or same as samplesheet

    main:

    if (params.skip_basecalling) {

        ch_fastq = input
            .branch { meta, files ->
            multiple: files.isDirectory() && files.listFiles().findAll { f -> f.name ==~ /.*\.(fastq|fq)(\.gz)?$/ }.size() > 1
            single: true
            }

        CAT_FASTQ(ch_fastq.multiple.map{ meta, files -> tuple([id: meta.id, single_end: true], files) })

        ch_input = CAT_FASTQ.out.reads.map { meta, file -> tuple(id:meta.id, file)}
            .mix(ch_fastq.single)

    } else if (!params.skip_mapping) {
        input
            .map { meta, reads ->
                def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_pass')
                tuple(new_meta, reads)
            }
            .set { ch_input }
    }

    if (params.skip_mapping) {

        ch_bam = cfdna_samplesheet
            .join(input)
            .map { meta, _purity, filter, bam, bai ->
                tuple(meta, bam, bai, filter) }
            .branch { tup ->
                to_filter: tup[3] == "yes"                                // Filter column true
                other: true                                             // Filter column false
            }

        ch_reads_nofilt = ch_bam.other
            .map { meta, bam, bai, _filter ->
                tuple(meta,bam,bai) }

        SAMTOOLS_FILTER_LEN(
            ch_bam.to_filter.map{meta,bam,bai,_filter ->
                def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_filt")
                tuple(new_meta,bam,bai)},
            [[],[],[]],
            [[],[]],
            [[],[]],
            'bai'
        )
        ch_filt_out = SAMTOOLS_FILTER_LEN.out.bam
            .join(SAMTOOLS_FILTER_LEN.out.bai)

        ch_reads_filtered = ch_filt_out
            .map { meta, bam, bai ->
                def meta_restore =  modifyMetaId(meta, 'remove_suffix', '', '', '_filt')
                tuple(meta_restore, bam, bai)}

        ch_reads_out_final = ch_reads_filtered
            .mix(ch_reads_nofilt)

    } else {
        ch_cfdna = cfdna_samplesheet
            .join(ch_input)
            .map { meta, _purity, filter, reads ->
                tuple(meta, reads, filter) }

        ch_samplesheet_to_process = ch_cfdna.branch { tup ->
            to_filter: tup[2] == "yes"                                // Filter column true
            other: true                                             // Filter column false
        }

        ch_reads_tofilt = ch_samplesheet_to_process.to_filter
            .map { meta, reads, _filter ->
                def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_filt")
                tuple(new_meta, reads) }

        ch_reads_nofilt = ch_samplesheet_to_process.other
            .map { meta, reads, _filter ->
                tuple(meta, reads) }

        CHOPPER(
            ch_reads_tofilt,
            []
        )

        ch_filt_out = CHOPPER.out.fastq

        ch_reads_filtered = ch_filt_out
            .map { meta, reads ->
                def meta_restore =  modifyMetaId(meta, 'remove_suffix', '', '', '_filt')
                tuple(meta_restore, reads)}

        ch_reads_out_final = ch_reads_filtered
            .mix(ch_reads_nofilt)

    }

    if (!params.skip_mapping) {

        SEQKIT_STATS(ch_reads_out_final)

        ch_stats = SEQKIT_STATS.out.stats

    } else {
        ch_stats = channel.empty()
    }

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

    emit:
    reads     = ch_reads_out_final
    stats     = ch_stats
    versions  = topic_versions_string                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
