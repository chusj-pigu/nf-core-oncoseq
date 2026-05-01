
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ       } from '../../../modules/local/cfdna_specific/main.nf'
include { SAMTOOLS_FILTER } from '../../../modules/local/samtools/main.nf'
include { CHOPPER_LENGTH  } from '../../../modules/local/chopper/main.nf'
include { SEQKIT_STATS    } from '../../../modules/local/seqkit/main.nf'
include { modifyMetaId    } from '../utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN MAPPING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow READS_FILTER {
    // Input channels:
    //   fastq_ch: Channel of tuples [meta, reads] (reads can be file or directory)
    //   ref:      Channel of tuples [meta, ref, ref_fasta, ref_fai]
    take:
    cfdna_samplesheet    // Channel: subset of input samplesheet or demux_samplesheet
    input               // Channel: output of basecalling or same as samplesheet
    max_len             // Channel : from params.max_length
    minqs               // Channel: from params.minqs

    main:

    if (params.skip_basecalling) {
        // If the input is a directory of fastq, merge them before processing them
             input
                .map { meta, reads ->
                    // Ensure 'reads' is a list and flatten it
                    def dir_list = reads instanceof List ? reads.flatten() : [reads]
                    def dir = file(dir_list[0])

                    if (dir.isDirectory()) {
                        // Collect all FASTQ files with common extensions
                        def files = dir.listFiles().findAll { f ->
                            f.name ==~ /.*\.(fastq|fq)(\.gz)?$/
                        }
                        return tuple(meta, files, 'multi')
                    } else {
                        return tuple(meta, [dir], 'single')
                    }
                }
                .set { ch_inter_fastq }

            ch_inter_fastq
                .branch { list ->
                    single: list[2] == 'single'
                    multi: list[2] == 'multi'
                }
                .set { ch_sep_fastq }

            CAT_FASTQ(ch_sep_fastq.multi.map{meta, files, _type -> tuple(meta, files)})

            ch_input = CAT_FASTQ.out.merged_fq
                .mix(ch_sep_fastq.single
                    .map{meta, files, _type -> tuple(meta, files)})

        } else {
            input
                .map { meta, reads ->
                    def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_pass')
                    tuple(new_meta, reads)
                }
                .set { ch_input }
        }

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
        .combine(max_len)
        .combine(minqs)

    ch_reads_nofilt = ch_samplesheet_to_process.other
        .map { meta, reads, _filter ->
            tuple(meta, reads) }

    ch_reads_tofilt_count = ch_samplesheet_to_process.to_filter
        .map { _meta, reads, _filter ->
            reads }
        .count()

    if (params.skip_mapping) {

        SAMTOOLS_FILTER(ch_reads_tofilt.map{meta,reads,_len,_qs -> tuple(meta,reads)})
        ch_filt_out = SAMTOOLS_FILTER.out.bamfile
        ch_versions = SAMTOOLS_FILTER.out.versions

    } else {

        CHOPPER_LENGTH(ch_reads_tofilt)
        ch_filt_out = CHOPPER_LENGTH.out.reads
        ch_versions = CHOPPER_LENGTH.out.versions

    }


    if (ch_reads_tofilt_count != 0) {
        ch_reads_filtered = ch_filt_out
            .map { meta, reads ->
                def meta_restore =  modifyMetaId(meta, 'remove_suffix', '', '', '_filt')
                tuple(meta_restore, reads)}
    } else {
        ch_reads_filtered = channel.empty()
        ch_versions = channel.empty()
    }

    ch_reads_out_final = ch_reads_filtered
        .mix(ch_reads_nofilt)

    if (!params.skip_mapping) {

        SEQKIT_STATS(ch_reads_out_final)
        ch_versions = ch_versions
            .mix(SEQKIT_STATS.out.versions)
        ch_stats = SEQKIT_STATS.out.stats

    } else {
        ch_stats = channel.empty()
    }

    emit:
    reads     = ch_reads_out_final
    stats     = ch_stats
    versions  = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
