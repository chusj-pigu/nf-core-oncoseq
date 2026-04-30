
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ      } from '../../../modules/local/cfdna_specific/main.nf'
include { CHOPPER_LENGTH } from '../../../modules/local/chopper/main.nf'
include { SEQKIT_STATS   } from '../../../modules/local/seqkit/main.nf'
include { modifyMetaId   } from '../utils_nfcore_oncoseq_pipeline'

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
    fastq               // Channel: output of basecalling or same as samplesheet
    max_len             // Channel : from params.max_length
    minqs               // Channel: from params.minqs

    main:

    if (params.skip_basecalling) {
        // If the input is a directory of fastq, merge them before processing them
             fastq
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

            ch_fastq = CAT_FASTQ.out.merged_fq
                .mix(ch_sep_fastq.single
                    .map{meta, files, _type -> tuple(meta, files)})

        } else {
            fastq
                .map { meta, reads ->
                    def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_pass')
                    tuple(new_meta, reads)
                }
                .set { ch_fastq }
        }

    ch_cfdna = cfdna_samplesheet
        .join(ch_fastq)
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

    CHOPPER_LENGTH(ch_reads_tofilt)

    ch_reads_tofilt_count = ch_samplesheet_to_process.to_filter
        .map { _meta, reads, _filter ->
            reads }
        .count()

    if (ch_reads_tofilt_count != 0) {
        ch_reads_filtered = CHOPPER_LENGTH.out.reads
            .map { meta, reads ->
                def meta_restore =  modifyMetaId(meta, 'remove_suffix', '', '', '_filt')
                tuple(meta_restore, reads)}

        ch_versions = CHOPPER_LENGTH.out.versions
    } else {
        ch_reads_filtered = Channel.empty()
        ch_versions = Channel.empty()
    }

    ch_reads_out_final = ch_reads_filtered
        .mix(ch_reads_nofilt)

    SEQKIT_STATS(ch_reads_out_final)

    ch_versions = ch_versions
        .mix(SEQKIT_STATS.out.versions)

    emit:
    reads     = ch_reads_out_final
    stats     = SEQKIT_STATS.out.stats
    versions  = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
