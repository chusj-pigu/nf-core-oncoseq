
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
    samplesheet         // Channel: input samplesheet
    fastq               // Channel: output of basecalling or same as samplesheet
    max_len             // Channel : from params.max_length
    minqs               // Channel: from params.minqs

    main:

    if (params.skip_basecalling) {
        // If the input is a directory of fastq, merge them before processing them
             samplesheet
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
                .mix(ch_inter_fastq.single.map{meta, files, _type -> tuple(meta, files)})

        } else {
            fastq
                .map { meta, reads ->
                    def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_pass')
                    tuple(new_meta, reads)
                }
                .set { ch_fastq }
        }

    if (params.demux) {
        ch_samplesheet_cfdna = samplesheet
            .join(ch_fastq)
            .map { meta, _input, _kit, _purity, filter, reads ->
                tuple(meta, reads, filter) }
    } else {
        ch_samplesheet_cfdna = samplesheet
            .join(ch_fastq)
            .map { meta, _input, _purity, filter, reads ->
                tuple(meta, reads, filter) }
    }

    ch_samplesheet_to_process = ch_samplesheet_cfdna.branch {
        to_filter: it[2]                                // Filter column true
        no_filter: !it[2]                               // Filter column false
    }

    ch_reads_tofilt = ch_samplesheet_to_process.to_filter
        .map { meta, reads, _filter ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_filt")
            tuple(new_meta, reads) }
        .combine(max_len)
        .combine(minqs)

    ch_reads_nofilt = ch_samplesheet_to_process.no_filter
        .map { meta, reads, _filter ->
            tuple(meta, reads) }

    CHOPPER_LENGTH(ch_reads_tofilt)

    ch_chopper_result = CHOPPER_LENGTH.out.reads.ifEmpty(null)

    if (ch_chopper_result != null) {
        SEQKIT_STATS(ch_chopper_result)
        ch_reads_filtered = ch_chopper_result
            .map { meta, reads ->
                def meta_restore =  modifyMetaId(meta, 'remove_suffix', '', '', '_filt')
                tuple(meta_restore, reads)}

        ch_versions = CHOPPER_LENGTH.out.version
            .mix(SEQKIT_STATS.out.version)
    } else {
        ch_reads_filtered = ch_chopper_result
        ch_versions = Channel.empty()
    }

    ch_reads_out_final = ch_reads_filtered
        .mix(ch_reads_nofilt)

    emit:
    reads     = ch_reads_out_final
    versions  = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
