/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DORADO_BASECALL                           } from '../../../modules/local/dorado/main.nf'
include { DORADO_DEMULTIPLEX                        } from '../../../modules/local/dorado/main.nf'
include { SAMTOOLS_QSFILTER                         } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_PASS } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_FAIL } from '../../../modules/local/samtools/main.nf'
include { SEQKIT_STATS as SEQKIT_STATS_PASS         } from '../../../modules/local/seqkit/main.nf'
include { SEQKIT_STATS as SEQKIT_STATS_FAIL         } from '../../../modules/local/seqkit/main.nf'
include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                      } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { modifyMetaId                              } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BASECALL_MULTIPLEX {

    //TODO Add reports for read stats figure and tables

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_demux       // channel : demux samplesheet read in from --demux_samplesheet
    main:

    ch_versions = channel.empty()

    ch_demux_corrected_barcode = ch_samplesheet
        .combine(ch_demux, by:0)
        .map { meta, _pod5, _ubam, _model, _modif, sample, barcode, kit ->
            def new_barcode = kit + '_' + barcode
            tuple(meta,sample,new_barcode) }

    DORADO_BASECALL(ch_samplesheet)

    ch_to_dmux = ch_demux
        .map { meta, _sample, _barcode, kit ->
            tuple(meta,kit) }
        .join(DORADO_BASECALL.out.ubam)

    DORADO_DEMULTIPLEX(ch_to_dmux)

     split_bams_ch = DORADO_DEMULTIPLEX.out.demux_ubam
        .flatMap { meta, bams ->
            bams.collect { bam -> tuple(meta, bam) }                    // Split output of Dorado multiplex in a multiple tuples of one bam per barcode
        }
        .map { meta, bam ->
            def barcode = bam.baseName.replaceAll(/^[^_]*_/, '')
            def new_meta = [ id: "${meta.id}_${barcode}" ]
            tuple(new_meta, bam)
        }

    SAMTOOLS_QSFILTER(split_bams_ch)

     // Get the sample_ids from the demux_samplesheet that specify each barcode to each sample_id
    ch_new_sample_ids_pass = ch_demux_corrected_barcode
        .map { meta, sample_id, barcode ->
            def new_meta = meta.id + '_' + barcode
            tuple(new_meta, sample_id)    // flatten key to string
        }
        .combine(
            SAMTOOLS_QSFILTER.out.ubam_pass.map { meta, bam ->
                tuple(meta.id, bam)      // also flatten key
            },
            by: 0
        )
        .map { _meta, sampleid, ubam ->                         // Get rid of barcodes here and use real sample_id as meta
            tuple(id: sampleid, ubam) }
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')       // Add pass to meta in tuples for output naming
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_${meta_suffix}")
            tuple(new_meta, ubam)
        }

    ch_new_sample_ids_fail = ch_demux_corrected_barcode
        .map { meta, sample_id, barcode ->
            def new_meta = meta.id + '_' + barcode
            tuple(new_meta, sample_id)    // flatten key to string
        }
        .combine(
            SAMTOOLS_QSFILTER.out.ubam_fail.map { meta, bam ->
                tuple(meta.id, bam)      // also flatten key
            },
            by: 0
        )
        .map { _meta, sampleid, ubam ->                         // Get rid of barcodes here and use real sample_id as meta
            tuple(id: sampleid, ubam) }
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')       // Add fail to meta in tuples for output naming
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_${meta_suffix}")
            tuple(new_meta, ubam)
            }

    SAMTOOLS_TOFASTQ_PASS(ch_new_sample_ids_pass)
    SAMTOOLS_TOFASTQ_FAIL(ch_new_sample_ids_fail)

    SEQKIT_STATS_PASS(SAMTOOLS_TOFASTQ_PASS.out.fq)              // Read stats for passed reads
    SEQKIT_STATS_FAIL(SAMTOOLS_TOFASTQ_FAIL.out.fq)              // Reads stats for failed reads

    ch_fastq = SAMTOOLS_TOFASTQ_PASS.out.fq
        .map { meta, reads ->
            def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_pass')
            tuple(new_meta, reads)}


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_versions = DORADO_BASECALL.out.versions
        .mix(DORADO_DEMULTIPLEX.out.versions)
        .mix(SAMTOOLS_QSFILTER.out.versions)
        .mix(SAMTOOLS_TOFASTQ_PASS.out.versions)
        .mix(SAMTOOLS_TOFASTQ_FAIL.out.versions)
        .mix(SEQKIT_STATS_PASS.out.versions)
        .mix(SEQKIT_STATS_FAIL.out.versions)


    emit:
    fastq          = ch_fastq
    stats_pass     = SEQKIT_STATS_PASS.out.stats        // TODO: QUARTO REPORT
    stats_fail     = SEQKIT_STATS_FAIL.out.stats        // TODO: QUARTO REPORT
    versions       = ch_versions              // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
