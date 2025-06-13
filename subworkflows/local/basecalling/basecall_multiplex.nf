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

    ch_versions = Channel.empty()

    DORADO_BASECALL(ch_samplesheet)

    DORADO_DEMULTIPLEX(DORADO_BASECALL.out.ubam)

    demultiplex_out_ch = DORADO_DEMULTIPLEX.out.demux_ubam

    split_bams_ch = demultiplex_out_ch
        .groupTuple()
        .flatMap{ meta, ubam_files ->
            ubam_files.collect { ubam ->
            // Extract barcodes from ubam file name
                def ubam_base = ubam.baseName.replaceFirst(/^.*?(barcode\d{2}|unclassified).*$/, '$1')
            // Create a tuple with the extracted barcode and ubam file
                tuple(id:ubam_base, ubam)
                }
            }

    SAMTOOLS_QSFILTER(split_bams_ch)

    // Get the sample_ids from the demux_samplesheet that specify each barcode to each sample_id
    ch_new_sample_ids_pass = ch_demux
        .combine(SAMTOOLS_QSFILTER.out.ubam_pass, by:0)
        .map { barcode, sample, ubam ->                         // Get rid of barcodes here and use real sample_id as meta
            tuple([id: sample], ubam) }
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')       // Add pass to meta in tuples for output naming
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_${meta_suffix}")
            tuple(new_meta, ubam)
            }

    ch_new_sample_ids_fail = ch_demux
        .combine(SAMTOOLS_QSFILTER.out.ubam_fail, by:0)
        .map { barcode, sample, ubam ->                         // Get rid of barcodes here and use real sample_id as meta
            tuple([id: sample], ubam) }
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')       // Add fail to meta in tuples for output naming
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_${meta_suffix}")
            tuple(new_meta, ubam)
            }

    SAMTOOLS_TOFASTQ_PASS(ch_new_sample_ids_pass)
    SAMTOOLS_TOFASTQ_FAIL(ch_new_sample_ids_fail)

    SEQKIT_STATS_PASS(SAMTOOLS_TOFASTQ_PASS.out.fq)              // Read stats for passed reads
    SEQKIT_STATS_FAIL(SAMTOOLS_TOFASTQ_FAIL.out.fq)              // Reads stats for failed reads


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
    fastq          = SAMTOOLS_TOFASTQ_PASS.out.fq
    stats_pass     = SEQKIT_STATS_PASS.out.stats        // TODO: QUARTO REPORT
    stats_fail     = SEQKIT_STATS_FAIL.out.stats        // TODO: QUARTO REPORT
    versions       = ch_versions              // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
