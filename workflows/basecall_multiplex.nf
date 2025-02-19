/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DORADO_BASECALL                           } from '../modules/local/dorado/main.nf'
include { DORADO_DEMULTIPLEX                        } from '../modules/local/dorado/main.nf'
include { SAMTOOLS_QSFILTER                         } from '../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_PASS } from '../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_FAIL } from '../modules/local/samtools/main.nf'
include { NANOPLOT_FASTQ                            } from '../modules/local/nanoplot/main.nf'
include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BASECALL_MULTIPLEX {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_demux       // channel : demux samplesheet read in from --demux_samplesheet
    main:

    ch_versions = Channel.empty()

    DORADO_BASECALL(ch_samplesheet)

    DORADO_DEMULTIPLEX(DORADO_BASECALL.out.ubam)

    demultiplex_out_ch = DORADO_DEMULTIPLEX.out.demux_ubam      // Capture the output channel

    split_bams_ch = demultiplex_out_ch
        .groupTuple()
        .flatMap{ meta, ubam_files ->
            ubam_files.collect { ubam ->
            // Use the baseName of the BAM file
                def ubam_base = ubam.baseName.replaceFirst(/^.*?(barcode\d{2}|unclassified).*$/, '$1')
            // Create a tuple with the extracted baseName and BAM file
                tuple(meta, ubam_base, ubam)
                }
            }

    SAMTOOLS_QSFILTER(split_bams_ch)

    // We remove the higher level sample_id used for basecalling
    ch_barcode_qsfilt_pass = SAMTOOLS_QSFILTER.out.ubam_pass
        .map { _meta, barcode, ubam ->
            tuple(barcode, ubam) }

    ch_barcode_qsfilt_fail = SAMTOOLS_QSFILTER.out.ubam_fail
        .map { _meta, barcode, ubam ->
            tuple(barcode, ubam) }

    // Get the sample_ids from the demux_samplesheet that specify each barcode to each sample_id
    ch_new_sample_ids_pass = ch_demux
        .combine(ch_barcode_qsfilt_pass, by:0)
        .map { barcode, sample, ubam ->
            tuple([id: sample], barcode, ubam) }

    ch_new_sample_ids_fail = ch_demux
        .combine(ch_barcode_qsfilt_fail, by:0)
        .map { barcode, sample, ubam ->
            tuple([id: sample], barcode, ubam) }

    SAMTOOLS_TOFASTQ_PASS(ch_new_sample_ids_pass)
    SAMTOOLS_TOFASTQ_FAIL(ch_new_sample_ids_fail)

    NANOPLOT_FASTQ(SAMTOOLS_TOFASTQ_PASS.out.fq)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'oncoseq_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }



    emit:
    fastq          = SAMTOOLS_TOFASTQ_PASS.out.fq
    nanoplot       = NANOPLOT_FASTQ.out.nanoplot
    versions       = ch_collated_versions              // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
