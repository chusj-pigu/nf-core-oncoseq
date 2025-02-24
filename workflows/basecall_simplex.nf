/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DORADO_BASECALL                           } from '../modules/local/dorado/main.nf'
include { SAMTOOLS_QSFILTER                         } from '../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_PASS } from '../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_FAIL } from '../modules/local/samtools/main.nf'
include { NANOPLOT_UBAM                             } from '../modules/local/nanoplot/main.nf'
include { QUARTO_FIGURE                             } from '../modules/local/quarto/main.nf'
include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../subworkflows/local/utils_nfcore_oncoseq_pipeline'

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

    ch_versions = Channel.empty()

    DORADO_BASECALL(ch_samplesheet)

    ch_basecall_out = DORADO_BASECALL.out.ubam
        .map { meta, ubam ->
            tuple(meta, meta.id, ubam)         // Make a mock barcode variable as sample_id (meta) to regularize with demultiplex workflow
            }

    NANOPLOT_UBAM(DORADO_BASECALL.out.ubam)    // Run Nanoplot on the output of Dorado directly using Min QS of 10

    SAMTOOLS_QSFILTER(ch_basecall_out)


    SAMTOOLS_TOFASTQ_PASS(SAMTOOLS_QSFILTER.out.ubam_pass)
    SAMTOOLS_TOFASTQ_FAIL(SAMTOOLS_QSFILTER.out.ubam_fail)

    // Gather the nanoplot tables to input to R script
    ch_read_stats = NANOPLOT_UBAM.out.txt
        .flatMap { meta, tables ->
            tables.collect { table ->
                tuple(meta, table) // Emit a tuple for each file path
            }
        }
        .map { meta, table ->
                def lines = table.splitText()
                def num_reads = lines[5].split(/\s+/)[3].replaceAll(',', '').toDouble()
                def tot_bases = lines[8].split(/\s+/)[2].replaceAll(',', '').toDouble()
                tuple(meta, num_reads, tot_bases)
            }
        .collectFile { item ->
            def sample_id = item[0].id // Extract 'sample1' from [id: 'sample1']
            [ "read_stats.txt", sample_id + '\t' + item[1] + '\t' + item[2] +  '\n']
        }

   //QUARTO_FIGURE

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
    nanoplot       = NANOPLOT_UBAM.out.figure
    versions       = ch_collated_versions              // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
