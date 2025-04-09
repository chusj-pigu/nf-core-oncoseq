/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MINIMAP2_ALIGN         } from '../../../modules/local/minimap2/main.nf'
include { SAMTOOLS_TOBAM         } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_SORT_INDEX    } from '../../../modules/local/samtools/main.nf'
include { CRAMINO_STATS          } from '../../../modules/local/cramino/main.nf'
include { QUARTO_TABLE           } from '../../../modules/local/quarto/main.nf'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAPPING {

    //TODO Add reports for coverage stats figure ?

    take:
    fastq_ch  // channel: from basecalling workflow or from --fastq if --skip_mapping is used
    ref_ch   // channel: from path read from params.ref or used directly on the command line with --genome GRCh38 for example with AWS
    main:

    ch_versions = Channel.empty()

    ch_ref = ref_ch
        .map { ref, _ref_idx -> ref }
    // Before mapping, remove suffix "pass" in input from qc filtering in basecalling workflow:
    ch_mapping_in = fastq_ch
        .map { meta, reads ->
            def meta_prefix = meta.id.replace('_pass', '')
            tuple(id:meta_prefix, reads)
            }
        .combine(ch_ref)

    MINIMAP2_ALIGN(ch_mapping_in)

    SAMTOOLS_TOBAM(MINIMAP2_ALIGN.out.sam)
    SAMTOOLS_SORT_INDEX(SAMTOOLS_TOBAM.out.bamfile)

    CRAMINO_STATS(SAMTOOLS_SORT_INDEX.out.sortedbamidx)

    // // Take mean coverage only from summary file of mosdepth to reduce file size loaded into R:
    // ch_mosdepth_coverage = MOSDEPTH_GENERAL.out.summary
    //     .map { meta, table ->
    //     // Read the file content as a list of lines
    //         def lines = table.readLines()
    //         def coverage = lines[-1].tokenize('\t')[3].toDouble()    // Last line and only take mean coverage column (4th)
    //         tuple(meta, coverage)
    //     }
    //     .collectFile { item ->
    //         def sample_id = item[0].id // Extract 'sample1' from [id: 'sample1']
    //         [ "coverage.txt", sample_id + '\t' + item[1] + '\n']
    //     }

    // QUARTO_TABLE( ch_mosdepth_coverage_table,
    //     "Mean coverage",
    //     "T",
    //     "Mosdepth_coverage",
    //     "mosdepth-general")

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_versions = MINIMAP2_ALIGN.out.versions
        .mix(SAMTOOLS_TOBAM.out.versions)
        .mix(SAMTOOLS_SORT_INDEX.out.versions)
        .mix(CRAMINO_STATS.out.versions)



    emit:
    bam              = SAMTOOLS_SORT_INDEX.out.sortedbamidx
    coverage         = CRAMINO_STATS.out.stats              // TODO: QUARTO REPORT
    versions         = ch_versions              // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
