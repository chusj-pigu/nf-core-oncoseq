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

    // Before mapping, remove suffix "pass" in input from qc filtering in basecalling workflow:
    ch_mapping_in = fastq_ch
        .map { meta, reads ->
            def meta_prefix = meta.id.replace('_pass', '')
            tuple(id:meta_prefix, reads)
            }
        .combine(ref_ch)

    MINIMAP2_ALIGN(ch_mapping_in)

    SAMTOOLS_TOBAM(MINIMAP2_ALIGN.out.sam)
    SAMTOOLS_SORT_INDEX(SAMTOOLS_TOBAM.out.bamfile)

    CRAMINO_STATS(SAMTOOLS_SORT_INDEX.out.sortedbamidx)
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
    bam              = SAMTOOLS_SORT_INDEX.out.sortedbamidx
    coverage         = CRAMINO_STATS.out.stats                  // TODO: QUARTO REPORT
    versions         = ch_collated_versions              // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
