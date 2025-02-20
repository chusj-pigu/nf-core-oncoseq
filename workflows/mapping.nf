/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MINIMAP2_ALIGN         } from '../modules/local/minimap2/main.nf'
include { SAMTOOLS_TOBAM         } from '../modules/local/samtools/main.nf'
include { SAMTOOLS_SORT          } from '../modules/local/samtools/main.nf'
include { SAMTOOLS_INDEX         } from '../modules/local/samtools/main.nf'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAPPING {

    take:
    fastq_ch  // channel: from basecalling workflow or from --fastq if --skip_mapping is used
    ref_ch   // channel: from path read from params.ref or used directly on the command line with --genome GRCh38 for example with AWS
    main:

    ch_versions = Channel.empty()

    ch_mapping_in = fastq_ch
        .combine(ref_ch)

    MINIMAP2_ALIGN(ch_mapping_in)

    SAMTOOLS_TOBAM(MINIMAP2_ALIGN.out.sam)
    SAMTOOLS_SORT(SAMTOOLS_TOBAM.out.bamfile)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sortedbam)

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
    versions       = ch_collated_versions              // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
