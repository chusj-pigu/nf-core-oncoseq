/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_SPLIT_BY_BED  } from '../../../modules/local/samtools/main.nf'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COVERAGE_SEPARATE {

    //TODO Add reports for coverage stats figure ?

    take:
    bam                     // channel: from mapping workflow, tuple with bam and bai
    bed                     // channel: from path read from params.bed, bed file used for adaptive sampling
    main:

    ch_versions = Channel.empty()

    ch_split_in = bam
        .combine(bed)                       // Keep padding in bed file

    SAMTOOLS_SPLIT_BY_BED(ch_split_in)

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
    versions         = ch_collated_versions              // channel: [ path(versions.yml) ]

}
