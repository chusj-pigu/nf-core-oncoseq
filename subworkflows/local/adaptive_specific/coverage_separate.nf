/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_SPLIT_BY_BED          } from '../../../modules/local/samtools/main.nf'
include { CRAMINO_STATS as CRAMINO_BG    } from '../../../modules/local/cramino/main.nf'
include { CRAMINO_STATS as CRAMINO_PANEL } from '../../../modules/local/cramino/main.nf'
include { MOSDEPTH_ADAPTIVE              } from '../../../modules/local/mosdepth/main.nf'
include { REMOVE_PADDING                 } from '../../../modules/local/adaptive_specific/main.nf'
include { paramsSummaryMap               } from 'plugin/nf-schema'
include { paramsSummaryMultiqc           } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML         } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText         } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

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

    // For now, we keep padding in bed file
    ch_split_in = params.adaptive_samplesheet == null ?
        bam.combine(bed.map{ meta,bed,padding -> bed }) :  // Combine if same bed file is used for all samples (drop proxy meta value and padding)
        bam.join(bed.map{ meta,bed,padding -> tuple(meta,bed) })       // Join if different bed files are used depending on samples (drop padding)

    SAMTOOLS_SPLIT_BY_BED(ch_split_in)

    CRAMINO_BG(SAMTOOLS_SPLIT_BY_BED.out.bg)
    CRAMINO_PANEL(SAMTOOLS_SPLIT_BY_BED.out.panel)

    // Remove padding from bed file for further coverage computations
    REMOVE_PADDING(bed)

    // Create channels for running mosdepth with different filters for each sample

    def addFilterAndCombine(meta, bam, bai, filter_type, bed_channel, flag, qual) {
        def meta_filt = meta.id + '_' + filter_type  // Add filter type to meta
        def filtered_channel = Channel.from(tuple(id: meta_filt, bam, bai))

    // Conditional logic to combine or join with bed channel
        def result_channel = params.adaptive_samplesheet == null ?
            filtered_channel.combine(bed_channel.flatten.last()) :              // Drop proxy meta value
            filtered_channel.join(bed_channel)

    // Add flag and qual to the resulting tuples
        result_channel.map { meta_filt, bam, bai, bed ->
            tuple(meta_filt, bam, bai, bed, flag, qual)
        }
    }

    // All alignments (No filters on MAPQ)
    ch_nofilt = addFilterAndCombine(SAMTOOLS_SPLIT_BY_BED.out.panel, 'nofilter', REMOVE_PADDING.out.bed, 1540, 0)
    // Primary alignments only
    ch_primary = addFilterAndCombine(SAMTOOLS_SPLIT_BY_BED.out.panel, 'primary', REMOVE_PADDING.out.bed, 1796, 0)
    // Unique alignments only (MAPQ = 60)
    ch_unique = addFilterAndCombine(SAMTOOLS_SPLIT_BY_BED.out.panel, 'unique', REMOVE_PADDING.out.bed, 1796, 60)

    mosdepth_in = ch_nofilt
        .mix(ch_primary,ch_unique)

    MOSDEPTH_ADAPTIVE(mosdepth_in)

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
