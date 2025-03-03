/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLAIRS_TO_CALL         } from '../../../modules/local/clairsto/main.nf'
include { SAMTOOLS_FAIDX         } from '../../../modules/local/samtools/main.nf'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLAIRS_TO_CALLING {

    //TODO Add reports for coverage stats figure ?

    take:
    bam  // channel: from mapping workflow (tuple include bai)
    ref_ch   // channel: from path read from params.ref or used directly on the command line with --genome GRCh38 for example with AWS
    chr_list    // channel: list of chromosomes to include for variant calling read from params.chr_list
    model       // channel: basecalling model
    main:

    ch_versions = Channel.empty()

    ch_faidx_in = bam
        .combine(ref_ch)
        .first()
        .map { meta, _bamfile, _bai, ref ->
            tuple(meta,ref) }                   // We only need to run it once because we use the same reference for all samples

    SAMTOOLS_FAIDX(ch_faidx_in)

    ref_idx_ch = SAMTOOLS_FAIDX.out.fasta_index
        .map { _meta, fasta_index -> fasta_index }           // Remove meta from tuple so we can join it with all samples

    ch_input_clairs = bam
        .combine(ref_ch)
        .combine(ref_idx_ch)
        .combine(chr_list)
        .combine(model)

    CLAIRS_TO_CALL(ch_input_clairs)

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
