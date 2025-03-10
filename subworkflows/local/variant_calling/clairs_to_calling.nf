/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLAIRS_TO_CALL         } from '../../../modules/local/clairsto/main.nf'
include  { SNPEFF_ANNOTATE        } from '../../../modules/local/snpeff/main.nf'
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

    // Choose database for annotation according to reference
    ch_databases = ref_ch.map { path ->
        if (path.contains('hg38') || path == params.genomes.GRCh38.fasta) {
            'GRCh38.p14'
        } else if (path.contains('hg19') || path == params.genomes.GRCh37.fasta) {
            'GRCh37.p13'
        } else {
            throw new IllegalArgumentException("Currently, this workflow only supports Small variant calling with reference genomes GRCh38 or GRCh37")
        }
    }

    ch_snp_annotate = CLAIRS_TO_CALL.out.vcf
        .combine(ch_databases)

    SNPEFF_ANNOTATE(ch_snp_annotate)


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
