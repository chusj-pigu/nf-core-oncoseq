/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CLAIRS_TO_CALL               } from '../../../modules/local/clairsto/main.nf'
include { SAMTOOLS_FAIDX               } from '../../../modules/local/samtools/main.nf'
include { SNPEFF_ANNOTATE              } from '../../../modules/local/snpeff/main.nf'
include { SNPSIFT_ANNOTATE             } from '../../../modules/local/snpeff/main.nf'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

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
    clinic_database
    main:

    ch_versions = Channel.empty()

    ch_faidx_in = bam
        .combine(ref_ch)
        .first()
        .map { meta, _bamfile, _bai, ref ->
            tuple(meta,ref) }                   // We only need to run it once because we use the same reference for all samples

    if (ref_ch.first() == ref_ch.last()) {

        SAMTOOLS_FAIDX(ch_faidx_in)

        ref_idx_ch = SAMTOOLS_FAIDX.out.fasta_index
            .map { _meta, fasta_index -> fasta_index }           // Remove meta from tuple so we can join it with all samples
    } else {
        ref_idx_ch = ref_ch.last()
    }

    ch_input_clairs = bam
        .combine(ref_ch.first())
        .combine(ref_idx_ch)
        .combine(chr_list)
        .combine(model)

    CLAIRS_TO_CALL(ch_input_clairs)

    // Branch ref channel to create database channel
    ch_databases = ref_ch.first().branch {
        hg38: it.name.matches('(?i).*(hg38|GRCh38).*')
            return 'GRCh38.p14'
        hg19: it.name.matches('(?i).*(hg19|GRCh37).*')
            return 'GRCh37.p13'
        other: true
            return 'Error'
    }

    // Generate error if reference is not hg38 nor hg19
    ch_error = ch_databases.other.map {
        throw new IllegalArgumentException("Unsupported reference genome: ${it.name}. Currently, only hg38/GRCh38 and hg19/GRCh37 are supported.")
    }

    // Function to add type (indel/snv) to meta and for annotation
    def annotateMeta = { meta, output, type ->
        def meta_vcf = meta.id + "_${type}"
        tuple(id:meta_vcf, output)
    }

    if (ch_databases.hg38) {
        ch_snp_annotate = CLAIRS_TO_CALL.out.indel
            .map { meta, output -> annotateMeta(meta, output, 'indel') }
            .mix(CLAIRS_TO_CALL.out.snv
                .map { meta, output -> annotateMeta(meta, output, 'snv') })
            .combine(ch_databases.hg38)
        SNPEFF_ANNOTATE(ch_snp_annotate)
    } else if (ch_databases.hg19) {
        ch_snp_annotate = CLAIRS_TO_CALL.out.indel
            .map { meta, output -> annotateMeta(meta, output, 'indel') }
            .mix(CLAIRS_TO_CALL.out.snv
                .map { meta, output -> annotateMeta(meta, output, 'snv') })
            .combine(ch_databases.hg19)
        SNPEFF_ANNOTATE(ch_snp_annotate)
    }

    ch_database = clinic_database.toList()

    ch_snpsift_annotate  = SNPEFF_ANNOTATE.out.vcf
        .combine(ch_database)
    
    SNPSIFT_ANNOTATE(ch_snpsift_annotate)

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
