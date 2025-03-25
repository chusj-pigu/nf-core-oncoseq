/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DORADO_BASECALL                           } from '../../../modules/local/dorado/main.nf'
include { SAMTOOLS_QSFILTER                         } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_PASS } from '../../../modules/local/samtools/main.nf'
include { SAMTOOLS_TOFASTQ as SAMTOOLS_TOFASTQ_FAIL } from '../../../modules/local/samtools/main.nf'
include { SEQKIT_STATS as SEQKIT_STATS_PASS         } from '../../../modules/local/seqkit/main.nf'
include { SEQKIT_STATS as SEQKIT_STATS_FAIL         } from '../../../modules/local/seqkit/main.nf'
include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                      } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

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

    SAMTOOLS_QSFILTER(DORADO_BASECALL.out.ubam)

    // Add pass and fail to meta in tuples for output naming

    ch_ubam_pass = SAMTOOLS_QSFILTER.out.ubam_pass
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')
            def meta_full   = meta.id + '_' + meta_suffix
            tuple(id:meta_full, ubam)
            }

    ch_ubam_fail = SAMTOOLS_QSFILTER.out.ubam_fail
        .map { meta, ubam ->
            def meta_suffix = ubam.baseName.tokenize('_')[-1].replace('.bam', '')
            def meta_full   = meta.id + '_' + meta_suffix
            tuple(id:meta_full, ubam)
            }

    SAMTOOLS_TOFASTQ_PASS(ch_ubam_pass)
    SAMTOOLS_TOFASTQ_FAIL(ch_ubam_fail)

    SEQKIT_STATS_PASS(SAMTOOLS_TOFASTQ_PASS.out.fq)              // Read stats for passed reads
    SEQKIT_STATS_FAIL(SAMTOOLS_TOFASTQ_FAIL.out.fq)              // Reads stats for failed reads


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
    stats_pass     = SEQKIT_STATS_PASS.out.stats            // TODO: QUARTO REPORT
    stats_fail     = SEQKIT_STATS_FAIL.out.stats            // TODO: QUARTO REPORT
    versions       = ch_collated_versions              // channel: [ path(versions.yml) ]

}
