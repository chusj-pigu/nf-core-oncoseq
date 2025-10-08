/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oncoseq: Sturgeon subworkflow
    - Classifier for CNS tumor, run on methylation data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MODKIT_PILEUP        } from '../../../modules/local/modkit/main.nf'
include { SAMTOOLS_INDEX       } from '../../../modules/local/samtools/main.nf'
include { MODKIT_ADJUST        } from '../../../modules/local/modkit/main.nf'
include { STURGEON_INPUT_TOBED } from '../../../modules/local/sturgeon/main.nf'
include { STURGEON_PREDICT     } from '../../../modules/local/sturgeon/main.nf'
include { modifyMetaId         } from '../utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN MAPPING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow STURGEON {
    // Input channels:
    //   fastq_ch: Channel of tuples [meta, reads] (reads can be file or directory)
    //   ref:      Channel of tuples [meta, ref, ref_fasta, ref_fai]
    take:
    bam         // Channel: from mapping workflow, bam + bai

    main:

    MODKIT_ADJUST(bam)

    ch_bam_adj = MODKIT_ADJUST.out.bam_adj
        .map { meta, bam_adj ->
            def meta_mod = modifyMetaId(meta, 'add_suffix', '', '', "_mod-adj")
            tuple(meta_mod, bam_adj)
            }

    SAMTOOLS_INDEX(ch_bam_adj)

    MODKIT_PILEUP(SAMTOOLS_INDEX.out.bamfile_index)

    STURGEON_INPUT_TOBED(MODKIT_PILEUP.out.bedmethyl)

    STURGEON_PREDICT(STURGEON_INPUT_TOBED.out.dir)

    // Collect versions from all modules
    ch_versions = MODKIT_ADJUST.out.versions
        .mix(SAMTOOLS_INDEX.out.versions)
        .mix(MODKIT_PILEUP.out.versions)
        .mix(STURGEON_INPUT_TOBED.out.versions)
        .mix(STURGEON_PREDICT.out.versions)

    emit:
    plot     = STURGEON_PREDICT.out.dir                // TODO: Quarto report
    versions = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
