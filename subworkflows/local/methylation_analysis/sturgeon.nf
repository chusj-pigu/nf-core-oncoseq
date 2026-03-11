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
    ref

    main:

    MODKIT_ADJUST(bam)

    ch_bam_adj = MODKIT_ADJUST.out.bam_adj
        .map { meta, bam_adj ->
            def meta_mod = modifyMetaId(meta, 'add_suffix', '', '', "_mod-adj")
            tuple(meta_mod, bam_adj)
            }

    SAMTOOLS_INDEX(ch_bam_adj)

    ch_in_modkit_pileup = ref 
        .map { meta, ref_type, fa, index ->
        def meta_mod = modifyMetaId(meta, 'add_suffix', '', '', "_mod-adj")
        tuple(meta_mod, fa)
        }
        .join(SAMTOOLS_INDEX.out.bamfile_index)
        .map { meta, fa, bam, bai ->
        tuple(meta, bam, bai, fa)}

    MODKIT_PILEUP(ch_in_modkit_pileup)

    STURGEON_INPUT_TOBED(MODKIT_PILEUP.out.bedmethyl)

    // Remove mod-adj from meta :

    ch_predict_in = STURGEON_INPUT_TOBED.out.dir
        .map { meta, dir ->
        def meta_restore = modifyMetaId(meta, 'remove_suffix', '', '', "_mod-adj")
        tuple(meta_restore,dir)
        }

    STURGEON_PREDICT(ch_predict_in)

    ch_type = Channel.of("Sturgeon")

    sturgeon_plot = STURGEON_PREDICT.out.pdf
        .combine(ch_type)

    sturgeon_csv = STURGEON_PREDICT.out.pred
        .combine(ch_type)

    // Collect versions from all modules
    ch_versions = MODKIT_ADJUST.out.versions
        .mix(SAMTOOLS_INDEX.out.versions)
        .mix(MODKIT_PILEUP.out.versions)
        .mix(STURGEON_INPUT_TOBED.out.versions)
        .mix(STURGEON_PREDICT.out.versions)

    emit:
    plot     = sturgeon_plot                        // TODO: Quarto report
    pred     = sturgeon_csv
    versions = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
