/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oncoseq: Mapping Subworkflow
    - Handles mapping of reads to reference using minimap2 and downstream processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MARLIN_PILEUP  } from '../../../modules/local/marlin/main.nf'             // modkit mpileup
include { MARLIN_MERGE   } from '../../../modules/local/marlin/main.nf'
include { MARLIN_PREDICT } from '../../../modules/local/marlin/main.nf'
include { MARLIN_PLOT    } from '../../../modules/local/marlin/main.nf'
include { modifyMetaId   } from '../utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN MAPPING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow MARLIN {
    // Input channels:
    //   fastq_ch: Channel of tuples [meta, reads] (reads can be file or directory)
    //   ref:      Channel of tuples [meta, ref, ref_fasta, ref_fai]
    take:
    bam         // Channel: from mapping workflow
    ref       // Channel: from input samplesheet

    main:
    
    ch_ref_pileup = ref
        .map { meta, ref_id, _ref_fasta, _ref_fai ->
            tuple(meta, ref_id) }

    ch_pileup_in = bam
        .join(ch_ref_pileup)
    
    MARLIN_PILEUP(ch_pileup_in)

    ch_merge_in = MARLIN_PILEUP.out.bedmethyl
        .join(ch_ref_pileup)

    MARLIN_MERGE(ch_merge_in)

    MARLIN_PREDICT(MARLIN_MERGE.out.merged_bedmethyl)

    MARLIN_PLOT(MARLIN_PREDICT.out.pred)

    // Collect versions from all modules
    ch_versions = MARLIN_PILEUP.out.versions
        .mix(MARLIN_PREDICT.out.versions)

    emit:
    plot     = MARLIN_PLOT.out.pred_pdf                 // TODO: Quarto report
    versions = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
