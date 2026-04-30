/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { TIDEHUNTER_FASTA } from '../../../modules/local/tidehunter/main.nf'
include { TIDEHUNTER_TAB   } from '../../../modules/local/tidehunter/main.nf'
include { SEQKIT_REPLACE   } from '../../../modules/local/seqkit/main.nf'
include { PIGZ_FA          } from '../../../modules/local/cfdna_specific/main.nf'
include { modifyMetaId     } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TIDEHUNTER_CONCENSUS {

    //TODO Add reports for coverage stats figure ?

    take:
    reads

    main:

    ch_reads = reads
        .map { meta, reads ->
        def new_meta = modifyMetaId(meta, 'add_suffix', '', '', "_cons")
        tuple(new_meta, reads)}

    TIDEHUNTER_FASTA(reads)

    SEQKIT_REPLACE(TIDEHUNTER_FASTA.out.fa)

    PIGZ_FA(SEQKIT_REPLACE.out.fasta)

    ch_fasta_out = PIGZ_FA.out.fa
        .map { meta, reads ->
        def meta_restore = modifyMetaId(meta, 'remove_suffix', '', '', "_cons") 
        tuple(meta_restore, reads)}

    TIDEHUNTER_TAB(reads)

    ch_versions = TIDEHUNTER_FASTA.out.versions
        .mix(TIDEHUNTER_TAB.out.versions)
        .mix(SEQKIT_REPLACE.out.versions)

    emit:
    fastq               = ch_fasta_out
    summary             = TIDEHUNTER_TAB.out.tsv     // TODO: Quarto report
    versions            = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
