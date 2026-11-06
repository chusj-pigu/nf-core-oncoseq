/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NASVAR_PIPELINE   } from '../../../modules/local/nasvar/main.nf'
include { modifyMetaId      } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW: ONLY COMPATIBLE WITH HG38 AND HS1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NASVAR {

    take:
    bam        // channel: from mapping workflow, includes index
    bed        // snp vcf from Clair3
    bed_nopadding
    ref         // reference channel with index
    ref_index

    main:
    ch_in = bam
        .join(bed)
        .join(bed_nopadding)
        .map { meta, bamfile, bai, bedfile, _padding, _lowfid, bednopad ->
            tuple(meta, bamfile, bai, bedfile, bednopad)}

    NASVAR_PIPELINE(
        ch_in,
        ref,
        ref_index
    )

    ch_versions = NASVAR_PIPELINE.out.versions

    emit:
    json                   = NASVAR_PIPELINE.out.json
    baf_plot               = NASVAR_PIPELINE.out.karyotype_baf_svg
    versions               = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
