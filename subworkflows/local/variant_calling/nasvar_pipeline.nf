/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { NASVAR_PIPELINE   } from '../../../modules/local/nasvar/main.nf'
include { HTSLIB_BGZIPTABIX } from '../../../modules/nf-core/htslib/bgziptabix'
include { SAMTOOLS_FAIDX    } from '../../../modules/nf-core/samtools/faidx/main.nf'
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

    // nasvar only accepts unzipped reference fasta:
    ch_ref = ref
        .branch { _meta, fasta ->
            bgzipped: fasta.name.endsWith(".gz")
            unzipped: true
        }

    // Only keep the index for references that are already unzipped
    ch_fai_unzipped = ch_ref.unzipped
        .join(ref_index)
        .map { meta, _fasta, fai -> tuple(meta, fai) }

    ch_ref_to_uncompress = ch_ref.bgzipped
        .join(ref_index)
        .map { meta, fasta, fai -> tuple(meta, fasta, fai, []) }

    HTSLIB_BGZIPTABIX(
        ch_ref_to_uncompress,
        "decompress",
        false,
        "fa"
    )

    SAMTOOLS_FAIDX(
        HTSLIB_BGZIPTABIX.out.output.map { meta, fasta -> tuple(meta, fasta, []) },
        false
    )

    ch_ref_final = HTSLIB_BGZIPTABIX.out.output
        .mix(ch_ref.unzipped)

    ch_fai_final = SAMTOOLS_FAIDX.out.fai
        .mix(ch_fai_unzipped)

    NASVAR_PIPELINE(
        ch_in,
        ch_ref_final,
        ch_fai_final
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
