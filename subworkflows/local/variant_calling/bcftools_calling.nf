/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BCFTOOLS_MPILEUP             } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_CALL                } from '../../../modules/local/bcftools/main.nf'
include { SNPEFF_ANNOTATE              } from '../../../modules/local/snpeff/main.nf'
include { SNPSIFT_ANNOTATE             } from '../../../modules/local/snpeff/main.nf'
include { BGZIP_VCF                    } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX               } from '../../../modules/local/bcftools/main.nf'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { modifyMetaId                 } from '../utils_nfcore_oncoseq_pipeline' // Function to modify meta IDs
include { SUBCHROM_CALL_PANEL          } from '../../../modules/local/subchrom/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BCFTOOLS_CALLING {

    //TODO Add reports for coverage stats figure ?

    take:
    bam  // channel: from mapping workflow (tuple include bai)
    ref     // channel: from input samplesheet
    clinic_database
    ch_panel_bin // subchrom panel bin file
    main:

    ch_versions = Channel.empty()

    ch_ref = ref
        .map { meta, _ref, ref_fasta, ref_fai ->
            tuple(meta, ref_fasta, ref_fai) }

    ch_mpileup_in = bam
        .join(ch_ref)

    BCFTOOLS_MPILEUP(ch_mpileup_in)

    BCFTOOLS_CALL(BCFTOOLS_MPILEUP.out.bcf)

    ch_ref_type = ref
        .map { meta, refid, _ref_fasta, _ref_fai ->
            tuple(meta, refid)
            }

    // Branch ref channel to create database channel
    ch_databases = ch_ref_type.branch {
        hg38: { meta, refid -> refid.matches('hg38|GRCh38') }
        hg19: { meta, refid -> refid.matches('hg19|GRCh37') }
        other: true
            return 'Error'
    }

    ch_databases_hg38 = ch_databases.hg38
        .map { meta, _refid -> tuple(meta, 'GRCh38.p14') }
    ch_databases_hg19 = ch_databases.hg19
        .map { meta, _refid -> tuple(meta, 'GRCh37.p13') }

    ch_databases_ref = ch_databases_hg38
        .mix(ch_databases_hg19)

    ch_snp_annotate = BCFTOOLS_CALL.out.vcf
        .join(ch_databases_ref)
        .map { meta, output, database ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_baseline_snp')
            tuple(new_meta, output, database)
        }

    SNPEFF_ANNOTATE(ch_snp_annotate)

    ch_clin_db = clinic_database.toSortedList()

    ch_snpsift_annotate  = SNPEFF_ANNOTATE.out.vcf
        .combine(ch_clin_db)

    SNPSIFT_ANNOTATE(ch_snpsift_annotate)

    // Add clinvar into meta_id of SNPSIFT output

    ch_snipsift_out = SNPSIFT_ANNOTATE.out.vcf
        .map { meta, vcf ->
            def meta_clinvar = modifyMetaId(meta, 'add_suffix', '', '', '_clinvar')
            tuple(meta_clinvar, vcf)
            }

    ch_vcf_final = SNPEFF_ANNOTATE.out.vcf
        .mix(ch_snipsift_out)

    // Compress and index vcf :

    BGZIP_VCF(ch_vcf_final)

    BCFTOOLS_INDEX(BGZIP_VCF.out.vcf_gz)

    ch_subchrom_in = bam.join(BCFTOOLS_CALL.out.vcf).join(ref).join(ch_panel_bin)
        .map {meta, sc_bam, sc_bai, sc_vcf, sc_refid, sc_refpath, _reffai, sc_panelbed ->
            tuple(meta, sc_bam, sc_bai, sc_vcf, sc_refid, sc_refpath, sc_panelbed)
        }

    ch_subchrom_in.view()

    SUBCHROM_CALL_PANEL(ch_subchrom_in)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_versions = BCFTOOLS_MPILEUP.out.versions
            .mix(BCFTOOLS_CALL.out.versions)
            .mix(SNPEFF_ANNOTATE.out.versions)
            .mix(BGZIP_VCF.out.versions)
            .mix(BCFTOOLS_INDEX.out.versions)

    emit:
    vcf              = BCFTOOLS_INDEX.out.vcf_tbi
    versions         = ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
