/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SNIFFLES_CALL                     } from '../../../modules/local/sniffles/main.nf'
include { SNPEFF_ANNOTATE                   } from '../../../modules/local/snpeff/main.nf'
include { SEVERUS_TUMOR_PHASED              } from '../../../modules/local/severus/main.nf'
include { SEVERUS_TUMOR_UNPHASED            } from '../../../modules/local/severus/main.nf'
include { STELLERATOR                       } from '../../../modules/local/stellerator/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HG38 } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HG19 } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_HS1  } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_FILTERVEP              } from '../../../modules/nf-core/ensemblvep/filtervep/main.nf'
include { BCFTOOLS_SORT                     } from '../../../modules/local/bcftools/main.nf'
include { BCFTOOLS_INDEX                    } from '../../../modules/local/bcftools/main.nf'
include { BGZIP_VCF                         } from '../../../modules/local/bcftools/main.nf'
include { modifyMetaId                      } from '../utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SV_CALLING {

    take:
    bam        // channel: from phasing workflow, includes index
    ref         // reference channel with index
    snp_vcf        // channel: from phasing workflow, phased snv vcf
    vep_cache

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SNIFLLES2
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Sniffles2
    ch_in_sniffles = bam
        .join(ref)

    SNIFFLES_CALL(ch_in_sniffles)

    BCFTOOLS_SORT(SNIFFLES_CALL.out.vcf)

    ch_ref_type = ref
        .map { meta, refid, _ref_fasta, _ref_fai ->
            tuple(meta, refid) }

    // Branch ref channel to create database channel
    ch_databases = ch_ref_type.branch {
        hg38: { _meta, refid -> refid.matches('hg38|GRCh38') }
        hg19: { _meta, refid -> refid.matches('hg19|GRCh37') }
        other: true
            return 'Error'
    }

    ch_databases_hg38 = ch_databases.hg38
        .map { meta, _refid -> tuple(meta, 'GRCh38.p14') }
    ch_databases_hg19 = ch_databases.hg19
        .map { meta, _refid -> tuple(meta, 'GRCh37.p13') }

    ch_databases_ref = ch_databases_hg38
        .mix(ch_databases_hg19)

    ch_sv_annotate = BCFTOOLS_SORT.out.vcf
        .join(ch_databases_ref)
        .map { meta, output, database ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_sv_sniffles')
            tuple(new_meta, output, database) }

    SNPEFF_ANNOTATE(ch_sv_annotate)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SEVERUS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (params.realtime || params.cfdna) {

        ch_severus_in = bam
            .join(ch_ref_type)
        SEVERUS_TUMOR_UNPHASED(ch_severus_in)
        ch_vcf_severus = SEVERUS_TUMOR_UNPHASED.out.vcf
        ch_severus_versions = SEVERUS_TUMOR_UNPHASED.out.versions
    } else {

        ch_vcf_snp = snp_vcf
            .map { meta, vcf_file ->
                def new_meta = modifyMetaId(meta, 'replace', '_germline_snp_snpeff', '', '')
                tuple(new_meta, vcf_file) }

        ch_severus_in = bam
            .join(ch_vcf_snp)
            .join(ch_ref_type)

        SEVERUS_TUMOR_PHASED(ch_severus_in)
        ch_vcf_severus = SEVERUS_TUMOR_PHASED.out.vcf
        ch_severus_versions = SEVERUS_TUMOR_PHASED.out.versions
    }

    ch_vep = ch_vcf_severus
        .join(ch_ref_type)
        .branch { meta, vcf, genome ->
            hg38: genome == "hg38"
            return tuple(modifyMetaId(meta, 'add_suffix', '', '', '_sv_severus'),vcf,[])
            hg19: genome == "hg19"
            return tuple(modifyMetaId(meta, 'add_suffix', '', '', '_sv_severus'),vcf,[])
            hs1: genome == "hs1"
                return tuple(modifyMetaId(meta, 'add_suffix', '', '', '_sv_severus'),vcf,[])
            }

    ch_fasta = ch_vcf_severus
        .join(ref)
        .branch { meta, _vcf, genome, ref_fasta, _ref_index ->
            hg38: genome == "hg38"
            return tuple(meta,ref_fasta)
            hg19: genome == "hg19"
            return tuple(meta,ref_fasta)
            hs1: genome == "hs1"
                return tuple(meta,ref_fasta)
            }

    def vep_cache_resolved = vep_cache ? vep_cache : []

    ENSEMBLVEP_HG38(
        ch_vep.hg38,
        "GRCh38",
        "homo_sapiens",
        params.vep_version,
        vep_cache_resolved,
        ch_fasta.hg38,
        []
    )

    ENSEMBLVEP_HG19(
        ch_vep.hg19,
        "GRCh37",
        "homo_sapiens",
        params.vep_version,
        vep_cache_resolved,
        ch_fasta.hg19,
        []
    )

    ENSEMBLVEP_HS1(
        ch_vep.hs1,
        "CHM13",
        "homo_sapiens",
        params.vep_version,
        vep_cache_resolved,
        ch_fasta.hs1,
        []
    )

    ch_vep_to_filter = ENSEMBLVEP_HG38.out.vcf
        .mix(ENSEMBLVEP_HG19.out.vcf)
        .mix(ENSEMBLVEP_HS1.out.vcf)

    ENSEMBLVEP_FILTERVEP(
        ch_vep_to_filter,
        []
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     STELLERATOR
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_genes = channel.fromPath(params.fusion_targets, checkIfExists:true)
        .splitCsv()
        .map { row ->
            tuple(row[0], row[1])
            }

    ch_in_stellerator = bam
        .join(ch_ref_type)
        .combine(ch_genes)

    STELLERATOR(ch_in_stellerator)

    ch_out_stellerator = STELLERATOR.out.table
        .map { meta, table ->
            def table_resolved = table.readLines().size > 1 ? table : "negative"
            tuple(meta, table_resolved) }
        .unique()
        .groupTuple()
        .branch { meta, tables ->
            empty: tables.size() == 1
                return [meta, tables, "stellerator"]
            pos:   true
                return [meta, tables.findAll { it != "negative" }] }

    ch_pos_stellerator = ch_out_stellerator.pos
    ch_empty_stellerator = ch_out_stellerator.empty


    ch_to_bgzip = SNPEFF_ANNOTATE.out.vcf
        .mix(ENSEMBLVEP_FILTERVEP.out.output)

    BGZIP_VCF(ch_to_bgzip)
    BCFTOOLS_INDEX(BGZIP_VCF.out.vcf_gz)


    ch_versions = SNIFFLES_CALL.out.versions
        .mix(SNPEFF_ANNOTATE.out.versions)
        .mix(BCFTOOLS_SORT.out.versions)
        .mix(BGZIP_VCF.out.versions)
        .mix(BCFTOOLS_INDEX.out.versions)
        .mix(ch_severus_versions)

    emit:
    vcf              = BCFTOOLS_INDEX.out.vcf_tbi
    stellerator      = ch_pos_stellerator
    empty_calls      = ch_empty_stellerator
    versions         = ch_versions

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
