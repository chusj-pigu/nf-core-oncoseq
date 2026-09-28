/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SNIFFLES                          } from '../../../modules/nf-core/sniffles'
include { BCFTOOLS_INDEX                    } from '../../../modules/nf-core/bcftools/index'
include { SEVERUS_TUMOR                     } from '../../../modules/local/severus/main.nf'
include { HTSLIB_BGZIPTABIX as BGZIP_INTER  } from '../../../modules/nf-core/htslib/bgziptabix'

include { BCFTOOLS_VIEW                     } from '../../../modules/nf-core/bcftools/view'
include { SNPEFF_ANNOTATE                   } from '../../../modules/local/snpeff'
include { STELLERATOR                       } from '../../../modules/local/stellerator/main.nf'
include { ENSEMBLVEP_VEP                    } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_FILTERVEP              } from '../../../modules/nf-core/ensemblvep/filtervep/main.nf'
include { HTSLIB_BGZIPTABIX as BGZIP_FINAL  } from '../../../modules/nf-core/htslib/bgziptabix'
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

    def genomeMap = [
        'hg38'   : 'GRCh38',
        'GRCh38' : 'GRCh38',
        'hg19'   : 'GRCh37',
        'GRCh37' : 'GRCh37',
        'hs1'    : 'CHM13',
        'CHM13'  : 'CHM13',
        't2t'    : 'CHM13'
    ]

    ch_vcf_annotated = channel.empty()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SNIFLLES2
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_ref = ref
        .map { genome, ref_fasta, _ref_faidx ->
            tuple(genome,ref_fasta)}

    SNIFFLES(
        bam,
        ch_ref,
        [[],[]],
        true,
        false

    )

    BCFTOOLS_INDEX(SNIFFLES.out.vcf)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SEVERUS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_vcf_snp = snp_vcf
        .map { meta, vcf_file ->
            def new_meta = modifyMetaId(meta, 'replace', '_germline_snp_snpeff', '', '')
            tuple(new_meta, vcf_file) }

    ch_severus_in = bam
        .join(ch_vcf_snp)
        .combine(ch_ref)
        .map { meta, bamfile, bai, vcf, genome, _ref_fasta ->
            def vcf_phased = vcf ?: []
            tuple(meta,bamfile,bai,vcf_phased,genome.id)}

    SEVERUS_TUMOR(ch_severus_in)

    BGZIP_INTER (
        SEVERUS_TUMOR.out.vcf.map { meta, vcf -> tuple(meta,vcf, [], [])},
        "compress",
        true,
        "vcf"
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     STELLERATOR
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_genes = channel.fromPath(params.fusion_targets, checkIfExists:true)

    ch_in_stellerator = bam
        .combine(ch_ref)
        .map { meta, bamfile, bai, genome, _ref_fasta ->
            tuple(meta, bamfile, bai, genome.id )}
        .combine(ch_genes)

    STELLERATOR(ch_in_stellerator)

    ch_out_stellerator = STELLERATOR.out.vcf
        .map { meta, vcf ->
            tuple(modifyMetaId(meta, 'add_suffix', '', '', '_sv_stellerator'),vcf)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ANNOTATION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_vcf_to_annotate = SNIFFLES.out.vcf
        .join(BCFTOOLS_INDEX.out.index)
        .map { meta, vcf, tbi -> tuple(modifyMetaId(meta, 'add_suffix', '', '', '_sv_sniffles'),vcf, tbi) }
        .mix(BGZIP_INTER.out.output.join(BGZIP_INTER.out.index)
            .map { meta, vcf, tbi -> tuple(modifyMetaId(meta, 'add_suffix', '', '', '_sv_severus'),vcf, tbi) } )


    if (params.genome == "hg38" || params.genome == "hg19") {
        ch_snpeff_db = ch_ref
            .combine(ch_vcf_to_annotate)
            .map { genome, _ref_fasta, meta, vcf, _tbi ->
                def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_snpeff')
                def db = genome.id == "hg38" ? 'GRCh38.p14' : 'GRCh37.p13'
                tuple(new_meta, vcf, db) }

        // Annotate VCF
        SNPEFF_ANNOTATE(ch_snpeff_db)

        ch_vcf_annotated = SNPEFF_ANNOTATE.out.vcf
    }


    // Count variants remaining to avoid Ensembl annotation on empty vcf throwing error
    BCFTOOLS_VIEW(
        ch_vcf_to_annotate,
        [],
        [],
        []
    )

    ch_bcftools_count = BCFTOOLS_VIEW.out.vcf
        .branch { meta, vcf ->
            pos: vcf.size() > 0
                return meta
            empty: true
                return meta
            }

    ch_vcf_empty = ch_bcftools_count.empty
        .join(ch_vcf_to_annotate)

    ch_vep = ch_bcftools_count.pos
        .join(ch_vcf_to_annotate)
        .map { meta, vcf, tbi ->
           def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_vep')
           tuple(new_meta, vcf, tbi) }

    def vep_cache_resolved = vep_cache ? vep_cache : []

    ch_ref_resolved_vep = ch_vep
        .combine(ch_ref)
        .map { meta, vcf, tbi, genome, ref_fasta ->
            def genome_res = genomeMap[genome.id]
            tuple(meta, vcf, tbi, genome_res, ref_fasta) }

    ENSEMBLVEP_VEP(
        ch_vep,
        ch_ref_resolved_vep.map { _meta, _vcf, _tbi, genome, _ref_fasta -> genome },
        "homo_sapiens",
        params.vep_version,
        vep_cache_resolved,
        ch_ref_resolved_vep.map { _meta, _vcf, _tbi, genome, ref_fasta -> tuple(id:genome,ref_fasta) },
        []
    )

    ENSEMBLVEP_FILTERVEP(
        ENSEMBLVEP_VEP.out.vcf,
        []
    )

    ch_vcf_annotated_tobgzip = ch_vcf_annotated
        .mix(ENSEMBLVEP_FILTERVEP.out.output)
        .mix(ch_vcf_empty)

    BGZIP_FINAL (
        ch_vcf_annotated_tobgzip.map { meta, vcf -> tuple(meta,vcf, [], [])},
        "compress",
        true,
        "vcf"
    )

    ch_vcf_final = BGZIP_FINAL.out.output
        .join(BGZIP_FINAL.out.index)
        .view()


    ch_versions = SEVERUS_TUMOR.out.versions
        .mix(STELLERATOR.out.versions)

    emit:
    vcf              = ch_vcf_final
    stellerator      = ch_out_stellerator
    versions         = ch_versions

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
