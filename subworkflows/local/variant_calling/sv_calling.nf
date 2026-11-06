/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { WGET as WGET_TR_BED               } from '../../../modules/nf-core/wget/main'
include { WGET as WGET_PON                  } from '../../../modules/nf-core/wget/main'
include { SNIFFLES                          } from '../../../modules/nf-core/sniffles'
include { BCFTOOLS_INDEX                    } from '../../../modules/nf-core/bcftools/index'
include { SEVERUS_TUMOR_SINGLE              } from '../../../modules/local/severus/main.nf'
include { HTSLIB_BGZIPTABIX                 } from '../../../modules/nf-core/htslib/bgziptabix'
include { STELLERATOR                       } from '../../../modules/local/stellerator/main.nf'

// functions
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { modifyMetaId                 } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'
include { findCached                   } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SV_CALLING {

    take:
    bam        // channel: from phasing workflow, includes index
    ref         // reference channel with index
    genome
    snp_vcf        // channel: from phasing workflow, phased snv vcf

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SNIFLLES2
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Download the tandem repeat bed file

    ch_td_repeat = channel.of("trf.bed")
        .map { type ->
            def url = genome == "hg38" ? "https://github.com/PacificBiosciences/pbsv/raw/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed"
                    : genome == "hg19" ? "https://github.com/PacificBiosciences/pbsv/raw/master/annotations/human_hs37d5.trf.bed"
                    : "https://github.com/PacificBiosciences/pbsv/raw/master/annotations/human_chm13v2.0_maskedY_rCRS.trf.bed"

            // 4th element: local copy, or null if not in the cache
            tuple(id: url.tokenize('/')[-1], url, type, findCached(url))
        }
        .branch { _meta, _url, _type, local_file ->
            cached:   local_file != null
            download: true
        }

    WGET_TR_BED(
        ch_td_repeat.download.map { meta, url, type, _local_file -> tuple(meta, url, type) }
    )

    ch_trf_cached = ch_td_repeat.cached.map { meta, _url, _type, local_file -> tuple(meta, local_file) }

    ch_trf_bed = WGET_TR_BED.out.outfile
        .mix(ch_trf_cached)
        .first()

    SNIFFLES(
        bam.map { meta, bam_file, index ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_sv_sniffles')
            tuple(new_meta, bam_file, index)
        },
        ref,
        ch_trf_bed,
        true,
        false
    )

    ch_vcf_sniffles = SNIFFLES.out.vcf
        .join(SNIFFLES.out.tbi)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SEVERUS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Download pon files if not present in ref cache

    ch_pon_url = channel.of("tsv.gz")
        .map { type ->
            def url = genome == "hg38" ? "https://github.com/KolmogorovLab/Severus/raw/main/pon/PoN_1000G_hg38_extended.tsv.gz"
                    : genome == "hs1" ? "https://github.com/KolmogorovLab/Severus/raw/main/pon/PoN_1000G_chm13.tsv.gz" :
                    "https://github.com/KolmogorovLab/Severus/raw/main/pon/PoN_1000G_hg19.tsv.gz"

            tuple(id: url.tokenize('/')[-1], url, type, findCached(url))
        }
        .branch { _meta, _url, _type, local_file ->
            cached:   local_file != null
            download: true
        }

    ch_pon_to_download = ch_pon_url.download
        .filter { meta, _url, _type, _local_file -> meta.id != "hg19" }
        .map { meta, url, type, _local -> tuple(meta, url, type) }

    WGET_PON(
        ch_pon_to_download
    )

    ch_pon_cached = ch_pon_url.cached.map { meta, _url, _type, local_file -> tuple(meta, local_file) }

    ch_pon = WGET_PON.out.outfile
        .mix(ch_pon_cached)
        .first()

    ch_vcf_snp = snp_vcf
        .map { meta, vcf_file ->
            def new_meta = meta.id.replaceAll(/_snp_phased$/, '')
            tuple(id:new_meta, vcf_file) }
        .mix(bam.map { meta, _bamfile, _bai -> tuple(meta, []) })
        .groupTuple()
        .map { meta, vcfs -> tuple(meta, vcfs.find { it -> it } ?: []) }

    ch_severus_in = bam
        .join(ch_vcf_snp)
        .map { meta, bamfile, bai, vcf ->
            def vcf_phased = vcf ?: []
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_sv_severus')
            tuple(new_meta,bamfile,bai,vcf_phased)}

    SEVERUS_TUMOR_SINGLE(
        ch_severus_in,
        ref,
        ch_trf_bed,
        ch_pon
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     STELLERATOR
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_genes = channel.fromPath(params.fusion_targets, checkIfExists:true).first()

    ch_bam_in_stellerator = bam
        .map { meta, bam_file, bai ->
        tuple(modifyMetaId(meta, 'add_suffix', '', '', '_sv_stellerator'),bam_file,bai)
        }

    STELLERATOR(
        ch_bam_in_stellerator,
        genome,
        ch_genes
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     POST-PROCESSING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_vcf_to_idx_gz = SEVERUS_TUMOR_SINGLE.out.vcf.map { meta, vcf -> tuple(meta,vcf, [], [])}
        .mix(STELLERATOR.out.vcf.map { meta, vcf -> tuple(meta,vcf, [], [])})

    HTSLIB_BGZIPTABIX (
        ch_vcf_to_idx_gz,
        "compress",
        true,
        "vcf"
    )

    ch_vcf_final = HTSLIB_BGZIPTABIX.out.output
        .join(HTSLIB_BGZIPTABIX.out.index)
        .mix(ch_vcf_sniffles)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collect software version information from all processes for inclusion in reports
        and pipeline provenance tracking. Essential for reproducibility.
    */
    ch_versions = SEVERUS_TUMOR_SINGLE.out.versions
        .mix(STELLERATOR.out.versions)

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)

    emit:
    vcf_to_annotate  = ch_vcf_final
    versions         = ch_collated_versions

}
