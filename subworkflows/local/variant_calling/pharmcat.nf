/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { WGET                                   } from '../../../modules/nf-core/wget/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_SNP } from '../../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_VIEW as BCFTOOLS_COUNT        } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_INDEX                         } from '../../../modules/nf-core/bcftools/index/main'
include { PHARMCAT_VCFPREPROCESSOR               } from '../../../modules/nf-core/pharmcat/vcfpreprocessor/main'
include { PHARMCAT_MATCHER                       } from '../../../modules/nf-core/pharmcat/matcher/main'
include { PHARMCAT_PHENOTYPER                    } from '../../../modules/nf-core/pharmcat/phenotyper/main'
include { PHARMCAT_REPORTER                      } from '../../../modules/nf-core/pharmcat/reporter/main'

// Functions
include { modifyMetaId             } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline/main'
include { findCached               } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline/main'
include { softwareVersionsToYAML   } from '../../../subworkflows/nf-core/utils_nfcore_pipeline/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW: ONLY COMPATIBLE WITH HG38
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHARMCAT {

    take:
    vcf        // Annotated and merged snp vcf file
    ref         // reference channel with index
    ref_index

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DOWNLOAD REQUIRED FILES IF NOT PRESENT IN REF_CACHE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_pgx_vcf = channel.of("vcf.bgz")
        .map { type ->
            def url = "https://github.com/PharmGKB/PharmCAT/raw/development/pharmcat_positions.vcf.bgz"
            tuple(id: "pharmcat_positions", url, type, findCached(url))
        }
        .branch { _meta, _url, _type, local_file ->
            cached:   local_file != null
            download: true
        }

    ch_pgx_vcf_index = channel.of("vcf.bgz.csi")
        .map { type ->
            def url = "https://github.com/PharmGKB/PharmCAT/raw/development/pharmcat_positions.vcf.bgz.csi"
            tuple(id: "pharmcat_positions", url, type, findCached(url))
        }
        .branch { _meta, _url, _type, local_file ->
            cached:   local_file != null
            download: true
        }

    ch_to_download = ch_pgx_vcf.download
        .mix(ch_pgx_vcf_index.download)
        .map { meta, url, type, _local_file -> tuple(meta, url, type) }

    WGET(ch_to_download)

    // Cached files, reshaped to (meta, file) like the WGET output
    ch_pgx_cached = ch_pgx_vcf.cached
        .mix(ch_pgx_vcf_index.cached)
        .map { meta, _url, _type, local_file -> tuple(meta, local_file) }

    ch_pgx_joined = WGET.out.outfile
        .mix(ch_pgx_cached)
        .groupTuple(size: 2)
        .map { meta, files ->
            def csi_pgx = files.find { it -> it.name.endsWith('.csi') }
            def vcf_pgx = files.find { it -> !it.name.endsWith('.csi') }
            tuple(meta, vcf_pgx, csi_pgx)
        }
        .first()

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PHARMCAT PIPELINE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Filter out low quality genotypes first
    ch_vcf_tofilt = vcf
        .map { meta, vcf_file, index ->
        def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_pharmcat')
        tuple(new_meta, vcf_file, index)}

    BCFTOOLS_FILTER_SNP(
        ch_vcf_tofilt
    )

    ch_vcf_to_preprocess = BCFTOOLS_FILTER_SNP.out.vcf
        .join(BCFTOOLS_FILTER_SNP.out.index)

    PHARMCAT_VCFPREPROCESSOR(
        ch_vcf_to_preprocess,
        ref,
        ref_index,
        ch_pgx_joined,
        [[], [], []]
    )

    BCFTOOLS_INDEX(
        PHARMCAT_VCFPREPROCESSOR.out.preprocessed_vcf
    )

    ch_preprocessed_vcf = PHARMCAT_VCFPREPROCESSOR.out.preprocessed_vcf
        .join(BCFTOOLS_INDEX.out.index)

    BCFTOOLS_COUNT(
        ch_preprocessed_vcf,
        [],
        [],
        []
    )

    ch_bcftools_count = BCFTOOLS_COUNT.out.vcf
        .branch { meta, vcf_file ->
            pos: vcf_file.size() > 0
                return meta
            empty: true
                return meta
            }
    ch_vcf_to_preprocess = ch_bcftools_count.pos
        .join(ch_preprocessed_vcf)

    PHARMCAT_MATCHER(
        ch_vcf_to_preprocess,
        []
    )

    PHARMCAT_PHENOTYPER(
        PHARMCAT_MATCHER.out.matcher_json.map { meta, json -> tuple(meta, json, [])}
    )

    PHARMCAT_REPORTER(
        PHARMCAT_PHENOTYPER.out.phenotyper_json
    )

     /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Collect software version information from all processes for inclusion in reports
        and pipeline provenance tracking. Essential for reproducibility.
    */

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

    ch_collated_versions = softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)

    emit:
    versions               = ch_collated_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
