/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This section imports all necessary modules, subworkflows, and functions required
    for variant calling with ClairS-TO, annotation with SNPEff/SNPSift, and VCF
    processing with BCFtools.
*/
include { BCFTOOLS_VIEW as BCFTOOLS_COUNT } from '../../../modules/nf-core/bcftools/view'
include { HTSLIB_BGZIPTABIX               } from '../../../modules/nf-core/htslib/bgziptabix'
include { ENSEMBLVEP_VEP                  } from '../../../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_FILTERVEP            } from '../../../modules/nf-core/ensemblvep/filtervep/main.nf'

// functions
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { modifyMetaId                 } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCF_ANNOTATION {

    //TODO Add reports for coverage stats figure and quality metrics visualization

    take:
    vcf          // channel: from variant calling subworkflow containing meta, vcf and its index (excluding stellerator vcf)
    ref           // channel: reference genome information from input samplesheet (tuple containing [refid, fasta])
    vep_cache
    genome
    main:

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     VCF ANNOTATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    // Count variants remaining to avoid Ensembl annotation if vcf empty
    BCFTOOLS_COUNT(
        vcf,
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

    ch_vep = ch_bcftools_count.pos
        .join(vcf)
        .map { meta, vcf_file, _vcf_index ->
            def meta_annotate = modifyMetaId(meta, 'add_suffix', '', '', '_vep')
            tuple(meta_annotate, vcf_file, [])
        }


    def vep_cache_resolved = vep_cache ? vep_cache : []

    // Convert genome id to ensembl genome ID
    def genomeMap = [
        'hg38'   : 'GRCh38',
        'hg19'   : 'GRCh37',
        'hs1'    : 'T2T-CHM13v2.0'
    ]

    def genome_resolved = genomeMap[genome]

    // Resolve species nomenclature according to vep cache

    def speciesMap = [
        'GRCh38'        : 'homo_sapiens',
        'GRCh37'        : 'homo_sapiens',
        'T2T-CHM13v2.0' : 'homo_sapiens_gca009914755v4'
    ]

    def species_resolved = speciesMap[genome_resolved]

    ENSEMBLVEP_VEP(
        ch_vep,
        genome_resolved,
        species_resolved,
        params.vep_version,
        vep_cache_resolved,
        ref,
        []
    )

    ENSEMBLVEP_FILTERVEP(
        ENSEMBLVEP_VEP.out.vcf,
        []
    )

    ch_vep_final = ENSEMBLVEP_VEP.out.vcf
        .join(ENSEMBLVEP_VEP.out.tbi)

    ch_vepfilt_to_bgzip = ENSEMBLVEP_FILTERVEP.out.output
        .map { meta, vcf_file ->
            def new_meta = modifyMetaId(meta, 'add_suffix', '', '', '_filt')
            tuple(new_meta, vcf_file, [], [])}

    HTSLIB_BGZIPTABIX(
        ch_vepfilt_to_bgzip,
        "compress",
        true,
        "vcf"
    )

    // Make a channel for empty vcf

    ch_sv_empty = ch_bcftools_count.empty
        .join(vcf)
        .map { meta, vcf_file, _vcf_index ->
            def meta_annotate = modifyMetaId(meta, 'add_suffix', '', '', '_vep_filt_empty')
            tuple(meta_annotate, vcf_file, [])
        }

    ch_vcf_filt_final = HTSLIB_BGZIPTABIX.out.output
        .join(HTSLIB_BGZIPTABIX.out.index)
        .mix(ch_sv_empty)

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
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Define the output channels that will be available to downstream processes.
        These include the annotated and indexed VCF files, version information,
        and any error channels for proper error handling.
    */

    emit:
    vcf_vep          = ch_vep_final
    vcf_filt         = ch_vcf_filt_final
    versions         = ch_collated_versions                       // Output collected version information for MultiQC and reports

}
