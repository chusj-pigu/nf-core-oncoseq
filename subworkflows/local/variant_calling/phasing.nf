/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { WHATSHAP_PHASE    } from '../../../modules/nf-core/whatshap/phase/main'
include { WHATSHAP_HAPLOTAG } from '../../../modules/nf-core/whatshap/haplotag/main'
include { WHATSHAP_STATS    } from '../../../modules/nf-core/whatshap/stats/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PHASED    } from '../../../modules/nf-core/samtools/index/main'

// functions
include { softwareVersionsToYAML       } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'
include { modifyMetaId                 } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Phase variants using WhatsHap
//
workflow PHASING_VARIANTS {


    take:
    bam         // channel: BAM files with index from mapping workflow [meta, bam, bai]
    ref         // channel: reference genome files [meta, ref_dir, ref_fasta]
    fai         // channel: reference genome index
    vcf         // channel: VCF files to be phased [meta, vcf, vcf_tbi]

    main:

    ch_vcf_to_phase = vcf
        .join(bam.map { meta, bamfile, bai ->
            def meta_joined = modifyMetaId(meta, 'add_suffix', '', '', '_snp')
            tuple(meta_joined, bamfile, bai)})
        .map { meta, vcf_file, tbi, bamfile, bai ->
            def meta_phased = modifyMetaId(meta, 'add_suffix', '', '', '_phased')
            tuple(meta_phased, vcf_file, tbi, bamfile, bai, [])}

    ch_ref = ref
        .join(fai)
        .first()

    WHATSHAP_PHASE(
        ch_vcf_to_phase,
        ch_ref
    )

    ch_phased_vcf = WHATSHAP_PHASE.out.vcf
        .join(WHATSHAP_PHASE.out.tbi)
        .join(bam.map { meta,bamfile, bai ->
            def meta_phased = modifyMetaId(meta, 'add_suffix', '', '', '_snp_phased')
            tuple(meta_phased, bamfile, bai)})

    WHATSHAP_HAPLOTAG(
        ch_phased_vcf,
        ref,
        fai,
        false
    )

    // Index the haplotype-tagged BAM files
    SAMTOOLS_INDEX_PHASED(WHATSHAP_HAPLOTAG.out.bam)

    ch_haptag_bam = WHATSHAP_HAPLOTAG.out.bam
        .join(SAMTOOLS_INDEX_PHASED.out.index)

    // Calculate phasing statistics (block lengths, switch errors, etc.)
    WHATSHAP_STATS(
        WHATSHAP_PHASE.out.vcf,
        false,
        true,
        false
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
    haptag_bam       = ch_haptag_bam                            // Haplotype-tagged BAM files with index
    phased_vcf       = WHATSHAP_PHASE.out.vcf                   // Phased VCF files
    phased_tbi       = WHATSHAP_PHASE.out.tbi
    versions         = ch_collated_versions                    // Software versions for reporting

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW SUMMARY & IDENTIFIED IMPROVEMENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WORKFLOW FLOW:
Input: BAM files, VCF files (regular + ClinVar), Reference genome
1. Separate regular variants from ClinVar variants
2. Phase variants using WhatsHap
3. Index phased VCF files
4. Create haplotype-tagged BAM files
5. Generate phasing statistics
Output: Haplotype-tagged BAMs, Phased VCFs, Statistics

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
