/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oncoseq: Mapping Subworkflow
    - Handles mapping of reads to reference using minimap2 and downstream processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MINIMAP2_ALIGN                         } from '../../../modules/nf-core/minimap2/align/main'         // minimap2 alignment
include { SAMTOOLS_INDEX                         } from '../../../modules/nf-core/samtools/index/main'         // Index BAM
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_CHUNK } from '../../../modules/nf-core/samtools/merge/main'         // Merge BAMs
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_FINAL } from '../../../modules/nf-core/samtools/merge/main'         // Merge BAMs
include { CRAMINO                                } from '../../../modules/nf-core/cramino/main'
include { SEQKIT_STATS                           } from '../../../modules/nf-core/seqkit/stats/main'
// functions
include { modifyMetaId           } from '../utils_nfcore_oncoseq_pipeline'
include { softwareVersionsToYAML } from '../../../subworkflows/nf-core/utils_nfcore_pipeline' // Version reporting

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN MAPPING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow MAPPING {
    // Input channels:
    //   in_ch: Channel of tuples [meta, reads] (reads can be file or directory)
    //   ref:      Channel of tuples [meta, ref, ref_fasta, ref_fai]
    take:
    in_ch       // Channel: from basecalling workflow, from --fastq if --skip_mapping is used, or from input samplesheet if skip_mapping is used
    ref       // Channel: from input samplesheet
    fai


    main:

    // Only expand in_ch if skip_basecalling is true
    if (params.skip_basecalling) {
        in_ch
            .map { meta, reads ->
                // Ensure 'reads' is a list and flatten it
                def dir_list = reads instanceof List ? reads.flatten() : [reads]
                def dir = file(dir_list[0])

                if (dir.isDirectory()) {
                    // Collect all FASTQ files with common extensions
                    def files = dir.listFiles().findAll { f ->
                        f.name ==~ /.*\.(fastq|fq)(\.gz)?$/
                    }
                    return tuple(meta, files)
                } else {
                    return tuple(meta, [dir])
                }
            }
            .set { in_ch }

            if (!params.cfdna) {
                SEQKIT_STATS(in_ch)
                ch_seqkit_out = SEQKIT_STATS.out.stats
            } else {
                ch_seqkit_out = channel.empty()
            }
        } else {
            ch_seqkit_out = channel.empty()
        }

    // Merge bams if multiple bam are provided when skip_mapping is used
    if (params.skip_mapping) {
        in_ch
            .flatMap { meta, bams -> resolveBamInput(meta, bams) }
            .branch { _meta, type, _bam, _bai ->
                single:   type == 'single'
                multi:    type == 'multi'
                to_index: type == 'to_index'
            }
            .set { bams_chunk_sep }


            SAMTOOLS_MERGE_CHUNK(
                bams_chunk_sep.multi.map { meta_chunk, _type, chunk, chunk_index -> tuple(meta_chunk, chunk, chunk_index) },
                ref.join(fai).map { meta, fasta, index -> tuple(meta, fasta, index, []) }.first(),
                'bai'
                )

            SAMTOOLS_MERGE_CHUNK.out.bam
                .join(SAMTOOLS_MERGE_CHUNK.out.index)
                .map { meta, bam, bai ->
                    def meta_restore = modifyMetaId(meta, 'remove_suffix', '', '', /_chunk\d+$/)
                    tuple(meta_restore, bam, bai)
                }
                .groupTuple()           // Group by original meta.id
                .set { final_merge_ch }

            // Count resulting bams
            ch_count = final_merge_ch
                .map { meta, bam_list, bai_list -> tuple(meta, bam_list.size(), bam_list, bai_list) }

            // Separate samples that require further merging from those that have only one resulting merged bam
            ch_count
                .branch { _meta, bam_count, _bam_list, _bai_list ->
                    single_bam: bam_count == 1
                    multi_bam: bam_count > 1
                }
                .set { ch_bam_merged }

            // Make another round of merging from those merged chunks
            SAMTOOLS_MERGE_FINAL(
                ch_bam_merged.multi_bam.map{ meta, _size, bam_list, bai_list -> tuple(meta, bam_list, bai_list) },
                ref.join(fai).map { meta, fasta, index -> tuple(meta, fasta, index, []) }.first(),
                'bai'
            )

            ch_to_index = bams_chunk_sep.to_index
                .map{ meta, _type, bam, _bai -> tuple(meta, bam) }

            ch_single_bam = bams_chunk_sep.single
                .map { meta, _type, bam, bai ->
                tuple(meta, bam, bai) }

            SAMTOOLS_INDEX(ch_to_index)

            ch_bam_indexed = ch_to_index
                .join(SAMTOOLS_INDEX.out.index)
                .mix(ch_bam_merged.single_bam.map { meta, _bam_count, bam_list, bai_list -> tuple(meta, bam_list, bai_list) })
                .mix(SAMTOOLS_MERGE_FINAL.out.bam.join(SAMTOOLS_MERGE_FINAL.out.index))

            bam_ch = ch_bam_indexed
                .mix(ch_single_bam)

            CRAMINO(bam_ch)

        } else {

            // Run minimap2 alignment
            MINIMAP2_ALIGN(
                in_ch,
                ref,
                true,
                'bai',
                false,
                false
            )

            bam_ch = MINIMAP2_ALIGN.out.bam
                .join(MINIMAP2_ALIGN.out.index)

            // Compute coverage stats
            CRAMINO(bam_ch)

        }

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
    bam      = bam_ch                               // Final sorted BAM with index
    coverage = CRAMINO.out.stats                   // Coverage stats
    seqkit   = ch_seqkit_out                       // Input fastq stats
    versions = ch_collated_versions               // All tool versions
}

// Look for sample.bam.bai or sample.bai next to the BAM
def findBamIndex(bam) {
    return [
        bam.resolveSibling("${bam.name}.bai"),
        bam.resolveSibling("${bam.baseName}.bai")
    ].find { it -> it.exists() }
}

def resolveBamInput(meta, bams) {
    // Normalise: a path, a list of paths, or a directory
    def inputs = (bams instanceof Collection ? bams.flatten() : [bams]).collect { it -> file(it) }

    def bam_files = inputs.collectMany { f ->
        f.isDirectory()
            ? f.listFiles().findAll { it -> it.name.endsWith('.bam') }.sort()
            : [f]
    }

    if (bam_files.isEmpty()) {
        throw new IllegalArgumentException("No BAM files found for sample ${meta.id}")
    }

    // [bam, bai-or-null] pairs
    def pairs = bam_files.collect { bam -> [bam, findBamIndex(bam)] }

    // Single BAM
    if (pairs.size() == 1) {
        def (bam, bai) = pairs[0]
        return [ tuple(meta, bai ? 'single' : 'to_index', bam, bai ?: []) ]
    }

    // Multiple BAMs: every one needs an index
    def missing = pairs.findAll { it -> !it[1] }.collect { it -> it[0].name }
    if (missing) {
        throw new IllegalArgumentException(
            "Missing BAM index for sample ${meta.id}: ${missing.join(', ')}")
    }

    // Chunks of 5, transposed into (bam list, bai list)
    return pairs.collate(5).withIndex().collect { chunk, i ->
        def meta_chunk = modifyMetaId(meta, 'add_suffix', '', '', "_chunk${i + 1}")
        tuple(meta_chunk, 'multi', chunk.collect { it -> it[0] }, chunk.collect { it -> it[1] })
    }
}
