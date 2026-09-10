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
include { MINIMAP2_ALIGN                         } from '../../../modules/local/minimap2/main.nf'         // minimap2 alignment
include { SAMTOOLS_TOBAM                         } from '../../../modules/local/samtools/main.nf'         // Convert SAM to BAM
include { SAMTOOLS_SORT as SAMTOOLS_SORT    } from '../../../modules/local/samtools/main.nf'         // Sort BAM
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SUB  } from '../../../modules/local/samtools/main.nf'         // Index BAM
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FULL } from '../../../modules/local/samtools/main.nf'         // Index BAM
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_CHUNK } from '../../../modules/local/samtools/main.nf'         // Merge BAMs
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_FINAL } from '../../../modules/local/samtools/main.nf'         // Merge BAMs
include { SAMTOOLS_TOFASTQ       } from '../../../modules/local/samtools/main.nf'
include { CRAMINO_STATS          } from '../../../modules/local/cramino/main.nf'
include { SEQKIT_STATS           } from '../../../modules/local/seqkit/main.nf'          // Coverage stats
include { modifyMetaId           } from '../utils_nfcore_oncoseq_pipeline'
include { QUARTO_TABLE           } from '../../../modules/local/quarto/main.nf'           // Reporting (optional)
include { paramsSummaryMap       } from 'plugin/nf-schema'                                // Parameter summary
include { paramsSummaryMultiqc   } from '../../../subworkflows/nf-core/utils_nfcore_pipeline' // MultiQC summary
include { softwareVersionsToYAML } from '../../../subworkflows/nf-core/utils_nfcore_pipeline' // Version reporting
include { methodsDescriptionText } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline' // Methods for MultiQC

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


    main:
    ch_versions = Channel.empty() // For collecting version info

    ch_ref = ref
        .map { meta, genome, ref_fasta, _ref_fai ->
            tuple(meta, genome, ref_fasta) }

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
                ch_seqkit_out = Channel.empty()
            }
        } else {
            ch_seqkit_out = Channel.empty()
        }

    // Merge bams if multiple bam are provided when skip_mapping is used
    if (params.skip_mapping) {
        in_ch
            .flatMap { meta, bams ->
                // Ensure 'bams' is a list and flatten it
                def dir_list = bams instanceof List ? bams.flatten().sort() : [bams]
                def dir = file(dir_list[0])

                if (dir.isDirectory()) {
                    bams = dir.listFiles().findAll { f -> f.name ==~ /.*\.bam$/ }

                    if (bams.size() == 1) {
                        def bam_single = bams
                        def bai = dir.listFiles().findAll { f -> f.name ==~ /.*\.bai$/ }
                        def type = bai.size() > 0 ? 'single' : 'to_index'
                        return [tuple(meta, type, bams.flatten(), bai.flatten())]
                    }
                    else {
                        // Case 2: Multiple BAMs → split into chunks of 20 for merging
                        def counter = 0
                        return bams.collate(5).collect { chunk ->
                            counter += 1
                            def meta_chunk = modifyMetaId(meta, 'add_suffix', '', '', "_chunk${counter}")
                            tuple(meta_chunk, 'multi', chunk)
                        }
                    }
                }
                else {
                    def bai = "${bams}.bai"
                    def type = bai.size() > 0 ? 'single' : 'to_index'
                    def index = bai.size() > 0 ? bai : 'index'
                    def bam_file = bams
                    return [tuple(meta, type, bam_file, index)]   // Only bam file is provided
                }
            }
            .set { bam_chunks_ch }

            // Process only chunks (tuples with list of bams)
            bam_chunks_ch
                .branch { list ->
                    single: list[1] == 'single'
                    multi: list[1] == 'multi'
                    to_index: list[1] == 'to_index'
                }
                .set { bams_chunk_sep }

            SAMTOOLS_MERGE_CHUNK(bams_chunk_sep.multi.map{ meta_chunk, _type, chunk -> tuple(meta_chunk, chunk) } )

            SAMTOOLS_MERGE_CHUNK.out.bamfile
                .map { meta, bam ->
                    def meta_restore = modifyMetaId(meta, 'remove_suffix', '', '', /_chunk\d+$/)
                    tuple(meta_restore, bam)
                }
                .groupTuple()           // Group by original meta.id
                .set { final_merge_ch }

            // Count resulting bams
            ch_count = final_merge_ch
                .map { meta, bam_list -> tuple(meta, bam_list.size(), bam_list) }

            // Separate samples that require further merging from those that have only one resulting merged bam
            ch_count
                .branch { meta, bam_count, bam_list ->
                    single_bam: bam_count == 1
                    multi_bam: bam_count > 1
                }
                .set { ch_bam_merged }

            // Make another round of merging from those merged chunks
            SAMTOOLS_MERGE_FINAL(ch_bam_merged.multi_bam.map{ meta, size, bam_list -> tuple(meta, bam_list) })

            ch_to_index = SAMTOOLS_MERGE_FINAL.out.bamfile
                .mix(ch_bam_merged.single_bam.map{ meta, size, bam_list -> tuple(meta, bam_list) })
                .mix(bams_chunk_sep.to_index.map{ meta, _type, bam, bai -> tuple(meta, bam) })

            ch_single_bam = bams_chunk_sep.single
                .map { meta, _type, bam, bai ->
                tuple(meta, bam, bai) }

            SAMTOOLS_INDEX_FULL(ch_to_index)
            bam_ch = SAMTOOLS_INDEX_FULL.out.bamfile_index
                .mix(ch_single_bam)

            CRAMINO_STATS(bam_ch)

            ch_versions = SAMTOOLS_MERGE_CHUNK.out.versions
                .mix(SAMTOOLS_INDEX_FULL.out.versions)
                .mix(CRAMINO_STATS.out.versions)

        } else {

            // Prepare mapping input: clean up meta.id and join with reference
            ch_mapping_in = in_ch
                .join(ch_ref)
                .map { meta, fastq, genome, ref_fasta ->
                    def meta_ref = modifyMetaId(meta, 'add_suffix', '', '', "_${genome}")
                    tuple(meta_ref, fastq, ref_fasta)
                    }

            // Run minimap2 alignment
            MINIMAP2_ALIGN(ch_mapping_in)

            // Convert SAM to BAM
            SAMTOOLS_TOBAM(MINIMAP2_ALIGN.out.sam)
            // Sort and index BAM
            SAMTOOLS_SORT(SAMTOOLS_TOBAM.out.bamfile)
            SAMTOOLS_INDEX_FULL(SAMTOOLS_SORT.out.sortedbam)

            // Restore meta ID by removing ref id

            ch_ref_id = ch_ref
                .map { meta, genome, _ref_fasta ->
                    def meta_ref = modifyMetaId(meta, 'add_suffix', '', '', "_${genome}")
                    tuple(meta_ref, meta.id) }

            bam_ch = SAMTOOLS_INDEX_FULL.out.bamfile_index
                .join(ch_ref_id)
                .map { meta, bam, bai, meta_restore ->
                    tuple(id:meta_restore, bam, bai) }

            // Compute coverage stats
            CRAMINO_STATS(bam_ch)

            // Collect versions from all modules
            ch_versions = MINIMAP2_ALIGN.out.versions
                .mix(SAMTOOLS_TOBAM.out.versions)
                .mix(SAMTOOLS_SORT.out.versions)
                .mix(SAMTOOLS_INDEX_FULL.out.versions)
                .mix(CRAMINO_STATS.out.versions)
        }

    emit:
    bam      = bam_ch                                   // Final sorted BAM with index
    coverage = CRAMINO_STATS.out.stats                // Coverage stats
    seqkit   = ch_seqkit_out                          // Input fastq stats
    versions = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
