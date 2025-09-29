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
include { SAMTOOLS_SORT                          } from '../../../modules/local/samtools/main.nf'         // Sort BAM
include { SAMTOOLS_INDEX                         } from '../../../modules/local/samtools/main.nf'         // Index BAM
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_CHUNK } from '../../../modules/local/samtools/main.nf'         // Merge BAMs
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_FINAL } from '../../../modules/local/samtools/main.nf'         // Merge BAMs
include { CRAMINO_STATS          } from '../../../modules/local/cramino/main.nf'          // Coverage stats
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
        }

    // Merge bams if multiple bam are provided when skip_mapping is used
    if (params.skip_mapping) {
        in_ch
            .flatMap { meta, bams ->
                // Ensure 'bams' is a list and flatten it
                def dir_list = bams instanceof List ? bams.flatten().sort() : [bams]
                def dir = file(dir_list[0])

                bams = dir.listFiles().findAll { f -> f.name ==~ /.*\.bam$/ }
                if (bams.size() == 1) {
                    def bam_single = bams
                    def bai = dir.listFiles().findAll { f -> f.name ==~ /.*\.bai$/ }
                    return [tuple(meta, 'single', bams.flatten(), bai.flatten())]
                } else {
                    // Case 2: Multiple BAMs → split into chunks of 20 for merging
                    def counter = 0
                    return bams.collate(20).collect { chunk ->
                        counter++
                        def meta_chunk = modifyMetaId(meta, 'add_suffix', '', '', "_chunk${counter}")
                        tuple(meta_chunk, 'multi', chunk)
                    }
                }
            }
            .set { bam_chunks_ch }

            // Process only chunks (tuples with list of bams)
            bam_chunks_ch
                .branch { list ->
                    single: list[1] == 'single'
                    multi: list[1] == 'multi'
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

            SAMTOOLS_INDEX(ch_to_index)
            bam_ch = SAMTOOLS_INDEX.out.bamfile_index
                .mix(bams_chunk_sep.single.map{ meta, _type, bam, bai -> tuple(meta, bam, bai)})
            CRAMINO_STATS(bam_ch)

            ch_versions = SAMTOOLS_MERGE_CHUNK.out.versions
                .mix(SAMTOOLS_INDEX.out.versions)
                .mix(CRAMINO_STATS.out.versions)

        } else {
            // Prepare reference channel: extract meta and fasta path
            ch_ref = ref
                .map { meta, _ref, ref_fasta, _ref_fai ->
                    tuple(meta, ref_fasta) }

            // Prepare mapping input: clean up meta.id and join with reference
            ch_mapping_in = in_ch
                .map { meta, reads ->
                    def new_meta = modifyMetaId(meta, 'remove_suffix', '', '', '_pass')
                    tuple(new_meta, reads)
                }
                .join(ch_ref)

            // Run minimap2 alignment
            MINIMAP2_ALIGN(ch_mapping_in)

            // Convert SAM to BAM
            SAMTOOLS_TOBAM(MINIMAP2_ALIGN.out.sam)
            // Sort and index BAM
            SAMTOOLS_SORT(SAMTOOLS_TOBAM.out.bamfile)
            SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sortedbam)
            bam_ch = SAMTOOLS_INDEX.out.bamfile_index

            // Compute coverage stats
            CRAMINO_STATS(bam_ch)

            // Collect versions from all modules
            ch_versions = MINIMAP2_ALIGN.out.versions
                .mix(SAMTOOLS_TOBAM.out.versions)
                .mix(SAMTOOLS_SORT.out.versions)
                .mix(SAMTOOLS_INDEX.out.versions)
                .mix(CRAMINO_STATS.out.versions)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        (Optional) Example for coverage reporting with QUARTO_TABLE
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ch_mosdepth_coverage = MOSDEPTH_GENERAL.out.summary
        .map { meta, table ->
            def lines = table.readLines()
            def coverage = lines[-1].tokenize('\t')[3].toDouble()
            tuple(meta, coverage)
        }
        .collectFile { item ->
            def sample_id = item[0].id
            [ "coverage.txt", sample_id + '\t' + item[1] + '\n']
        }
    QUARTO_TABLE( ch_mosdepth_coverage_table, "Mean coverage", "T", "Mosdepth_coverage", "mosdepth-general")
    */

    emit:
    bam      = bam_ch                                   // Final sorted BAM with index
    coverage = CRAMINO_STATS.out.stats                // Coverage stats
    versions = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
