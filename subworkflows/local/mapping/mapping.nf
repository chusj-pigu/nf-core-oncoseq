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
include { MINIMAP2_ALIGN         } from '../../../modules/local/minimap2/main.nf'         // minimap2 alignment
include { SAMTOOLS_TOBAM         } from '../../../modules/local/samtools/main.nf'         // Convert SAM to BAM
include { SAMTOOLS_SORT          } from '../../../modules/local/samtools/main.nf'         // Sort BAM
include { SAMTOOLS_INDEX         } from '../../../modules/local/samtools/main.nf'         // Index BAM
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
    //   fastq_ch: Channel of tuples [meta, reads] (reads can be file or directory)
    //   ref:      Channel of tuples [meta, ref, ref_fasta, ref_fai]
    take:
    fastq_ch  // Channel: from basecalling workflow or from --fastq if --skip_mapping is used
    ref       // Channel: from input samplesheet



    main:
    ch_versions = Channel.empty() // For collecting version info

    // Only expand fastq_ch if skip_basecalling is true
    if (params.skip_basecalling) {
        fastq_ch
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
            .set { fastq_ch }
    }


    // Prepare reference channel: extract meta and fasta path
    ch_ref = ref
        .map { meta, _ref, ref_fasta, _ref_fai ->
            tuple(meta, ref_fasta) }

    // Prepare mapping input: clean up meta.id and join with reference
    ch_mapping_in = fastq_ch
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
    // Compute coverage stats
    CRAMINO_STATS(SAMTOOLS_INDEX.out.bamfile_index)

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

    // Collect versions from all modules
    ch_versions = MINIMAP2_ALIGN.out.versions
        .mix(SAMTOOLS_TOBAM.out.versions)
        .mix(SAMTOOLS_SORT.out.versions)
        .mix(SAMTOOLS_INDEX.out.versions)
        .mix(CRAMINO_STATS.out.versions)

    emit:
    bam      = SAMTOOLS_INDEX.out.bamfile_index   // Final sorted BAM with index
    coverage = CRAMINO_STATS.out.stats                // Coverage stats
    versions = ch_versions                            // All tool versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
