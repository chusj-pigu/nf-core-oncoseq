/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_TEXT } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE } from '../../../modules/local/quarto/main.nf'
include { modifyMetaId } from '../../../subworkflows/local/utils_nfcore_oncoseq_pipeline'


workflow WGS_REPORT {

    take:
    ch_cramino
    ch_seqkit
    ch_minqs

    main:

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SECTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_section_inputs = Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSTRUCT QUARTO TABLE FOR CRAMINO STATS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_cramino_stats = ch_cramino
        .map { meta, table ->
            def lines = table.readLines()
            def nreads = (lines[3].tokenize('\t')[1].toDouble() / 1000000).toString()
            def nbases = lines[4].tokenize('\t')[1]
            def cov = lines[5].tokenize('\t')[1]
            def n50 = lines[7].tokenize('\t')[1]
            tuple(meta.id, nreads, nbases, cov, n50)
        }

    if (params.skip_mapping) {
        ch_stats = ch_cramino_stats
            .collectFile { table ->
                def content = ["Metric" + '\t' + "Value" + '\n' +
                    "Mean Coverage" + '\t' + table[3] + 'X' + '\n' +
                    "N reads (M)" + '\t' + table[1] + '\n' +
                    "N bases (GB)" + '\t' + table[2] + '\n' +
                    "N50" + '\t' + table[4]].join('\n')
                return [ "${table[0]}_wgs_stats.tsv", content + '\n' ]
            }
            .map { table ->
                def meta = table.name.replace('_wgs_stats.tsv', '')
                tuple(id:meta,table)
            }
            .combine(ch_minqs)
            .map { tuple ->
            // Extract the existing values from the tuple
            def (meta, file, minqs) = tuple

            // Transform chrom into two new variables
            def caption = "Summary statistics cramino for reads with phred QS > ${minqs}"
            def col_names = ""
            def section = "QC"
            def process = "stats-qc-${meta.id}"

            // Return a new tuple with the additional variables
            return [meta, file, caption, col_names, section, process ]
        }

    } else {
        ch_compared_stats = ch_seqkit
            .map { meta, table ->
                //Read the file content as a list of lines
                def lines = table.readLines()
                // Second line and 4th column, in million of reads
                def n_reads = (lines[1].tokenize('\t')[3].toDouble() / 1000000) .toString()
                // In billion of bases
                def n_bases = (lines[1].tokenize('\t')[4].toDouble() / 1000000000) .toString()
                def n50 =  lines[1].tokenize('\t')[12]
            tuple(meta.id,n_reads,n_bases,n50)
            }
            .join(ch_cramino_stats)
            .collectFile { table ->
                def content = ["Metric" + '\t' + "Pre-alignment" + '\t' + "Post alignment" + '\n' +
                    "Mean Coverage" + '\t' + "" + '\t' + table[6] + 'X' + '\n' +
                    "N reads (M)" + '\t' + table[1] + '\t' + table[4] + '\n' +
                    "N bases (GB)" + '\t' + table[2] + '\t' + table[5] + '\n' +
                    "N50" + '\t' + table[3] + '\t' + table[7]].join('\n')
                return [ "${table[0]}_wgs_stats.tsv", content + '\n' ]
            }
            .map { table ->
                def meta = table.name.replace('_wgs_stats.tsv', '')
                tuple(id:meta,table)
            }
            .combine(ch_minqs)

        ch_stats = ch_compared_stats
            .map { tuple ->
            // Extract the existing values from the tuple
            def (meta, file, minqs) = tuple

            // Transform chrom into two new variables
            def caption = "Summary statistics with seqkit and cramino for reads with phred QS > ${minqs}"
            def col_names = ""
            def section = "QC"
            def process = "stats-qc-${meta.id}"

            // Return a new tuple with the additional variables
            return [meta, file, caption, col_names, section, process ]
        }
    }

    QUARTO_TABLE(
        ch_stats
    )


    ch_section_inputs = QUARTO_TABLE.out.quarto_table


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSTRUCT QUARTO FIGURE FOR MPGI COVERAGE PLOT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


    ch_section_inputs = ch_section_inputs
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    QUARTO_SECTION(
        ch_section_inputs,
        "Summary statistics"
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    emit:
    sections = QUARTO_SECTION.out.quarto_section

}
