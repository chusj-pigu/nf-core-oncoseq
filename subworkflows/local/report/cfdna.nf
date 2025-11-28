/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_SECTION as QUARTO_SECTION_QC  } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION as QUARTO_SECTION_CNV } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_QC    } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE as QUARTO_FIGURE_CNV   } from '../../../modules/local/quarto/main.nf'
include { CONVERT_PDF_PNG                      } from '../../../modules/local/magick/main.nf'
include { QUARTO_TABLE_COLNAMES                } from '../../../modules/local/quarto/main.nf'


workflow CFNDA_REPORT {

    take:
    ch_id
    ch_seqkit
    ch_cramino
    ch_ichor_fig

    main:

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSTRUCT QUARTO TABLE FOR CRAMINO and NANOPLOT STATS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_id_table = ch_id
        .map { meta, project ->
        tuple(meta.id, project) }
    ch_cramino_txt = ch_cramino
        .map { meta, table ->
            def lines = table.readLines()
            def nreads = (lines[3].tokenize('\t')[1].toDouble() / 1000000).toString()
            def nbases = lines[4].tokenize('\t')[1]
            def cov = lines[5].tokenize('\t')[1]
            def n50 = lines[7].tokenize('\t')[1]
            tuple(meta.id, nreads, nbases, cov, n50)
        }

    ch_stats = ch_seqkit
        .map { meta, table ->
            //Read the file content as a list of lines
            def lines = table.readLines()
            // Second line and 4th column, in million of reads
            def n_reads = (lines[1].tokenize('\t')[3].toDouble() / 1000000) .toString()
            // In billion of bases
            def n_bases = (lines[1].tokenize('\t')[4].toDouble() / 1000000000) .toString()
            def n50 =  lines[1].tokenize('\t')[12]
            tuple(meta.id, n_reads, n_bases, n50) }
        .join(ch_cramino_txt)
        .join(ch_id_table)
        .map { meta, nreads_filt, nbases_filt, n50_filt, nreads_al, nbases_al, cov, n50_al, project ->
            tuple(project, meta, nreads_filt, nreads_al, nbases_filt, nbases_al, cov, n50_filt, n50_al) }
        .collectFile { table ->
            def content = [table[1] + '\t' + table[2] + '\t' + table[3] +
             '\t' + table[4] + '\t' + table[5] + '\t' + table[6] + '\t' +
             table[7] + '\t' + table[8] ].join('\n')
            return [ "${table[0]}_cfdna_stats.tsv", content + '\n' ]
        }
        .map { table ->
            def meta = table.name.replace('_cfdna_stats.tsv', '')
            tuple(id:meta,table)
            }                                                             // Add back meta to the new table

    ch_quarto_table = ch_stats
        .map { meta, table ->
            def caption = "Stats of filtered and aligned reads"
            def col_names = "Sample, N reads filtered (M), N reads aligned (M), N bases filtered (GB), N bases faligned (GB), Coverage (X), N50 reads filtered, N50 reads aligned"
            def section = "QC"
            def process = "stats_qc-${meta.id}"
            tuple(meta, table, caption, col_names, section, process)
        }

    QUARTO_TABLE_COLNAMES(
        ch_quarto_table
    )


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSTRUCT QUARTO FIGURE FOR MPGI COVERAGE PLOT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // ch_qc_files = ch_nanoplot_fig
    //     .map { meta, fig ->
    //         def filtered = fig.find { it.name.endsWith('WeightedLogTransformed_HistogramReadlength.svg') }
    //         tuple(meta, filtered)
    //     }
    //     .map { tuple ->
    //         // Extract the existing values from the tuple
    //         def (meta, file) = tuple

    //         // Transform chrom into two new variables
    //         def caption = "Read length distribution for ${meta.id}"
    //         def section = "QC"
    //         def process = "length-plot-${meta.id}"

    //         // Return a new tuple with the additional variables
    //         return [meta, file, caption, section, process ]
    //     }


    // QUARTO_FIGURE_QC(
    //     ch_qc_files
    // )

    ch_section_qc = QUARTO_TABLE_COLNAMES.out.quarto_table
    //     .mix(QUARTO_FIGURE_QC.out.quarto_figure)

    ch_section_qc = ch_section_qc
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    QUARTO_SECTION_QC(
        ch_section_qc,
        "QC statistics"
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONSTRUCT QUARTO FIGURE AND SECTION FOR ICHORCNA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    CONVERT_PDF_PNG(ch_ichor_fig)

    ch_ichor_files = CONVERT_PDF_PNG.out.png
        .join(ch_id)
        .map { meta, file, project ->
            // Transform chrom into two new variables
            def caption = "IchorCNA best solution for ${meta.id}"
            def section = "CNV"
            def process = "ichor-plot-${meta.id}"

            // Return a new tuple with the additional variables
            tuple(id:project, file, caption, section, process )
        }

    QUARTO_FIGURE_CNV(
        ch_ichor_files
    )

    ch_section_cnv = QUARTO_FIGURE_CNV.out.quarto_figure
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    QUARTO_SECTION_CNV(
        ch_section_cnv,
        "CNV"
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    ch_sections = QUARTO_SECTION_QC.out.quarto_section
        .mix(QUARTO_SECTION_CNV.out.quarto_section)

    emit:
    sections = ch_sections

}
