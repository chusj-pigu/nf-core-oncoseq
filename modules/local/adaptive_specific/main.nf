process COVERAGE_PLOT {

    //TODO: SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/chusj-pigu/tidyverse:006154f90a8b1b1c8647a246a76cb0562517da61'
    label 'local'
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'

    tag "$meta.id"

    input:
    tuple val(meta),
        path(nofilt_bed),
        path(primary_bed),
        path(unique_bed),
        val(bg_cov),
        path(low_fidelity_list),
        path(header_sq)     // @SQ lines of the BAM header (SAMTOOLS_SPLITHEADER sequences)

    output:
    tuple val(meta),
        path("*.pdf"),
        emit: cov_plot_pdf
    tuple val(meta),
        path("*.svg"),
        emit: cov_plot_svg
    tuple val(meta),
        path("*.csv"),
        emit: cov_df
    tuple val(meta),
        path("*_chrom_map.tsv"),
        emit: chrom_map
    path "versions.yml",
        emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Map BAM chromosome names (UCSC / NCBI GRCh / NCBI T2T / Ensembl) to UCSC names
    chrom_name_map.sh < ${header_sq} > ${prefix}_chrom_map.tsv

    coverage_plot.R \\
        -n ${nofilt_bed} \\
        -p ${primary_bed} \\
        -u ${unique_bed} \\
        -b ${bg_cov} \\
        -l ${low_fidelity_list} \\
        -c ${prefix}_chrom_map.tsv \\
        -o ${prefix}_coverage_mapq.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1)

    END_VERSIONS
    """
}
