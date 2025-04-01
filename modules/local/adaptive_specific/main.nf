// Processes specific to adaptive sampling workflow that don't require containers to run

process REMOVE_PADDING {

    label 'local'
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'

    tag 'bed'

    input:
    tuple val(meta),
        path(bed),
        val(padding)

    output:
    tuple val(meta),
        path("*.bed"),
        emit: bed

    script:
    def prefix = task.ext.prefix ?: "${bed}.baseName"
    def output_name = "${prefix}_nopadding.bed"
    """
    awk -F'\\t' '
        BEGIN {OFS="\\t"}
        { \$2 = \$2 + $padding; \$3 = \$3 - $padding; print }
    ' $bed > $output_name
    """
}

process PIGZ_BED {

    //TODO: FIND A BETTER SOLUTION FOR THIS
    container 'ghcr.io/chusj-pigu/samtools:latest'
    label 'local'
    label 'process_low'
    label 'process_low_cpu'
    label 'process_very_low_memory'

    tag '$meta.id'

    input:
    tuple val(meta),
        path(bed)

    output:
    tuple val(meta),
        path("*.bed"),
        emit: bed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus
    """
    pigz -dc \\
        -p ${threads} \\
        ${bed} > ${prefix}.bed
    """
}

process COVERAGE_PLOT {

    //TODO: SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/chusj-pigu/tidyverse:latest'
    label 'local'
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'

    tag 'bed'

    input:
    tuple val(meta),
        path(nofilt_bed),
        path(primary_bed),
        path(unique_bed),
        val(bg_cov),
        path(low_fidelity_list)

    output:
    tuple val(meta),
        path("*.pdf"),
        emit: cov_plot

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    coverage_plot.R \\
        -n ${nofilt_bed} \\
        -p ${primary_bed} \\
        -u ${unique_bed} \\
        -b ${bg_cov} \\
        -l ${low_fidelity_list} \\
        -o ${prefix}_coverage_mapq.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1)

    END_VERSIONS
    """
}
