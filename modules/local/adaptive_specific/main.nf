// Processes specific to adaptive sampling workflow that don't require containers to run

process REMOVE_PADDING {

    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'

    tag 'bed'

    input:
    tuple val(meta),
        path(bed)
        val(padding)

    output:
    tuple val(meta),
        path("*.bed"),
        emit: bed

    script:
    def prefix = task.ext.prefix ?: "${bed}.baseName"
    """
    awk -F'\t' \\
        'BEGIN {OFS="\t"} \\
        { $2 = $2 + ${padding}; $3 = $3 - ${padding}; print }' \\
        ${bed} > ${prefix}_nopadding.bed
    """
}
