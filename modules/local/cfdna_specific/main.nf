process CAT_FASTQ {

    //TODO: FIND A BETTER SOLUTION FOR THIS
    label 'local'
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'
    label 'process_very_low_time'

    tag "$meta.id"

    input:
    tuple val(meta),
        path(fastq_dir)

    output:
    tuple val(meta),
        path("*.fq.gz"),
        emit: merged_fq

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${fastq_dir} > ${prefix}.fq.gz
    """
}

process PIGZ_FA {

    //TODO: FIND A BETTER SOLUTION FOR THIS
    container 'ghcr.io/chusj-pigu/samtools:latest'
    label 'process_low'
    label 'process_low_cpu'
    label 'process_very_low_memory'

    tag "$meta.id"

    input:
    tuple val(meta),
        path(fa)

    output:
    tuple val(meta),
        path("*fa.gz"),
        emit: fa

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus
    """
    pigz -c \\
        -p ${threads} \\
        ${fa} > ${prefix}.fa.gz
    """
}
