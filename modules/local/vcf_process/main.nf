process SV_PROCESS {

    //TODO: SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/chusj-pigu/tidyverse:latest'
    label 'local'
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'

    tag "$meta.id"

    input:
    tuple val(meta),
        path(filt_vcf)

    output:
    tuple val(meta),
        path("*.tsv"),
        emit: filt_tsv
    tuple val(meta),
        path("region_fusions.txt"),
        emit: fusion_txt,
        optional:true
    tuple val(meta),
        path("region_indel.txt"),
        emit: indel_txt,
        optional:true
    path "versions.yml",
        emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep -E 'HIGH|MODERATE' ${filt_vcf} > ${prefix}_sv_filt.tsv
    generate_sv_filt_regions.R \\
        --input ${prefix}_sv_filt.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1)

    END_VERSIONS
    """
}
