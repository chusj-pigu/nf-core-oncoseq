process STAGE_REFERENCE_FILES {
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'
    label 'process_very_low_time'

    tag "$label"

    input:
    tuple val(meta),
        val(label),
        path(ref_files)

    output:
    tuple val(meta),
        val(label),
        path('staged/*'),
        emit: staged
    path 'versions.yml',
        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cacheDir = params.reference_cache_dir
        ? params.reference_cache_dir.toString()
        : null
    def cacheLinks = cacheDir
        ? """
          mkdir -p '${cacheDir}'
          for file in staged/*; do
              ln -sfn "\$(readlink -f "\${file}")" '${cacheDir}/'\$(basename "\${file}")
          done
          """
        : ''

    """
    mkdir staged

    for file in ${ref_files}; do
        ln -sfn "\$(readlink -f "\${file}")" staged/\$(basename "\${file}")
    done

    ${cacheLinks}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staged: "true"
    END_VERSIONS
    """
}
