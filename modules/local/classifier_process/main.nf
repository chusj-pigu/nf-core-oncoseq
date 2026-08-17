process PARSE_JSON_COMBINED {

    //TODO: SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/chusj-pigu/tidyverse:006154f90a8b1b1c8647a246a76cb0562517da61'
    label 'local'
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'

    tag "$meta.id"

    input:
    tuple val(meta),
        path(json)

    output:
    tuple val(meta),
        path("*alma.tsv"),
        emit: alma
    tuple val(meta),
        path("*crossnn_capper_et_al.tsv"),
        emit: capper
    tuple val(meta),
        path("*crossnn_pancan_devel_v5i.tsv"),
        emit: pancan
    tuple val(meta),
        path("*lamprey.tsv"),
        emit: lamprey
    tuple val(meta),
        path("*marlin.tsv"),
        emit: marlin
    tuple val(meta),
        path("*mpact.tsv"),
        emit: mpact
    tuple val(meta),
        path("*nanomix_llse.tsv"),
        emit: nanomix
    tuple val(meta),
        path("*sturgeon_brainstem.tsv"),
        emit: sturgeon_brainstem
    tuple val(meta),
        path("*sturgeon_general.tsv"),
        emit: sturgeon_general
    tuple val(meta),
        path("*tucan.tsv"),
        emit: tucan
    path "versions.yml",
        emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_classifier_json.R \\
        --sample ${prefix} \\
        --json ${json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1)
        dplyr: "\$(echo 'cat(as.character(packageVersion(\"dplyr\")))' | R --vanilla --slave)"
        jsonlite: "\$(echo 'cat(as.character(packageVersion(\"jsonlite\")))' | R --vanilla --slave)"
        optparse: "\$(echo 'cat(as.character(packageVersion(\"optparse\")))' | R --vanilla --slave)"
        purrr: "\$(echo 'cat(as.character(packageVersion(\"purrr\")))' | R --vanilla --slave)"
        readr: "\$(echo 'cat(as.character(packageVersion(\"readr\")))' | R --vanilla --slave)"
        stringr: "\$(echo 'cat(as.character(packageVersion(\"stringr\")))' | R --vanilla --slave)"

    END_VERSIONS
    """
}
