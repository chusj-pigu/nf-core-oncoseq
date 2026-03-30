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
        path(filt_vcf),
        path(gene_list)

    output:
    tuple val(meta),
        path("*filt.tsv"),
        emit: filt_tsv
    tuple val(meta),
        path("*ids.txt"),
        emit: filt_ids
    tuple val(meta),
        path("*region_fusions.txt"),
        emit: fusion_txt,
        optional:true
    tuple val(meta),
        path("*region_indel.txt"),
        emit: indel_txt,
        optional:true
    tuple val(meta),
        path("*table_fusions.tsv"),
        emit: fusion_tsv,
        optional:true
    tuple val(meta),
        path("*table_indel.tsv"),
        emit: indel_tsv,
        optional:true
    tuple val(meta),
        path("*targets_nohit.txt"),
        emit: targets,
        optional:true
    path "versions.yml",
        emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Process only high and moderate effect mutations
    grep -E \\
        'HIGH|MODERATE' \\
        "${filt_vcf}" > \\
        "${prefix}_filt.tsv" || :

    # Save their IDs for filtering (Sniffles2.* patterns)
    grep -oE \
        'Sniffles2\\.[A-Z]+\\.[A-Za-z0-9_]+' \\
        "${prefix}_filt.tsv" > \\
        "${prefix}_filt_ids.txt" || :

    # Ensure placeholders exist if empty
    [ -f "${prefix}_filt.tsv" ] || touch "${prefix}_filt.tsv"
    [ -f "${prefix}_filt_ids.txt" ] || touch "${prefix}_filt_ids.txt"

    # Transform into figeno region input file
    generate_sv_filt_regions.R \\
        --input "${prefix}_filt.tsv" \\
        --target ${gene_list}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1)

    END_VERSIONS
    """
}

process QDNASEQ_PROCESS {

    label 'local'
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'

    tag "$meta.id"

    input:
    tuple val(meta),
        path(calls_bed),
        path(segs_bed)

    output:
    tuple val(meta),
        path("*CNVs"),
        emit: cnv_file
    tuple val(meta),
        path("*ratio.txt"),
        emit: ratio_file

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk -F'\\t' '
        NR == 1 { next }
        {
            print \$1, \$2, \$3, (\$5+2), (\$5 > 0 ? "gain" : "loss")
        }
    ' OFS='\\t' $calls_bed > ${prefix}_CNVs

    awk -F'\\t' '
        BEGIN {
            OFS="\\t"
            print "Chromosome","Start","Ratio"
        }
        NR > 1 {
            print \$1, \$2, (\$5+1)
        }
    ' $segs_bed > ${prefix}_ratio.txt
    """
}

process ENSEMBL_VEP_TABLE {

    label 'local'
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'

    tag "$meta.id"

    input:
    tuple val(meta),
        path(bed)

    output:
    tuple val(meta),
        path("*.csv"),
        emit: csv

    when:
    bed.size() > 0

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk -F'\t' '
    BEGIN { OFS="," }
    {
        split(\$7, transcripts, ",")          # CSQ field
        split(\$9, ad_vals, ",")             # AD field (last column)
        
        ad_alt = (length(ad_vals) >= 2 ? ad_vals[2]+0 : 0)

        for (i in transcripts) {
            split(transcripts[i], f, "|")

            if (ad_alt > 5 && \$8 > 20 && \$6 == "PASS") {
                print \$1,\$2,\$3,\$4,\$5, \
                    f[2],f[3],f[4],f[7], \
                    f[11],f[12],f[18],f[49], \
                    \$8,ad_alt
            }
        }
    }' ${bed} > ${prefix}_filt.csv
    """
}
