process SV_PROCESS {

    //TODO: SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/chusj-pigu/tidyverse:006154f90a8b1b1c8647a246a76cb0562517da61'
    label 'local'
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'

    tag "$meta.id"

    input:
    tuple val(meta),
        path(sniffles_vcf),
        path(severus_vcf),
        path(bed),
        path(gene_list),
        path(blacklist)

    output:
    tuple val(meta),
        path("*filt.tsv"),
        emit: filt_tsv
    tuple val(meta),
        path("*region_fusions.txt"),
        emit: fusion_txt,
        optional:true
    tuple val(meta),
        path("*region_other.txt"),
        emit: other_txt,
        optional:true
    tuple val(meta),
        path("*table_fusions.tsv"),
        emit: fusion_tsv,
        optional:true
    tuple val(meta),
        path("*table_other.tsv"),
        emit: other_tsv,
        optional:true
    tuple val(meta),
        path("*targets_nohit.txt"),
        emit: targets,
        optional:true
    tuple val(meta),
        path("*table_figeno.tsv"),
        emit: figeno_table,
        optional:true
    path "versions.yml",
        emit: versions

    script:
    def prefix_sniffles = task.ext.prefix ?: "${meta.id}_sniffles"
    def prefix_severus = task.ext.prefix ?: "${meta.id}_severus"
    """
    # Process only high and moderate effect mutations
    grep -v '^##' \\
        "${sniffles_vcf}" > \\
        "${prefix_sniffles}_filt.tsv" || :

    grep -v '^##' \\
        "${severus_vcf}" > \\
        "${prefix_severus}_filt.tsv" || :

    # Ensure placeholders exist if empty
    [ -f "${prefix_sniffles}_filt.tsv" ] || touch "${prefix_sniffles}_filt.tsv"
    [ -f "${prefix_severus}_filt.tsv" ] || touch "${prefix_severus}_filt.tsv"

    # Transform into figeno region input file
    generate_sv_filt_regions.R \\
        --target ${gene_list} \\
        --exclude ${blacklist} \\
        --panel ${bed}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1)

    END_VERSIONS
    """
}

process STELLERATOR_PROCESS {

    //TODO: SET FIXED VERSION WHEN PIPELINE IS STABLE
    container 'ghcr.io/chusj-pigu/tidyverse:006154f90a8b1b1c8647a246a76cb0562517da61'
    label 'local'
    label 'process_low'
    label 'process_single_cpu'
    label 'process_very_low_memory'

    tag "$meta.id"

    input:
    tuple val(meta),
        path(tables)

    output:
    tuple val(meta),
        path("*table_fusions.tsv"),
        emit: tsv
    tuple val(meta),
        path("*fusions.txt"),
        emit: figeno
    path "versions.yml",
        emit: versions

    script:
    """
    generate_sv_filt_regions_stellerator.R

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
        path(bed),
        path(list_exclude)

    output:
    tuple val(meta),
        path("*.csv"),
        emit: csv

    when:
    bed.size() > 0

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_depth_threshold = params.cfdna ? '' : "&& \$8 >= 20" // Only apply read depth filter for non-cfdna samples
    def pass_filter = prefix.contains('germline') ? '&& \$6 == "PASS"' : ''
    """
        awk -F'\t' '
        BEGIN {
            OFS=","
            print "CHROM,POS,REF,ALT,QUAL,SYMBOL,Consequence,IMPACT,CLIN_SIG,Feature,RefSeq_ID,HGVSc,HGVSp,Existing_variation,SIFT,PolyPhen,gnomADe_AF,Read_depth (DP),Variant_depth (AD)"
        }
        {
            split(\$7, transcripts, ",")
            split(\$9, ad_vals, ",")
            ad_alt = (length(ad_vals) >= 2 ? ad_vals[2]+0 : 0)
            for (i in transcripts) {
                split(transcripts[i], f, "|")
                if (ad_alt >= 5 ${read_depth_threshold} ${pass_filter}) {
                    print \$1,\$2,\$3,\$4,\$5, \
                    f[4],f[2],f[3],f[72],f[7],f[27], \
                    f[11],f[12],f[18],f[37],f[38],f[49], \
                    \$8,ad_alt
                }
            }
        }' ${bed} | grep -vF -f ${list_exclude} > ${prefix}_filt.csv
    """
}
