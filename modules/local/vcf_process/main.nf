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
        path(stellerator_vcf),
        path(bed),
        path(gene_list),
        path(blacklist)

    output:
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
    """
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
    def support = params.cfdna
        ? 2
        : (params.realtime != null && params.realtime <= 6 ? 0 : 3)
    """
    generate_sv_filt_regions_stellerator.R \\
        -m ${support}

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
    def read_depth_threshold = params.cfdna ? '' : "&& dp >= 20" // Only apply read depth filter for non-cfdna samples
    def pass_filter = prefix.contains('germline') ? '&& \$6 == "PASS"' : ''   // Keep Nonsomatic variants for somatic vcf
    """
    awk -F'\t' '
        BEGIN {
            OFS=","
            print "CHROM,POS,REF,ALT,QUAL,FILTER,SYMBOL,Consequence,IMPACT,CLIN_SIG,Feature,RefSeq_ID,HGVSc,HGVSp,Existing_variation,SIFT,PolyPhen,gnomADe_AF,Read_depth (DP),Variant_depth (AD),VAF (%)"
        }
        /^#/ { next }
        {
            dp = \$8+0
            split(\$9, ad_vals, ",")
            split(\$4, alts, ",")
            split(\$7, transcripts, ",")

            for (i = 1; i <= length(transcripts); i++) {
                split(transcripts[i], f, "|")

                if (f[1] != "" && f[1] != "-") {
                    # Allele is explicit: single match by value
                    alt_allele = f[1]
                    alt_idx = 0
                    for (a = 1; a <= length(alts); a++) {
                        if (alts[a] == alt_allele) { alt_idx = a; break }
                    }
                    ad_alt = (alt_idx > 0 ? ad_vals[alt_idx + 1]+0 : 0)
                    vaf = sprintf("%.2f", ad_alt / dp * 100)
                    if (ad_alt >= 5 ${read_depth_threshold} ${pass_filter}) {
                        print \$1,\$2,\$3,alt_allele,\$5,\$6,
                            f[4],f[2],f[3],f[72],f[7],f[27],
                            f[11],f[12],f[18],f[38],f[39],f[49],
                            dp,ad_alt,vaf
                    }
                } else {
                    # Allele is "-": VEP collapsed all ALTs into one annotation,
                    # so emit one row per ALT each with its own AD depth
                    for (a = 1; a <= length(alts); a++) {
                        ad_alt = ad_vals[a + 1]+0
                        vaf = sprintf("%.2f", ad_alt / dp * 100)
                        if (ad_alt >= 5 ${read_depth_threshold} ${pass_filter}) {
                            print \$1,\$2,\$3,alts[a],\$5,\$6,
                                f[4],f[2],f[3],f[72],f[7],f[27],
                                f[11],f[12],f[18],f[38],f[39],f[49],
                                dp,ad_alt,vaf
                        }
                    }
                }
            }
        }
    ' ${bed} | grep -vF -f ${list_exclude} > ${prefix}_filt.csv
    """
}
