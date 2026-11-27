#!/usr/bin/env bash

# Copyright (c) the nf-core community under an open-source MIT license.
#
# Build a chromosome name mapping (original -> UCSC style) from a SAM/BAM header.
#
# Usage: samtools view -H in.bam | chrom_name_map.sh > chrom_map.tsv
#        chrom_name_map.sh < sample_sequences.txt > chrom_map.tsv   (SAMTOOLS_SPLITHEADER @SQ output)
#
# Output: tab-separated, two columns (original_name, ucsc_name), one line per @SQ.
#         A leading "# convention: <X>" comment line reports the detected naming
#         convention: UCSC, NCBI_GRCh (NC_0000NN), NCBI_T2T (NC_0609NN),
#         ENSEMBL (1..22, X, Y, MT), MIXED or UNKNOWN.
#
# Recognised names:
#   NC_000001.x .. NC_000022.x -> chr1..chr22, NC_000023.x -> chrX, NC_000024.x -> chrY
#   NC_060925.x .. NC_060946.x -> chr1..chr22, NC_060947.x -> chrX, NC_060948.x -> chrY (T2T-CHM13v2.0)
#   NC_012920.x                -> chrM
#   1..22, X, Y, M, MT         -> chr1..chr22, chrX, chrY, chrM
#   chr*                       -> unchanged
# Anything else (unplaced scaffolds, decoys...) is mapped to itself.

set -euo pipefail

awk -F'\t' '
function sex(n) { return n == 23 ? "chrX" : (n == 24 ? "chrY" : "chr" n) }
$1 == "@SQ" {
    sn = ""
    for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) { sn = substr($i, 4); break }
    if (sn == "") next

    base = sn; sub(/\.[0-9]+$/, "", base)    # drop accession version
    ucsc = sn; style = ""

    if (base ~ /^NC_0000(0[1-9]|1[0-9]|2[0-4])$/) {
        ucsc = sex(substr(base, 8, 2) + 0); style = "NCBI_GRCh"
    } else if (base ~ /^NC_0609(2[5-9]|3[0-9]|4[0-8])$/) {
        ucsc = sex(substr(base, 4, 6) - 60924); style = "NCBI_T2T"
    } else if (base == "NC_012920") {
        ucsc = "chrM"
    } else if (sn ~ /^chr/) {
        style = "UCSC"
    } else if (sn ~ /^([0-9]|1[0-9]|2[0-2]|X|Y|MT|M)$/) {
        ucsc = (sn == "MT" ? "chrM" : "chr" sn); style = "ENSEMBL"
    }

    if (style != "") seen[style] = 1
    names[++n] = sn; map[n] = ucsc
}
END {
    k = 0; conv = "UNKNOWN"
    for (s in seen) { k++; conv = s }
    if (k > 1) conv = "MIXED"
    print "# convention: " conv
    for (i = 1; i <= n; i++) print names[i] "\t" map[i]
}
'
