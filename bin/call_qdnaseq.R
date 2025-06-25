#!/usr/bin/env Rscript

# Script to call CNVs with QDNAseq based on run_qdnaseq.R of wf-human-variation of epi2me https://github.com/epi2me-labs/wf-human-variation

suppressPackageStartupMessages({
    library(optparse)
    library(QDNAseq)
})

option_list <- list(
    make_option(c("-b", "--bam"),
        type = "character", default = NULL,
        help = "Path to the input bam file [default: %default]", metavar = "FILE"
    ),
    make_option(c("--prefix"),
        type = "character", default = NULL,
        help = "Output files prefix [default: %default]", metavar = "STRING"
    ),
    make_option(c("--binsize"),
        type = "integer", default = 500,
        help = "Bin size in Kbp [default: %default]", metavar = "NUMBER"
    ),
    make_option(c("-r", "--reference"),
        type = "character", default = NULL,
        help = "Reference type to use for GC/mappability bins, either hg38 or hg19 [default: %default]", metavar = "STRING"
    )
)

# ---- Parse options ----
opt <- parse_args(OptionParser(option_list = option_list))

if (opt$reference == "hg38" | opt$reference == "GRCh38") {
    library(QDNAseq.hg38)
    bins <- getBinAnnotations(binSize = opt$binsize, genome = "hg38")
} else if (opt$reference == "hg19" | opt$reference == "GRCh37") {
    library(QDNAseq.hg19)
    bins <- getBinAnnotations(binSize = opt$binsize, genome = "hg19")
}

## Primary analysis

read_counts <- binReadCounts(bins, bamfiles = opt$bam)

bed_raw <- paste(opt$prefix, "raw_bins.bed", sep = "_")
exportBins(read_counts, file = bed_raw, format = "bed", type = "copynumber")

# Apply filters and plot median read counts as a function of GC content and mappability
filt_read_counts <- applyFilters(read_counts, residual = TRUE, blacklist = TRUE)

# Estimate the GC / mappability correction
filt_read_counts <- estimateCorrection(filt_read_counts)

# Create copy numbers object
filt_read_counts <- applyFilters(filt_read_counts, chromosomes = NA)

# Apply correction for GC content and mappability
copy_number <- correctBins(filt_read_counts)

# Normalize to median
copy_number <- normalizeBins(copy_number)

# Smooth outliers
copy_number <- smoothOutlierBins(copy_number)

# Export bins
bed_norm <- paste(opt$prefix, "bins.bed", sep = "_")
exportBins(copy_number, bed_norm, format = "bed")

## Secondary analysis

cn_seg <- segmentBins(copy_number, transformFun = "sqrt")

cn_seg <- normalizeSegmentedBins(cn_seg)

cutoffDEL <- 0.5
cutoffLOSS <- 1.5
cutoffGAIN <- 2.5
cn_called <- callBins(cn_seg, method = "cutoff", cutoffs = log2(c(deletion = cutoffDEL, loss = cutoffLOSS, gain = cutoffGAIN, amplification = 10) / 2))



# Create PNG output
png_file <- paste(opt$prefix, "cov.png", sep = "_")
png(png_file)
# plot(cn_called)
plot(
  cn_called,
  delcol = "purple",
  losscol = "pink",
  gaincol = "lightblue",
  ampcol = "darkgreen",
  pointcol = "gray40",
  segcol = "orange"
)
dev.off()

noise_png_file <- paste(opt$prefix, "noise_plot.png", sep = "_")
png(noise_png_file)
noisePlot(filt_read_counts)
dev.off()

isobar_png_file <- paste(opt$prefix, "isobar_plot.png", sep = "_")
png(isobar_png_file)
isobarPlot(filt_read_counts)
dev.off()

# write outputs
bedout <- paste(opt$prefix, "calls.bed", sep = "_")
exportBins(cn_called, file = bedout, format = "bed", type = "calls")
vcfout <- paste(opt$prefix, "calls.vcf", sep = "_")
exportBins(cn_called, file = vcfout, format = "vcf", type = "calls")

segout <- paste(opt$prefix, "segs.seg", sep = "_")
exportBins(cn_called, file = segout, format = "seg", type = "segments")
segout <- paste(opt$prefix, "segs.bed", sep = "_")
exportBins(cn_called, file = segout, format = "bed", type = "segments")
segout <- paste(opt$prefix, "segs.vcf", sep = "_")
exportBins(cn_called, file = segout, format = "vcf", type = "segments")
