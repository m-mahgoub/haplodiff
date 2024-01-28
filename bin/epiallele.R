#!/usr/bin/env Rscript

library(epialleleR)
library(VariantAnnotation)
library(optparse)

source("~/.Rprofile")

# Define command line arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input BAM file path", metavar="file"),
  make_option(c("-g", "--genome"), type="character", default=NULL, help="Genome file path", metavar="file"),
  make_option(c("-v", "--vcf"), type="character", default=NULL, help="VCF file path", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output file path", metavar="file")
)

# Parse command line arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input)) {
  stop("Input BAM file is required", call. = FALSE)
}
if (is.null(opt$genome)) {
  stop("Genome file is required", call. = FALSE)
}
if (is.null(opt$vcf)) {
  stop("VCF file is required", call. = FALSE)
}
if (is.null(opt$output)) {
  stop("Output file is required", call. = FALSE)
}

# Assign arguments to variables
input.bam <- opt$input
genome <- opt$genome
vcf <- opt$vcf
output <- opt$output

# Proceed with the rest of your script
bam.data <- preprocessBam(input.bam, min.mapq = 1, skip.duplicates = TRUE, nthreads = 2, verbose = TRUE)

vcf.report <- generateVcfReport(
  bam=bam.data,
  vcf=vcf,
)

write.table(vcf.report, file=output, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
