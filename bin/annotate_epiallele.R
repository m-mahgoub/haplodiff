#!/usr/bin/env Rscript

library(optparse)
library(annotatr)
library(GenomicRanges)
library(VariantAnnotation)
library(rtracklayer)

source("~/.Rprofile")

# Setting up command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to the input file (epiallele report)", metavar="file"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to the output file (annotated report)", metavar="file")
)

# Parsing command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$input)) {
  stop("Error: An input file path is required", call. = FALSE)
}

if (is.null(opt$output)) {
  stop("Error: An output file path is required", call. = FALSE)
}

# Use the provided input and output paths
epiallele_report_path <- opt$input
annotated_report_output <- opt$output

# Read the file (assuming it's tab-separated and has a header)
epiallele_report <- read.table(epiallele_report_path, header = TRUE, sep = "\t")

# Data manipulation steps
epiallele_report$chromStart <- epiallele_report$range - 1
names(epiallele_report)[names(epiallele_report) == "range"] <- "chromEnd"
bed_data <- epiallele_report[, c("seqnames", "chromStart", "chromEnd", "name")]
write.table(bed_data, file = "meth.tmp.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Read regions and annotate
regions <- read_regions("meth.tmp.bed", genome = 'hg38', format = 'bed')
annots <- c("hg38_genes_promoters", "hg38_genes_1to5kb", "hg38_genes_exons", "hg38_genes_introns")
annotations <- build_annotations(genome = 'hg38', annotations = annots)
report_annotated <- annotate_regions(regions = regions, annotations = annotations, ignore.strand = TRUE, quiet = FALSE)

# Write output
df_report_annotated <- data.frame(report_annotated)
write.table(df_report_annotated, annotated_report_output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
