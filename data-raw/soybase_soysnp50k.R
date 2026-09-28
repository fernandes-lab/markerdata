#!/usr/bin/env Rscript
# data-raw/soybase_soysnp50k.R
#
# Build the *soybean* marker dataset (SoySNP50K_div) for the markerdata
# package, in the SAME numeric format as
# simplePHENOTYPES::SNP55K_maize282_maf04 (a data.frame whose first five
# columns are snp, allele, chr, pos, cm, followed by one -1/0/1 column per
# accession). markerdata does not ship this object; it publishes the packed
# result as a GitHub release asset that marker_data() downloads and caches
# on first use (see data-raw/release.R).
#
# Source: SoySNP50K iSelect BeadChip genotyping of the USDA Soybean Germplasm
# Collection -- 20,087 Glycine max / G. soja accessions x 42,509 SNPs.
#   Song Q. et al. (2013)  PLoS ONE 8:e54985                 (the array)
#   Song Q. & Hyten D.L. (2015) diversity dataset (SoyBase)  (these files)
# Files (Wm82.gnm2 assembly), verified reachable 2026-09-14:
#   https://data.soybase.org/Glycine/max/diversity/Wm82.gnm2.div.Song_Hyten_2015/
#     glyma.Wm82.gnm2.div.Song_Hyten_2015.hmp.gz   (HapMap, ~94 MB gz)
#     glyma.Wm82.gnm2.div.Song_Hyten_2015.vcf.gz   (VCF,    ~138 MB gz)
# License: public USDA / SoyBase data; please cite the papers above (see also
# data-raw/README.md -- redistribution terms for this dataset are not yet
# confirmed and should be checked before a CRAN release).
#
# The full HapMap is large, so this script (1) downloads it once and caches
# it under data-raw/soybase_cache/, (2) keeps a manageable slice of
# accessions (columns), (3) converts that slice to the numeric -1/0/1 format
# with simplePHENOTYPES::as_numeric(from = "hapmap"), (4) adds a synthetic
# genetic map (the SoyBase files carry no cM), (5) MAF-filters (>=0.05) and
# thins evenly across chromosomes to a curated n_snps-SNP subset, and (6)
# sanity-checks with as_population() and saves data-raw/built/SoySNP50K_div.rda.
#
# Usage:  Rscript data-raw/soybase_soysnp50k.R [n_accessions] [n_snps] [out.rda]
#   defaults: 300 accessions, 12000 SNPs -> data-raw/built/SoySNP50K_div.rda

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
n_acc  <- if (length(args) >= 1) as.integer(args[[1]]) else 300L
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/SoySNP50K_div.rda"

base <- "https://data.soybase.org/Glycine/max/diversity/Wm82.gnm2.div.Song_Hyten_2015"
url  <- file.path(base, "glyma.Wm82.gnm2.div.Song_Hyten_2015.hmp.gz")

cache <- "data-raw/soybase_cache"
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
gz <- file.path(cache, basename(url))
if (!file.exists(gz)) {
  message("Downloading SoySNP50K HapMap (~94 MB) ...")
  utils::download.file(url, gz, mode = "wb")
}

# Keep the 11 HapMap metadata columns + the first n_acc accession columns, so
# as_numeric() only converts a manageable slice. To pick specific PIs, replace the
# `keep` column indices; to convert ALL 20,087 accessions, use the .gz directly.
hmp_sub <- file.path(cache, sprintf("soysnp50k_%dacc.hmp.txt", n_acc))
if (!file.exists(hmp_sub)) {
  message("Subsetting to ", n_acc, " accessions ...")
  con_in <- gzfile(gz, "rt")
  hdr  <- strsplit(readLines(con_in, n = 1L), "\t", fixed = TRUE)[[1]]
  keep <- c(1:11, 11L + seq_len(min(n_acc, length(hdr) - 11L)))
  con_out <- file(hmp_sub, "wt")
  writeLines(paste(hdr[keep], collapse = "\t"), con_out)
  repeat {
    ln <- readLines(con_in, n = 20000L)
    if (!length(ln)) break
    f <- strsplit(ln, "\t", fixed = TRUE)
    writeLines(vapply(f, function(x) paste(x[keep], collapse = "\t"), character(1)),
               con_out)
  }
  close(con_in); close(con_out)
}

# HapMap (nucleotide/IUPAC codes) -> numeric (-1/0/1). as_numeric writes a sibling
# "<name>_numeric.txt" file (and returns NULL), matching the SP dataset's columns.
as_numeric(hmp_sub, from = "hapmap")
num_file <- sub("\\.hmp\\.txt$", "_numeric.txt", hmp_sub)
SoySNP50K_div <- utils::read.delim(num_file, check.names = FALSE,
                                   stringsAsFactors = FALSE)

# The SoyBase HapMap has no genetic map; build one from the physical positions so
# meiosis can be simulated (the engine errors clearly otherwise).
SoySNP50K_div$cm <- synthetic_map(SoySNP50K_div$chr, SoySNP50K_div$pos)

## finalize: MAF filter (>=0.05) + thin evenly across chromosomes to n_snps -
# This is the step that turns the full converted slice into the curated,
# shipped-size subset (mirrors the finalize step in the other data-raw/*.R
# scripts).
G   <- as.matrix(SoySNP50K_div[, 6:ncol(SoySNP50K_div), drop = FALSE])
p   <- (rowMeans(G) + 1) / 2
maf <- pmin(p, 1 - p)
keep <- which(maf >= 0.05)
message(sprintf("SNPs with MAF >= 0.05: %d / %d (chromosome labels: %d distinct).",
                length(keep), length(maf),
                length(unique(SoySNP50K_div$chr[keep]))))
ord <- keep[order(SoySNP50K_div$chr[keep], SoySNP50K_div$pos[keep])]
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]
SoySNP50K_div <- SoySNP50K_div[ord, , drop = FALSE]
rownames(SoySNP50K_div) <- NULL

## sanity check + save -------------------------------------------------------
stopifnot(all(unlist(SoySNP50K_div[, 6:ncol(SoySNP50K_div)]) %in% c(-1L, 0L, 1L)))
pop <- as_population(SoySNP50K_div,
                     individuals = seq_len(min(10L, ncol(SoySNP50K_div) - 5L)))
gt  <- table(unlist(SoySNP50K_div[, 6:ncol(SoySNP50K_div)]))
message(sprintf("Built SoySNP50K_div: %d SNPs x %d accessions (as_population OK: %s).",
                nrow(SoySNP50K_div), ncol(SoySNP50K_div) - 5L,
                inherits(pop, "Population")))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "))
save(SoySNP50K_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
