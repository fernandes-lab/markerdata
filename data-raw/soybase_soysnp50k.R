#!/usr/bin/env Rscript
# data-raw/soybase_soysnp50k.R
#
# Build a *soybean* genotype dataset for the Breeding Designer / simplePHENOTYPES,
# in the SAME numeric format as simplePHENOTYPES::SNP55K_maize282_maf04 (a
# data.frame whose first five columns are snp, allele, chr, pos, cm, followed by
# one -1/0/1 column per accession).
#
# Source: SoySNP50K iSelect BeadChip genotyping of the USDA Soybean Germplasm
# Collection -- 20,087 Glycine max / G. soja accessions x 42,509 SNPs.
#   Song Q. et al. (2013)  PLoS ONE 8:e54985                 (the array)
#   Song Q. & Hyten D.L. (2015) diversity dataset (SoyBase)  (these files)
# Files (Wm82.gnm2 assembly), verified reachable 2026-09-14:
#   https://data.soybase.org/Glycine/max/diversity/Wm82.gnm2.div.Song_Hyten_2015/
#     glyma.Wm82.gnm2.div.Song_Hyten_2015.hmp.gz   (HapMap, ~94 MB gz)
#     glyma.Wm82.gnm2.div.Song_Hyten_2015.vcf.gz   (VCF,    ~138 MB gz)
# License: public USDA / SoyBase data; please cite the papers above.
#
# The full HapMap is large, so this script keeps a subset of accessions (columns)
# to make a tractable example, converts it to the numeric -1/0/1 format with
# simplePHENOTYPES::as_numeric(from = "hapmap"), adds a synthetic genetic map (the
# SoyBase HapMap carries no cM), sanity-checks with as_population(), and saves an
# .RData you can feed the designer as founders data.
#
# Usage:  Rscript data-raw/soybase_soysnp50k.R [n_accessions] [out.RData]
#         (defaults: 300 accessions -> SoySNP50K_div.RData)

suppressMessages(library(simplePHENOTYPES))
args  <- commandArgs(trailingOnly = TRUE)
n_acc <- if (length(args) >= 1) as.integer(args[[1]]) else 300L
out   <- if (length(args) >= 2) args[[2]] else "SoySNP50K_div.RData"

base <- "https://data.soybase.org/Glycine/max/diversity/Wm82.gnm2.div.Song_Hyten_2015"
url  <- file.path(base, "glyma.Wm82.gnm2.div.Song_Hyten_2015.hmp.gz")

cache <- "soybase_cache"
dir.create(cache, showWarnings = FALSE)
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

# sanity check + save
pop <- as_population(SoySNP50K_div,
                     individuals = seq_len(min(10L, ncol(SoySNP50K_div) - 5L)))
message(sprintf("Built SoySNP50K_div: %d SNPs x %d accessions (as_population OK: %s).",
                nrow(SoySNP50K_div), ncol(SoySNP50K_div) - 5L,
                inherits(pop, "Population")))
save(SoySNP50K_div, file = out)
message("Saved -> ", normalizePath(out))

# ---- This is how the bundled dataset was built -----------------------------
# breedingDesigner already ships a curated subset as data(SoySNP50K_div) (12,000
# MAF>=0.05 SNPs x 300 accessions; see R/data.R), so the canvas offers it directly
# in the Genotype menu and the soybean templates use it. Run this script to rebuild
# or enlarge that subset (edit n_accessions / the MAF & thinning in the finalize
# step below), or to produce your own .RData:
#   run_design(design, data_objects = list(SoySNP50K_div = SoySNP50K_div))
#   # ...where a `founders` node has params$data == "SoySNP50K_div".
#
# To regenerate the bundled data/SoySNP50K_div.rda from this object:
#   d <- readRDS("soybase_cache/SoySNP50K_div_full.rds")
#   G <- as.matrix(d[,-(1:5)]); p <- (rowMeans(G)+1)/2; maf <- pmin(p,1-p)
#   keep <- which(maf >= 0.05); if (length(keep) > 12000)
#     keep <- keep[round(seq(1, length(keep), length.out = 12000))]
#   SoySNP50K_div <- d[keep, , drop = FALSE]; rownames(SoySNP50K_div) <- NULL
#   save(SoySNP50K_div, file = "data/SoySNP50K_div.rda", compress = "xz")
