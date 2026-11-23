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
# The full HapMap is large (20,087 accessions x 42,509 SNPs), so this script
# (1) downloads it once and caches it under data-raw/soybase_cache/, (2)
# streams it in chunks of SNP rows (never holding the whole text), converting
# each chunk to the numeric -1/0/1 format with
# simplePHENOTYPES::as_numeric(from = "hapmap"), (3) drops SNPs with more than
# MAX_MISS missing calls or MAF < 0.05 (computed on the observed calls), and
# imputes the remaining missing calls to the major homozygote (the same rule
# as data-raw/s1database_panels.R), (4) adds a synthetic genetic map (the
# SoyBase files carry no cM), (5) thins evenly across chromosomes to a curated
# n_snps-SNP subset, and (6) sanity-checks with as_population() and saves
# data-raw/built/SoySNP50K_div.rda. ALL accessions are kept by default.
#
# Usage:  Rscript data-raw/soybase_soysnp50k.R [n_accessions|all] [n_snps] [out.rda]
#   defaults: all accessions, 12000 SNPs -> data-raw/built/SoySNP50K_div.rda
#   (a number keeps only the first n_accessions accessions)

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
n_acc  <- if (length(args) >= 1 && args[[1]] != "all") as.integer(args[[1]]) else Inf
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/SoySNP50K_div.rda"

MAX_MISS <- 0.05  # drop SNPs missing in more than this fraction of accessions
CHUNK    <- 2500L # SNP rows converted per step

base <- "https://data.soybase.org/Glycine/max/diversity/Wm82.gnm2.div.Song_Hyten_2015"
url  <- file.path(base, "glyma.Wm82.gnm2.div.Song_Hyten_2015.hmp.gz")

cache <- "data-raw/soybase_cache"
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
gz <- file.path(cache, basename(url))
if (!file.exists(gz)) {
  message("Downloading SoySNP50K HapMap (~94 MB) ...")
  utils::download.file(url, gz, mode = "wb")
}

# Stream the HapMap in chunks of CHUNK SNP rows. Each chunk is written as a small
# HapMap file, converted by as_numeric() (-1/0/1, NA = missing call), filtered
# and stored as a raw matrix (genotype + 1) so the whole panel stays in memory
# at 1 byte per call. The map columns of every SNP are kept to build the
# synthetic cM map on the full SNP set before thinning.
con <- gzfile(gz, "rt")
hdr <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
keep_col <- seq_len(min(length(hdr), 11L + n_acc)) # 11 HapMap metadata columns + accessions
limited <- length(keep_col) < length(hdr)
tmp <- tempfile(fileext = ".hmp.txt")
map <- list(); geno <- list(); n_read <- 0L
repeat {
  ln <- readLines(con, n = CHUNK)
  if (!length(ln)) break
  if (limited) {
    ln <- vapply(strsplit(ln, "\t", fixed = TRUE),
                 function(x) paste(x[keep_col], collapse = "\t"), character(1))
  }
  h <- paste(hdr[keep_col], collapse = "\t")
  writeLines(c(h, ln), tmp)
  x <- as_numeric(tmp, from = "hapmap", to_r = TRUE, to_file = FALSE, verbose = FALSE)
  G <- as.matrix(x[, 6:ncol(x), drop = FALSE])
  miss <- rowMeans(is.na(G))
  p <- (rowMeans(G, na.rm = TRUE) + 1) / 2
  ok <- miss <= MAX_MISS & pmin(p, 1 - p) >= 0.05
  ok[is.na(ok)] <- FALSE
  map[[length(map) + 1L]] <- data.frame(snp = x$snp, allele = x$allele, chr = x$chr,
                                        pos = x$pos, cm = NA_real_, ok = ok,
                                        stringsAsFactors = FALSE)
  G <- G[ok, , drop = FALSE]
  if (anyNA(G)) { # impute missing calls to the major homozygote
    maj <- ifelse(rowSums(G == 1L, na.rm = TRUE) > rowSums(G == -1L, na.rm = TRUE), 1L, -1L)
    G[is.na(G)] <- maj[row(G)[is.na(G)]]
  }
  geno[[length(geno) + 1L]] <- matrix(as.raw(G + 1L), nrow(G), ncol(G),
                                      dimnames = list(NULL, colnames(G)))
  n_read <- n_read + length(ln)
  message(sprintf("  %d SNPs read, %d kept so far", n_read, sum(vapply(geno, nrow, 1L))))
}
close(con)
unlink(tmp)
map  <- do.call(rbind, map)
geno <- do.call(rbind, geno)
message(sprintf("SNPs kept (missing <= %.2f, MAF >= 0.05): %d / %d (chromosome labels: %d distinct).",
                MAX_MISS, nrow(geno), nrow(map), length(unique(map$chr[map$ok]))))

# The SoyBase HapMap has no genetic map; build one from the physical positions of
# ALL SNPs so meiosis can be simulated (the engine errors clearly otherwise).
map$cm <- synthetic_map(map$chr, map$pos)

## finalize: thin evenly across chromosomes to n_snps -------------------------
# This is the step that turns the full converted panel into the curated,
# released subset (mirrors the finalize step in the other data-raw/*.R
# scripts). `idx` maps rows of `geno` back to rows of `map`.
idx <- which(map$ok)
ord <- order(map$chr[idx], map$pos[idx])
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]
SoySNP50K_div <- cbind(map[idx[ord], c("snp", "allele", "chr", "pos", "cm")],
                       as.data.frame(matrix(as.integer(geno[ord, , drop = FALSE]) - 1L,
                                            length(ord), ncol(geno),
                                            dimnames = list(NULL, colnames(geno))),
                                     check.names = FALSE))
rownames(SoySNP50K_div) <- NULL
rm(geno)

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
