#!/usr/bin/env Rscript
# data-raw/mendeley_tomato_51k.R
#
# Build the *tomato* marker dataset (Tomato_div) for the markerdata package,
# in the SAME numeric format as simplePHENOTYPES::SNP55K_maize282_maf04 (snp,
# allele, chr, pos, cm, then one -1/0/1 column per accession). markerdata does
# not ship this object; it publishes the packed result as a GitHub release
# asset that marker_data() downloads and caches on first use (see
# data-raw/release.R).
#
# Source: "The 51K SNP data of two tomato collections" (Axiom 51K tomato array;
# TGC1 = 162 accessions, TGC2 = 191 accessions, 353 in all), Mendeley Data
#   https://data.mendeley.com/datasets/bxcpc274fh/1   doi:10.17632/bxcpc274fh.1
#   files: "TGC1 51K SNP chip genotypic data.xlsx" (28,752,970 bytes)
#          "TGC2 51K SNP chip genotypic data.xlsx" (34,835,699 bytes)
# Deposited by Sung-Chur Sim. License: CC BY 4.0 (Creative Commons Attribution
# 4.0 International, as shown on the Mendeley Data record).
# Paper using these data:
#   Yeon J., Nguyen T.T.P., Kim M., Sim S.-C. (2024) "Prediction accuracy of
#   genomic estimated breeding values for fruit traits in cultivated tomato
#   (Solanum lycopersicum L.)", BMC Plant Biology 24:222,
#   doi:10.1186/s12870-024-04934-8
# The genome assembly is not stated in the deposit or the paper (UNVERIFIED).
#
# Each sheet has columns "SNP ID", "Chr" (1-12), "Position (bp)", then one
# two-nucleotide call per accession ("CC", "CT", ..., "NN" = missing). The two
# files list the same 51,214 SNPs in the same order. The deposit carries no
# allele columns, so for each SNP the allele observed more often across the 353
# accessions is "A" (coded -1 when homozygous), the other is "B" (coded 1), and
# a heterozygote is 0; the `allele` column records them as "major/minor". SNPs
# with more than two observed alleles are dropped. Residual missing calls
# (about 3%) are imputed to the rounded per-SNP mean, then a synthetic genetic
# map is added, SNPs are MAF-filtered (>=0.05) and thinned evenly to a curated
# 12,000-SNP subset -> data-raw/built/Tomato_div.rda. Cultivated tomato is a
# selfing species, so accessions are mostly homozygous (additive schemes; cf.
# SoySNP50K_div).
#
# Usage:  Rscript data-raw/mendeley_tomato_51k.R [tomato_dir] [n_snps] [out.rda]
#   defaults: data-raw/tomato, 12000 SNPs -> data-raw/built/Tomato_div.rda.

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
dir_in <- if (length(args) >= 1) args[[1]] else "data-raw/tomato"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/Tomato_div.rda"
files  <- file.path(dir_in, c("TGC1 51K SNP chip genotypic data.xlsx",
                              "TGC2 51K SNP chip genotypic data.xlsx"))
stopifnot(all(file.exists(files)))

## 1. read both sheets as text ---------------------------------------------
x <- lapply(files, function(f)
  as.data.frame(readxl::read_excel(f, col_types = "text"), check.names = FALSE))
stopifnot(identical(x[[1]][, 1:3], x[[2]][, 1:3]))     # same SNPs, same order
Gn  <- as.matrix(do.call(cbind, lapply(x, function(d) d[, -(1:3), drop = FALSE])))
snp <- x[[1]][[1]]; chr <- x[[1]][[2]]
pos <- suppressWarnings(as.integer(x[[1]][[3]]))
rm(x); invisible(gc())
message(sprintf("xlsx: %d SNPs x %d accessions.", nrow(Gn), ncol(Gn)))
stopifnot(!anyDuplicated(colnames(Gn)), all(nzchar(colnames(Gn))),
          !anyDuplicated(snp))

## 2. recode two-nucleotide calls -> -1/0/1 (A = more frequent allele) -----
b1 <- substr(Gn, 1L, 1L); b2 <- substr(Gn, 2L, 2L)
miss <- !(b1 %in% c("A", "C", "G", "T")) | !(b2 %in% c("A", "C", "G", "T"))
dim(b1) <- dim(b2) <- dim(miss) <- dim(Gn)
b1[miss] <- NA; b2[miss] <- NA
cnt <- sapply(c("A", "C", "G", "T"),
              function(b) rowSums(b1 == b, na.rm = TRUE) + rowSums(b2 == b, na.rm = TRUE))
n_all <- rowSums(cnt > 0)
major <- c("A", "C", "G", "T")[max.col(cnt, ties.method = "first")]
cnt2  <- cnt; cnt2[cbind(seq_len(nrow(cnt)), max.col(cnt, ties.method = "first"))] <- -1
minor <- c("A", "C", "G", "T")[max.col(cnt2, ties.method = "first")]
M <- (b1 == minor) + (b2 == minor) - 1L                 # copies of minor allele - 1
M[miss] <- NA
dimnames(M) <- list(NULL, colnames(Gn))
rm(Gn, b1, b2, miss, cnt, cnt2); invisible(gc())

## 3. keep biallelic SNPs on chromosomes 1-12 ------------------------------
ok <- n_all == 2L & chr %in% as.character(1:12) & !is.na(pos)
message(sprintf("Dropped %d SNPs (non-biallelic or not on chr 1-12).", sum(!ok)))
G <- data.frame(snp = snp[ok], allele = paste(major[ok], minor[ok], sep = "/"),
                chr = sprintf("%02d", as.integer(chr[ok])), pos = pos[ok],
                cm = NA_real_, check.names = FALSE, stringsAsFactors = FALSE)
M <- M[ok, , drop = FALSE]
o <- order(G$chr, G$pos); G <- G[o, , drop = FALSE]; M <- M[o, , drop = FALSE]

## 4. impute residual missing (rounded per-SNP mean) -----------------------
storage.mode(M) <- "double"
na_frac <- mean(is.na(M))
if (na_frac > 0) {
  message(sprintf("Imputing %.3f%% missing calls to the rounded per-SNP mean.",
                  100 * na_frac))
  rmn <- rowMeans(M, na.rm = TRUE); rmn[is.na(rmn)] <- 0
  ix <- which(is.na(M), arr.ind = TRUE)
  M[ix] <- pmax(-1, pmin(1, round(rmn[ix[, 1]])))
}
storage.mode(M) <- "integer"

## 5. synthetic map + MAF filter + genome-wide thin ------------------------
G$cm <- synthetic_map(G$chr, G$pos)
p   <- (rowMeans(M) + 1) / 2
maf <- pmin(p, 1 - p)
keep <- which(maf >= 0.05)
message(sprintf("SNPs with MAF >= 0.05: %d / %d.", length(keep), length(maf)))
ord <- keep[order(G$chr[keep], G$pos[keep])]
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]

Tomato_div <- data.frame(G[ord, 1:5, drop = FALSE],
                         as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                         check.names = FALSE, stringsAsFactors = FALSE)
rownames(Tomato_div) <- NULL

## 6. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(Tomato_div[, 6:ncol(Tomato_div)]) %in% c(-1L, 0L, 1L)))
pop <- as_population(Tomato_div, individuals = seq_len(min(10L, ncol(Tomato_div) - 5L)))
gt  <- table(unlist(Tomato_div[, 6:ncol(Tomato_div)]))
message(sprintf("Built Tomato_div: %d SNPs x %d accessions (as_population OK: %s).",
                nrow(Tomato_div), ncol(Tomato_div) - 5L, inherits(pop, "Population")))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Tomato_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
