#!/usr/bin/env Rscript
# data-raw/dryad_salmon_kijas.R
#
# Build the Salmon_div marker dataset for the markerdata package, in the SAME numeric
# format as simplePHENOTYPES::SNP55K_maize282_maf04 (snp, allele, chr, pos, cm,
# then one -1/0/1 column per individual); packed and released as for the others.
#
# Source: farmed Tasmanian Atlantic salmon (Salmo salar), 777 fish, Affymetrix
# Axiom array; 106,492 SNPs in the deposited PLINK files (the record text says
# 218,132 genotyped before QC). Zenodo record 5013066 (Salmon_SNP.zip), mirror of
# Dryad doi:10.5061/dryad.pm354. License: CC0.
#   Kijas J, Elliot N, Kube P, Evans B, Botwright N, King H, Primmer CR, Verbyla K
#   (2017) Diversity and linkage disequilibrium in farmed Tasmanian Atlantic salmon.
#   Animal Genetics 48(2):237-241, doi:10.1111/age.12513
# Assembly not stated in the deposit (UNVERIFIED).
#
# PLINK binary input (read with plink2 --export A; allele = other/counted, the
# counted allele's copies minus 1 gives -1/0/1). Keeps autosomes only: ssa01-ssa29 (labels 1-29), all kept.
# Missing calls are imputed to the rounded per-SNP mean (as Cattle_div); SNPs are
# MAF-filtered (>= 0.05) and thinned evenly to at most 12,000 SNPs.
# Heterozygosity is present (outbred), so a dominance layer can be simulated.
#
# Usage:  Rscript data-raw/dryad_salmon_kijas.R [bed_prefix] [n_snps] [out.rda]
#   defaults: data-raw/salmon/x/Salmon_AllFinal, 12000 SNPs -> data-raw/built/Salmon_div.rda.

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
prefix <- if (length(args) >= 1) args[[1]] else "data-raw/salmon/x/Salmon_AllFinal"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/Salmon_div.rda"
stopifnot(file.exists(paste0(prefix, c(".bed", ".bim", ".fam"))))
autosomes <- 1:29

tmp <- file.path(tempdir(), "Salmon_div"); dir.create(tmp, showWarnings = FALSE)

## 1. sanitised PLINK copy (id columns of the .fam only) + plink2 export -----
bim <- data.table::fread(paste0(prefix, ".bim"), header = FALSE, data.table = FALSE,
                         colClasses = list(character = c(1, 2, 5, 6)))
stopifnot(sum(bim$V1 %in% autosomes & bim$V4 > 0) > 1000, !anyDuplicated(bim$V2))
fam <- data.table::fread(paste0(prefix, ".fam"), header = FALSE, data.table = FALSE,
                         colClasses = "character")
stopifnot(!anyDuplicated(fam$V2), all(nzchar(fam$V2)))
data.table::fwrite(data.frame(fam$V1, fam$V2, 0, 0, 0, -9), file.path(tmp, "p.fam"),
                   sep = " ", col.names = FALSE)
file.symlink(normalizePath(paste0(prefix, ".bed")), file.path(tmp, "p.bed"))
file.symlink(normalizePath(paste0(prefix, ".bim")), file.path(tmp, "p.bim"))
rc <- system2("plink2", c("--bfile", file.path(tmp, "p"), "--chr-set", "40", "no-xy", "--chr", paste(autosomes, collapse = ","),
                          "--maf", "0.05", "--export", "A", "--out", file.path(tmp, "e"),
                          "--memory", "6000"), stdout = FALSE)
stopifnot(rc == 0)
raw <- data.table::fread(file.path(tmp, "e.raw"), data.table = FALSE, sep = "\t")
iid <- as.character(raw$IID)
M <- t(as.matrix(raw[, -(1:6), drop = FALSE])); rm(raw); invisible(gc())
counted <- sub(".*_", "", rownames(M)); rownames(M) <- sub("_[^_]*$", "", rownames(M))
colnames(M) <- iid

## 2. map ------------------------------------------------------------------
ix <- match(rownames(M), bim$V2); stopifnot(!anyNA(ix))
b  <- bim[ix, ]
other <- ifelse(counted == b$V5, b$V6, b$V5)
G <- data.frame(snp = b$V2, allele = paste0(other, "/", counted), chr = as.integer(b$V1),
                pos = as.integer(b$V4), cm = NA_real_, check.names = FALSE,
                stringsAsFactors = FALSE)
message(sprintf("Autosomal SNPs with raw MAF >= 0.05: %d x %d individuals.", nrow(G), ncol(M)))

## 3. impute residual missing (rounded per-SNP mean), recode to -1/0/1 ------
M <- M - 1
na_frac <- mean(is.na(M))
if (na_frac > 0) {
  message(sprintf("Imputing %.3f%% missing calls to the rounded per-SNP mean.", 100 * na_frac))
  rmn <- rowMeans(M, na.rm = TRUE); rmn[is.na(rmn)] <- 0
  w <- which(is.na(M), arr.ind = TRUE)
  M[w] <- pmax(-1, pmin(1, round(rmn[w[, 1]])))
}
storage.mode(M) <- "integer"

## 4. synthetic map + MAF filter + genome-wide thin ------------------------
G$chr <- sprintf("%02d", G$chr)
G$cm  <- synthetic_map(G$chr, G$pos)
p   <- (rowMeans(M) + 1) / 2
maf <- pmin(p, 1 - p)
keep <- which(maf >= 0.05)
message(sprintf("SNPs with MAF >= 0.05: %d / %d.", length(keep), length(maf)))
ord <- keep[order(G$chr[keep], G$pos[keep])]
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]

Salmon_div <- data.frame(G[ord, 1:5, drop = FALSE],
                 as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                 check.names = FALSE, stringsAsFactors = FALSE)
rownames(Salmon_div) <- NULL

## 5. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(Salmon_div[, 6:ncol(Salmon_div)]) %in% c(-1L, 0L, 1L)))
pkgload::load_all(".", quiet = TRUE)
.check_marker(Salmon_div, "Salmon_div")
gt <- table(unlist(Salmon_div[, 6:ncol(Salmon_div)]))
message(sprintf("Built Salmon_div: %d SNPs x %d individuals.", nrow(Salmon_div), ncol(Salmon_div) - 5L))
message("Chromosomes kept: ", paste(unique(Salmon_div$chr), collapse = " "))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Salmon_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
