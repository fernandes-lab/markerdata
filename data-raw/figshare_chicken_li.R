#!/usr/bin/env Rscript
# data-raw/figshare_chicken_li.R
#
# Build the Chicken_div marker dataset for the markerdata package, in the SAME numeric
# format as simplePHENOTYPES::SNP55K_maize282_maf04 (snp, allele, chr, pos, cm,
# then one -1/0/1 column per individual); packed and released as for the others.
#
# Source: Figshare article 16627099 (GWAS_Chicken.bed/.bim/.fam), 839 F2 hens
# (White Leghorn x Dongxiang; per the related BMC Genomics 2021 paper), Affymetrix
# Axiom 600K chicken array, 435,867 SNPs. License: CC0.
#   Li Q (2021) GWAS_Chicken data. figshare. doi:10.6084/m9.figshare.16627099.v1
# Related paper: Li Q, Duan Z, Sun C, Zheng J, Xu G, Yang N (2021) BMC Genomics 22:786,
#   doi:10.1186/s12864-021-08103-1 (association with the deposit UNVERIFIED).
# Assembly not stated in the record; positions fit galGal4 (UNVERIFIED).
# NOTE: the deposited .fam has 16 columns (intensity/QC values in cols 3-16); only
# FID/IID are used. Individual ids are the CEL file names.
#
# PLINK binary input (read with plink2 --export A; allele = other/counted, the
# counted allele's copies minus 1 gives -1/0/1). Keeps autosomes only: autosomes labelled 1-28 (labels 29 and 30, 188 SNPs, dropped: not verifiable as autosomes; no Z/W labels exist in the .bim).
# Missing calls are imputed to the rounded per-SNP mean (as Cattle_div); SNPs are
# MAF-filtered (>= 0.05) and thinned evenly to at most 12,000 SNPs.
# Heterozygosity is present (outbred), so a dominance layer can be simulated.
#
# Usage:  Rscript data-raw/figshare_chicken_li.R [bed_prefix] [n_snps] [out.rda]
#   defaults: data-raw/chicken/GWAS_Chicken, 12000 SNPs -> data-raw/built/Chicken_div.rda.

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
prefix <- if (length(args) >= 1) args[[1]] else "data-raw/chicken/GWAS_Chicken"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/Chicken_div.rda"
stopifnot(file.exists(paste0(prefix, c(".bed", ".bim", ".fam"))))
autosomes <- 1:28

tmp <- file.path(tempdir(), "Chicken_div"); dir.create(tmp, showWarnings = FALSE)

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

Chicken_div <- data.frame(G[ord, 1:5, drop = FALSE],
                 as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                 check.names = FALSE, stringsAsFactors = FALSE)
rownames(Chicken_div) <- NULL

## 5. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(Chicken_div[, 6:ncol(Chicken_div)]) %in% c(-1L, 0L, 1L)))
pkgload::load_all(".", quiet = TRUE)
.check_marker(Chicken_div, "Chicken_div")
gt <- table(unlist(Chicken_div[, 6:ncol(Chicken_div)]))
message(sprintf("Built Chicken_div: %d SNPs x %d individuals.", nrow(Chicken_div), ncol(Chicken_div) - 5L))
message("Chromosomes kept: ", paste(unique(Chicken_div$chr), collapse = " "))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Chicken_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
