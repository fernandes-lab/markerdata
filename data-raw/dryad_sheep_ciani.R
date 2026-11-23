#!/usr/bin/env Rscript
# data-raw/dryad_sheep_ciani.R
#
# Build the Sheep_div marker dataset for the markerdata package, in the SAME numeric
# format as simplePHENOTYPES::SNP55K_maize282_maf04 (snp, allele, chr, pos, cm,
# then one -1/0/1 column per individual); packed and released as for the others.
#
# Source: Merino and Merino-derived sheep breeds + other breeds (37 populations, incl. a few
# wild Ovis), 671 animals, Illumina OvineSNP50, 45,847 SNPs. Dryad
# doi:10.5061/dryad.2p0qf (CIANI.rar; fetched from the identical Zenodo mirror
# record 4954171, md5 c8491556d590170015586cd85187e430 = Dryad digest). License: CC0
# (Dryad API: https://spdx.org/licenses/CC0-1.0.html).
#   Ciani E, Lasagna E, D'Andrea M, Alloggio I, Marroni F, Ceccobelli S, Delgado
#   Bermejo JV, Sarti FM, Kijas J, Lenstra JA, Pilla F (2015) Merino and
#   Merino-derived sheep breeds: a genome-wide intercontinental study. Genetics
#   Selection Evolution 47:64, doi:10.1186/s12711-015-0139-z
# Assembly not stated in the deposit (UNVERIFIED). Label 0 (282 unplaced SNPs) and 27
# (assumed X, UNVERIFIED) are dropped.
#
# PLINK binary input (read with plink2 --export A; allele = other/counted, the
# counted allele's copies minus 1 gives -1/0/1). Keeps autosomes only: autosomes 1-26.
# Missing calls are imputed to the rounded per-SNP mean (as Cattle_div); SNPs are
# MAF-filtered (>= 0.05) and thinned evenly to at most 12,000 SNPs.
# Heterozygosity is present (outbred), so a dominance layer can be simulated.
#
# Usage:  Rscript data-raw/dryad_sheep_ciani.R [bed_prefix] [n_snps] [out.rda]
#   defaults: data-raw/sheep/x/CIANI_GSE_g14091_MS196762, 12000 SNPs -> data-raw/built/Sheep_div.rda.

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
prefix <- if (length(args) >= 1) args[[1]] else "data-raw/sheep/x/CIANI_GSE_g14091_MS196762"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/Sheep_div.rda"
stopifnot(file.exists(paste0(prefix, c(".bed", ".bim", ".fam"))))
autosomes <- 1:26

tmp <- file.path(tempdir(), "Sheep_div"); dir.create(tmp, showWarnings = FALSE)

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

Sheep_div <- data.frame(G[ord, 1:5, drop = FALSE],
                 as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                 check.names = FALSE, stringsAsFactors = FALSE)
rownames(Sheep_div) <- NULL

## 5. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(Sheep_div[, 6:ncol(Sheep_div)]) %in% c(-1L, 0L, 1L)))
pkgload::load_all(".", quiet = TRUE)
.check_marker(Sheep_div, "Sheep_div")
gt <- table(unlist(Sheep_div[, 6:ncol(Sheep_div)]))
message(sprintf("Built Sheep_div: %d SNPs x %d individuals.", nrow(Sheep_div), ncol(Sheep_div) - 5L))
message("Chromosomes kept: ", paste(unique(Sheep_div$chr), collapse = " "))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Sheep_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
