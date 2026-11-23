#!/usr/bin/env Rscript
# data-raw/dryad_potato_endelman.R
#
# Build the *potato* marker dataset (Potato_div) for the markerdata package, in
# the SAME numeric format as simplePHENOTYPES::SNP55K_maize282_maf04 (snp, allele,
# chr, pos, cm, then one -1/0/1 column per clone). Published as a GitHub release
# asset (see data-raw/release.R).
#
# Source: 586 autotetraploid potato (Solanum tuberosum, 2n=4x=48) clones, SNP
# array, Zenodo record 6818259 = Dryad doi:10.5061/dryad.3n5tb2rk0, file
# marker_data.csv (README: "Allele dosages (encoded 0-4) ... 10,185 SNPs ...
# positions from the DMv4.03 reference genome"). License: CC0 (Zenodo metadata).
#   Caraza-Harter M.V., Endelman J.B. (2022) "The genetic architectures of vine
#     and skin maturity in tetraploid potato", Theor. Appl. Genet.,
#     doi:10.1007/s00122-022-04159-z   (data: Endelman & Caraza-Harter, Zenodo
#     record 6818259 / Dryad doi:10.5061/dryad.3n5tb2rk0)
# Population: GWAS panel of breeding clones (Wisconsin potato program).
#
# !! ALLELE DOSAGE IS LOST !! The package format holds only -1/0/1, so the
# dosage 0..4 is DIPLOIDIZED: dosage 0 (nulliplex) -> -1, dosage 4 (quadruplex)
# -> 1, dosages 1,2,3 -> 0. The "allele" column is "N/N" (the file does not give
# nucleotides). MAF is computed from the ORIGINAL dosage (mean dosage / 4).
# Chromosomes are 1..12 (zero-padded "01".."12"); no unplaced contigs in the file.
# No missing calls in the source; any NA would be imputed like the other scripts.
#
# Usage:  Rscript data-raw/dryad_potato_endelman.R [marker_data.csv] [n_snps] [out.rda]

suppressMessages({library(simplePHENOTYPES); pkgload::load_all(".", quiet = TRUE)})
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else "data-raw/potato/marker_data.csv"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/Potato_div.rda"
ploidy <- 4L
stopifnot(file.exists(infile))

## 1. read ----------------------------------------------------------------
v <- data.table::fread(infile, header = TRUE, sep = ",", data.table = FALSE)
D <- as.matrix(v[, -(1:3), drop = FALSE])           # dosage 0..4
storage.mode(D) <- "double"
stopifnot(all(D %in% c(0:ploidy, NA)), !anyDuplicated(colnames(D)), all(nzchar(colnames(D))))
chr <- suppressWarnings(as.integer(v$chrom)); pos <- as.integer(v$bp)
ok  <- !is.na(chr) & chr >= 1 & chr <= 12 & !is.na(pos)
message(sprintf("%d SNPs x %d clones; %d on chromosomes 1-12.", nrow(v), ncol(D), sum(ok)))
G <- data.frame(snp = v$marker[ok], allele = "N/N", chr = sprintf("%02d", chr[ok]),
                pos = pos[ok], cm = NA_real_, stringsAsFactors = FALSE)
D <- D[ok, , drop = FALSE]
o <- order(G$chr, G$pos); G <- G[o, ]; D <- D[o, , drop = FALSE]
G$snp <- make.unique(G$snp)

## 2. MAF from the original dosage; then diploidize -------------------------
maf <- pmin(rowMeans(D, na.rm = TRUE) / ploidy, 1 - rowMeans(D, na.rm = TRUE) / ploidy)
M <- D; M[!is.na(D)] <- ifelse(D[!is.na(D)] == 0, -1, ifelse(D[!is.na(D)] == ploidy, 1, 0))
message(sprintf("Heterozygous dosage classes collapsed to 0: %.1f%% of calls.",
                100 * mean(M == 0, na.rm = TRUE)))

## 3. impute residual missing (rounded per-SNP mean) ------------------------
na_frac <- mean(is.na(M))
if (na_frac > 0) {
  message(sprintf("Imputing %.3f%% missing calls.", 100 * na_frac))
  rmn <- rowMeans(M, na.rm = TRUE); rmn[is.na(rmn)] <- 0
  ix <- which(is.na(M), arr.ind = TRUE)
  M[ix] <- pmax(-1, pmin(1, round(rmn[ix[, 1]])))
}
storage.mode(M) <- "integer"

## 4. synthetic map + MAF filter + genome-wide thin -------------------------
G$cm <- synthetic_map(G$chr, G$pos)
keep <- which(maf >= 0.05)
message(sprintf("SNPs with MAF >= 0.05: %d / %d.", length(keep), length(maf)))
ord <- keep[order(G$chr[keep], G$pos[keep])]
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]

Potato_div <- data.frame(G[ord, 1:5, drop = FALSE],
                         as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                         check.names = FALSE, stringsAsFactors = FALSE)
rownames(Potato_div) <- NULL

## 5. contract check + save -------------------------------------------------
markerdata:::.check_marker(Potato_div, "Potato_div")
gt <- table(unlist(Potato_div[, 6:ncol(Potato_div)]))
message(sprintf("Built Potato_div: %d SNPs x %d clones.", nrow(Potato_div), ncol(Potato_div) - 5L))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Potato_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
