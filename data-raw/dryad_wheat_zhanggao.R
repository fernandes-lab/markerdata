#!/usr/bin/env Rscript
# data-raw/dryad_wheat_zhanggao.R
#
# Build a *hexaploid bread-wheat* marker dataset (WheatHex_div) in the SAME
# numeric format as simplePHENOTYPES::SNP55K_maize282_maf04 (snp, allele, chr,
# pos, cm, then one -1/0/1 column per line).
#
# STATUS: this is a SECOND wheat option, kept alongside the tetraploid TGC
# panel built by data-raw/t3_wheat_tgc.R (Wheat_div). WheatHex_div is not yet
# released or cataloged (it has no entry in inst/extdata/catalog.csv / no
# marker_data("WheatHex_div") support) -- add it to data-raw/catalog.R and run
# data-raw/release.R once it is ready to publish.
#
# Source: an Axiom-genotyped bread-wheat panel (159 lines, 53,063 SNPs with
# chromosome + position + Allele_A/Allele_B), deposited (CC0) at
#   Dryad doi:10.5061/dryad.4xgxd25ht  (original_genotype_data.xlsb)
#     Zhang & Gao (2024), Dryad -- "Original genotype data of 159 wheat samples".
# License: CC0 (public domain, Dryad); cite the Dryad deposit above.
# NOTE: this is a data-only deposit (no linked journal article); it was chosen as a
# no-account, single-file, MAPPED wheat set. A stronger breeding-program reference
# is the Gulf Atlantic Wheat Nursery panel (Sarinelli et al. 2019, TAG 132:1247;
# data on T3/Wheat, account required) -- swap to it if a cleaner citation is wanted.
#
# The .xlsb (Excel binary) sheet has columns probeset_id, Chromosome (1A..7D),
# position, Allele_A, Allele_B, then one two-nucleotide genotype per line
# ("TT"/"TG"/"GG", "NA" = missing). An embedded Python (pyxlsb) dumps it to TSV;
# R then recodes each SNP to -1/0/1 via its Allele_A/Allele_B, keeps the 21 nuclear
# chromosomes, adds a synthetic map, MAF-filters (>=0.05) and thins evenly to a
# curated 12,000-SNP subset -> data-raw/built/WheatHex_div.rda. Hexaploid wheat
# lines are largely homozygous (additive schemes; cf. SoySNP50K_div).
#
# Requires: python3 with `pyxlsb` (pip install --user pyxlsb).
#
# Usage:  Rscript data-raw/dryad_wheat_zhanggao.R [genotype.xlsb] [n_snps] [out.rda]
#   defaults: data-raw/wheat/original_genotype_data.xlsb, 12000 SNPs ->
#             data-raw/built/WheatHex_div.rda.

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else
  "data-raw/wheat/original_genotype_data.xlsb"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/WheatHex_div.rda"
stopifnot(file.exists(infile))

tmp <- tempdir()

## 1. dump the .xlsb sheet to TSV via pyxlsb --------------------------------
py  <- file.path(tmp, "xlsb_dump.py")
tsv <- file.path(tmp, "wheat.tsv")
writeLines(c(
  "import sys, csv, pyxlsb",
  "inx, out = sys.argv[1], sys.argv[2]",
  "wb = pyxlsb.open_workbook(inx)",
  "with wb.get_sheet(wb.sheets[0]) as ws, open(out, 'w', newline='') as fo:",
  "    w = csv.writer(fo, delimiter='\\t')",
  "    for row in ws.rows():",
  "        w.writerow(['' if c.v is None else c.v for c in row])"
), py)
st <- system2("python3", c(shQuote(py), shQuote(normalizePath(infile)), shQuote(tsv)))
stopifnot(st == 0L, file.exists(tsv))
dat <- data.table::fread(tsv, header = TRUE, sep = "\t", data.table = FALSE,
                         colClasses = "character")
message(sprintf("xlsb: %d SNPs x %d lines.", nrow(dat), ncol(dat) - 5L))

## 2. recode two-nucleotide genotypes -> -1/0/1 via Allele_A/Allele_B ------
# Columns by position (the deposit's 5th header is misspelled "Allel_B"):
# 1 probeset_id, 2 Chromosome, 3 position, 4 Allele_A, 5 Allele_B, 6.. genotypes.
A <- dat[[4]]; B <- dat[[5]]
Gn <- as.matrix(dat[, 6:ncol(dat), drop = FALSE])
aa <- paste0(A, A); bb <- paste0(B, B); h1 <- paste0(A, B); h2 <- paste0(B, A)
M <- matrix(NA_integer_, nrow(Gn), ncol(Gn), dimnames = list(NULL, colnames(Gn)))
M[Gn == aa] <- -1L                                # a1/a2/het recycle down columns
M[Gn == bb] <-  1L
M[Gn == h1 | Gn == h2] <- 0L

## 3. parse chr/pos; keep the 21 nuclear chromosomes (1A..7D) ---------------
chr <- dat[[2]]
pos <- suppressWarnings(as.integer(as.numeric(dat[[3]])))
biallelic <- !is.na(A) & !is.na(B) & nchar(A) == 1L & nchar(B) == 1L & A != B
ok  <- grepl("^[1-7][ABD]$", chr) & !is.na(pos) & biallelic
G <- data.frame(snp = dat[[1]][ok], allele = paste(A[ok], B[ok], sep = "/"),
                chr = chr[ok], pos = pos[ok], cm = NA_real_,
                check.names = FALSE, stringsAsFactors = FALSE)
M <- M[ok, , drop = FALSE]
o <- order(G$chr, G$pos); G <- G[o, , drop = FALSE]; M <- M[o, , drop = FALSE]
message(sprintf("Kept %d biallelic SNPs on chromosomes %s.",
                nrow(G), paste(sort(unique(G$chr)), collapse = ",")))

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

WheatHex_div <- data.frame(G[ord, 1:5, drop = FALSE],
                           as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                           check.names = FALSE, stringsAsFactors = FALSE)
rownames(WheatHex_div) <- NULL

## 6. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(WheatHex_div[, 6:ncol(WheatHex_div)]) %in% c(-1L, 0L, 1L)))
pop <- as_population(WheatHex_div, individuals = seq_len(min(10L, ncol(WheatHex_div) - 5L)))
gt  <- table(unlist(WheatHex_div[, 6:ncol(WheatHex_div)]))
message(sprintf("Built WheatHex_div: %d SNPs x %d lines (as_population OK: %s).",
                nrow(WheatHex_div), ncol(WheatHex_div) - 5L, inherits(pop, "Population")))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(WheatHex_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
