#!/usr/bin/env Rscript
# data-raw/usda_sugarcane_ratoon.R
#
# Build the *sugarcane* marker dataset (Sugarcane_div) for the markerdata
# package, in the SAME numeric format as simplePHENOTYPES::SNP55K_maize282_maf04
# (snp, allele, chr, pos, cm, then one -1/0/1 column per clone). Published as a
# GitHub release asset (see data-raw/release.R).
#
# Source: USDA Ag Data Commons doi:10.15482/USDA.ADC/1528986, file
# sugarcane.10501.SNPs_.432.Inds_.vcf (also figshare ndownloader file 43734678).
# License: public domain (USA.gov Public Domain Label 1.0, per the data.gov record).
#   Islam M.S., Corak K., McCord P., Hulse-Kemp A.M., Lipka A.E. (2023) "A first
#     look at the ability to use genomic prediction for improving the ratooning
#     ability of sugarcane", Front. Plant Sci. 14:1205999,
#     doi:10.3389/fpls.2023.1205999
# Population: 432 clones (414 USDA-ARS Canal Point FL breeding clones + 18
# commercial cultivars). Capture-Seq reads aligned to Sorghum bicolor v3.1 (no
# sugarcane reference available), so chr 1..10 and positions are SORGHUM
# coordinates. Sugarcane is a high polyploid (assumed decaploid in the paper).
#
# !! ALLELE DOSAGE IS NOT IN THE FILE AND IS LOST !! The VCF holds PSEUDO-DIPLOID
# GT calls (0/0, 0/1, 1/1; the authors lacked depth to estimate dosage), with
# REF = major allele. They are mapped like a diploid: 0/0 -> -1, 0/1 -> 0,
# 1/1 -> 1. MAF is computed from those alt-allele counts / 2. Multi-allelic sites
# (ALT with a comma) and any non 0/1 genotype are dropped/NA. The file has no
# missing calls; any NA would be imputed like the other scripts.
#
# Usage:  Rscript data-raw/usda_sugarcane_ratoon.R [file.vcf] [n_snps] [out.rda]

suppressMessages({library(simplePHENOTYPES); pkgload::load_all(".", quiet = TRUE)})
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else
  "data-raw/sugarcane/sugarcane.10501.SNPs_.432.Inds_.vcf"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/Sugarcane_div.rda"
stopifnot(file.exists(infile))

## 1. read VCF -------------------------------------------------------------
v <- data.table::fread(cmd = sprintf("grep -v '^##' %s", shQuote(normalizePath(infile))),
                       header = TRUE, sep = "\t", data.table = FALSE, colClasses = "character")
names(v)[1] <- "CHROM"
GT <- sub(":.*$", "", as.matrix(v[, 10:ncol(v), drop = FALSE]))
GT <- gsub("|", "/", GT, fixed = TRUE)
stopifnot(!anyDuplicated(colnames(GT)), all(nzchar(colnames(GT))))
message(sprintf("VCF: %d SNPs x %d clones.", nrow(v), ncol(GT)))

## 2. alt dosage (0/1/2, pseudo-diploid) ------------------------------------
D <- matrix(NA_real_, nrow(GT), ncol(GT), dimnames = list(NULL, colnames(GT)))
D[GT == "0/0"] <- 0; D[GT == "0/1" | GT == "1/0"] <- 1; D[GT == "1/1"] <- 2

## 3. biallelic SNPs on chromosomes 1..10 (all SNPs in the file) ------------
chr <- suppressWarnings(as.integer(v$CHROM)); pos <- suppressWarnings(as.integer(v$POS))
ok <- nchar(v$REF) == 1L & nchar(v$ALT) == 1L & !is.na(chr) & chr >= 1 & chr <= 10 & !is.na(pos)
message(sprintf("%d biallelic SNPs on chromosomes 1-10 (of %d); 0 on unplaced contigs.", sum(ok), nrow(v)))
G <- data.frame(snp = v$ID[ok], allele = paste(v$REF[ok], v$ALT[ok], sep = "/"),
                chr = sprintf("%02d", chr[ok]), pos = pos[ok], cm = NA_real_,
                stringsAsFactors = FALSE)
D <- D[ok, , drop = FALSE]
o <- order(G$chr, G$pos); G <- G[o, ]; D <- D[o, , drop = FALSE]
G$snp <- make.unique(G$snp)

## 4. MAF from original calls; recode -1/0/1; impute -----------------------
p <- rowMeans(D, na.rm = TRUE) / 2; maf <- pmin(p, 1 - p)
M <- D - 1
message(sprintf("Heterozygous calls (0) : %.1f%% of calls.", 100 * mean(M == 0, na.rm = TRUE)))
na_frac <- mean(is.na(M))
if (na_frac > 0) {
  message(sprintf("Imputing %.3f%% missing calls.", 100 * na_frac))
  rmn <- rowMeans(M, na.rm = TRUE); rmn[is.na(rmn)] <- 0
  ix <- which(is.na(M), arr.ind = TRUE)
  M[ix] <- pmax(-1, pmin(1, round(rmn[ix[, 1]])))
}
storage.mode(M) <- "integer"

## 5. synthetic map + MAF filter + genome-wide thin -------------------------
G$cm <- synthetic_map(G$chr, G$pos)
keep <- which(maf >= 0.05)
message(sprintf("SNPs with MAF >= 0.05: %d / %d.", length(keep), length(maf)))
ord <- keep[order(G$chr[keep], G$pos[keep])]
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]

Sugarcane_div <- data.frame(G[ord, 1:5, drop = FALSE],
                            as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                            check.names = FALSE, stringsAsFactors = FALSE)
rownames(Sugarcane_div) <- NULL

## 6. contract check + save -------------------------------------------------
markerdata:::.check_marker(Sugarcane_div, "Sugarcane_div")
gt <- table(unlist(Sugarcane_div[, 6:ncol(Sugarcane_div)]))
message(sprintf("Built Sugarcane_div: %d SNPs x %d clones.", nrow(Sugarcane_div), ncol(Sugarcane_div) - 5L))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Sugarcane_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
