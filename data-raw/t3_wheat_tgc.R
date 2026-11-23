#!/usr/bin/env Rscript
# data-raw/t3_wheat_tgc.R
#
# Build the *wheat* marker dataset (Wheat_div) for the markerdata package, in
# the SAME numeric format as simplePHENOTYPES::SNP55K_maize282_maf04 (snp, allele,
# chr, pos, cm, then one -1/0/1 column per accession). markerdata does not ship
# this object; it publishes the packed result as a GitHub release asset that
# marker_data() downloads and caches on first use (see data-raw/release.R).
#
# Source: the Tetraploid wheat Global Collection (TGC) -- wild and domesticated
# emmer, durum landraces and tetraploid subspecies -- genotyped on the wheat 90K
# iSelect SNP array, distributed via T3/Wheat (Triticeae Toolbox) as
#   2020-11-13_..._TGC-GTC_collection_SNP90K_mapped_filtered.vcf
# from
#   Maccaferri M., Harris N.S., Twardziok S.O., ... Cattivelli L., Tuberosa R.
#     (2019) "Durum wheat genome highlights past domestication signatures and
#     future improvement targets", Nature Genetics 51:885-895,
#     doi:10.1038/s41588-019-0381-3   (reference per USDA GrainGenes,
#     "Global Durum Genomic Resources").
# License: T3/Wheat public data (Toronto Statement: cite the data source and T3;
# USDA catalog label: public domain). Redistribution terms
# have not been confirmed with T3/Wheat -- see data-raw/README.md.
#
# The VCF (## meta, then a #CHROM header + one row per SNP) is recoded to -1/0/1
# (alt-allele dosage minus 1) from the GT field, restricted to biallelic SNPs on
# the nuclear chromosomes, given a synthetic map, MAF-filtered (>=0.05) and thinned
# evenly across chromosomes to a curated 12,000-SNP subset ->
# data-raw/built/Wheat_div.rda. Tetraploid wheat lines are largely homozygous
# (additive schemes; cf. SoySNP50K_div).
#
# POLYPLOIDY NOTE: durum/emmer wheat is an allotetraploid (2n=4x=28, AABB) with
# DISOMIC inheritance (homologues pair within a subgenome; homoeologues do not,
# under the Ph1 system). A subgenome-specific SNP on chromosome 1A therefore
# segregates as an ordinary diploid locus linked only to other 1A markers, which
# is exactly what this script's diploid-per-subgenome coding (VCF GT 0/0, 0/1,
# 1/1 -> -1/0/1, one genetic map per chromosome 1A..7B) represents: the model is
# not an approximation of the inheritance. The only approximation is at the
# assay level (the 90K array's probes are subgenome-specific).
#
# Usage:  Rscript data-raw/t3_wheat_tgc.R [collection.vcf(.gz)] [n_snps] [out.rda]
#   defaults: the first *.vcf(.gz) under data-raw/wheat/, 12000 SNPs ->
#             data-raw/built/Wheat_div.rda.

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else
  list.files("data-raw/wheat", "\\.vcf(\\.gz)?$", full.names = TRUE)[1]
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/Wheat_div.rda"
stopifnot(!is.na(infile), file.exists(infile))
message("VCF: ", infile, " (", round(file.size(infile)/1e6), " MB)")

## 1. read VCF (skip ## meta; #CHROM is the header) ------------------------
reader <- if (grepl("\\.gz$", infile)) {
  sprintf("gzip -dc %s | grep -v '^##'", shQuote(normalizePath(infile)))
} else {
  sprintf("grep -v '^##' %s", shQuote(normalizePath(infile)))
}
v <- data.table::fread(cmd = reader, header = TRUE, sep = "\t",
                       data.table = FALSE, colClasses = "character")
names(v)[1] <- "CHROM"
message(sprintf("VCF: %d SNPs x %d samples.", nrow(v), ncol(v) - 9L))

## 2. recode GT -> alt dosage -1/0/1 ---------------------------------------
GT <- as.matrix(v[, 10:ncol(v), drop = FALSE])
GT <- sub(":.*$", "", GT)                          # keep the GT subfield
GT <- gsub("|", "/", GT, fixed = TRUE)             # unphase
alt <- matrix(NA_integer_, nrow(GT), ncol(GT), dimnames = list(NULL, colnames(GT)))
alt[GT == "0/0"] <- 0L
alt[GT == "0/1" | GT == "1/0"] <- 1L
alt[GT == "1/1"] <- 2L
M <- alt - 1L                                      # -1/0/1 (dosage of ALT), NA else

## 3. assemble; keep biallelic SNPs on nuclear chromosomes -----------------
chr <- toupper(sub("^CHR", "", v$CHROM, ignore.case = TRUE))
pos <- suppressWarnings(as.integer(v$POS))
biallelic <- nchar(v$REF) == 1L & nchar(v$ALT) == 1L & !grepl(",", v$ALT)
ok  <- grepl("^[1-7][ABD]$", chr) & !is.na(pos) & biallelic
G <- data.frame(snp = ifelse(v$ID == "." | v$ID == "", paste0(chr, "_", pos), v$ID)[ok],
                allele = paste(v$REF[ok], v$ALT[ok], sep = "/"),
                chr = chr[ok], pos = pos[ok], cm = NA_real_,
                check.names = FALSE, stringsAsFactors = FALSE)
M <- M[ok, , drop = FALSE]
o <- order(G$chr, G$pos); G <- G[o, , drop = FALSE]; M <- M[o, , drop = FALSE]
G$snp <- make.unique(G$snp)                        # IWB probe IDs repeat across subgenomes
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

Wheat_div <- data.frame(G[ord, 1:5, drop = FALSE],
                        as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                        check.names = FALSE, stringsAsFactors = FALSE)
rownames(Wheat_div) <- NULL

## 6. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(Wheat_div[, 6:ncol(Wheat_div)]) %in% c(-1L, 0L, 1L)))
pop <- as_population(Wheat_div, individuals = seq_len(min(10L, ncol(Wheat_div) - 5L)))
gt  <- table(unlist(Wheat_div[, 6:ncol(Wheat_div)]))
message(sprintf("Built Wheat_div: %d SNPs x %d lines (as_population OK: %s).",
                nrow(Wheat_div), ncol(Wheat_div) - 5L, inherits(pop, "Population")))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Wheat_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
