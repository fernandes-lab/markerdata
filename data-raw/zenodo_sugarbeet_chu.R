#!/usr/bin/env Rscript
# data-raw/zenodo_sugarbeet_chu.R
#
# Build the *sugar beet* marker dataset (SugarBeet_div) for the markerdata
# package, in the SAME numeric format as simplePHENOTYPES::SNP55K_maize282_maf04
# (snp, allele, chr, pos, cm, then one -1/0/1 column per accession). markerdata
# does not ship this object; it publishes the packed result as a GitHub release
# asset that marker_data() downloads and caches on first use (see
# data-raw/release.R).
#
# Source: GBS genotypes of 1,928 publicly available US Beta vulgaris germplasm
# lines (148,317 SNPs in the file; imputed, MAF >= 0.05 per the file name),
# Zenodo record https://zenodo.org/records/11551725
#   file: GBS_imputed_MA05_100bp_082421.hmp.txt (578,988,346 bytes, HapMap)
# Depositor: Chenggen Chu (USDA-ARS). License: CC BY 4.0 (Zenodo record).
# Paper: Chu C., et al. (2023) Crop Science 63:2255-2273,
#   doi:10.1002/csc2.20978 -- "Potential of publicly available Beta vulgaris
#   germplasm for sustainable sugarbeet improvement indicated by combining
#   analysis of genetic diversity and historic resistance evaluation".
#   (Full author list: UNVERIFIED -- the Zenodo record lists only Chu.)
# The genome assembly is not stated in the file or the record (UNVERIFIED).
#
# HapMap columns: rs#, alleles ("A/G"), chrom (1-9 plus "SFD_*" scaffold labels),
# pos, then 11 standard metadata columns and one nucleotide/IUPAC call per
# accession (R/Y/K/M/S/W = heterozygote, N = missing). The file is streamed in
# chunks (it is too big to load naively). For each SNP the first listed allele
# is "A" (homozygote -> -1), the second "B" (-> 1), the IUPAC code of the pair
# is the heterozygote (0); anything else is missing. Non-ACGT (indel) or
# non-biallelic SNPs are dropped. Chromosome labels are kept as in the file,
# including the SFD_* scaffold labels, as soybean (SoySNP50K_div) keeps its
# SCAFFOLD_* labels. Residual missing calls are imputed to the rounded
# per-SNP mean, then a synthetic genetic map is added, SNPs are MAF-filtered
# (>=0.05) and thinned evenly to a curated 12,000-SNP subset ->
# data-raw/built/SugarBeet_div.rda. Sugar beet is an outcrossing species
# (cultivars, breeding lines and wild sea beet), so heterozygotes are present.
#
# Usage:  Rscript data-raw/zenodo_sugarbeet_chu.R [file.hmp.txt] [n_snps] [out.rda]
#   defaults: data-raw/sugarbeet/GBS_imputed_MA05_100bp_082421.hmp.txt, 12000
#             SNPs -> data-raw/built/SugarBeet_div.rda.

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else
  "data-raw/sugarbeet/GBS_imputed_MA05_100bp_082421.hmp.txt"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/SugarBeet_div.rda"
stopifnot(file.exists(infile))

## 1. stream the HapMap in chunks, recode -> -1/0/1, keep MAF >= 0.05 ------
iupac <- c(AC = "M", CA = "M", AG = "R", GA = "R", AT = "W", TA = "W",
           CG = "S", GC = "S", CT = "Y", TC = "Y", GT = "K", TG = "K")
con <- file(infile, "rt"); on.exit(close(con))
hdr <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
ind <- hdr[-(1:11)]
stopifnot(hdr[1:4] == c("rs#", "alleles", "chrom", "pos"),
          !anyDuplicated(ind), all(nzchar(ind)))
message(sprintf("HapMap: %d accessions.", length(ind)))

meta_l <- list(); M_l <- list(); n_read <- 0L; na_n <- 0; cell_n <- 0
repeat {
  ln <- readLines(con, n = 20000L)
  if (!length(ln)) break
  d <- data.table::fread(text = ln, header = FALSE, sep = "\t", quote = "",
                         colClasses = "character", data.table = FALSE)
  n_read <- n_read + nrow(d)
  al <- d[[2]]
  a1 <- substr(al, 1L, 1L); a2 <- substr(al, 3L, 3L)
  okb <- grepl("^[ACGT]/[ACGT]$", al) & a1 != a2
  Gn <- as.matrix(d[, -(1:11), drop = FALSE])
  het <- unname(iupac[paste0(a1, a2)])
  M <- matrix(NA_integer_, nrow(Gn), ncol(Gn))
  M[Gn == a1] <- -1L                       # a1/a2/het recycle down the columns
  M[Gn == a2] <-  1L
  M[Gn == het] <-  0L
  M[!okb, ] <- NA_integer_
  na_n <- na_n + sum(is.na(M[okb, , drop = FALSE])); cell_n <- cell_n + sum(okb) * ncol(M)
  p <- (rowMeans(M, na.rm = TRUE) + 1) / 2
  keep <- okb & !is.na(p) & pmin(p, 1 - p) >= 0.05   # pre-imputation MAF screen
  meta_l[[length(meta_l) + 1L]] <- data.frame(
    snp = d[[1]][keep], allele = al[keep], chr = d[[3]][keep],
    pos = as.integer(d[[4]][keep]), cm = NA_real_, stringsAsFactors = FALSE)
  M_l[[length(M_l) + 1L]] <- M[keep, , drop = FALSE]
  rm(d, Gn, M); invisible(gc())
}
G <- do.call(rbind, meta_l); M <- do.call(rbind, M_l); rm(meta_l, M_l)
colnames(M) <- ind
message(sprintf("Read %d SNPs; %d biallelic with MAF >= 0.05 (%.2f%% missing calls).",
                n_read, nrow(G), 100 * na_n / cell_n))
stopifnot(!anyNA(G$pos), !anyDuplicated(G$snp))

## 2. impute residual missing (rounded per-SNP mean) -----------------------
na_frac <- mean(is.na(M))
if (na_frac > 0) {
  message(sprintf("Imputing %.3f%% missing calls to the rounded per-SNP mean.",
                  100 * na_frac))
  rmn <- rowMeans(M, na.rm = TRUE); rmn[is.na(rmn)] <- 0
  ix <- which(is.na(M), arr.ind = TRUE)
  M[ix] <- as.integer(pmax(-1, pmin(1, round(rmn[ix[, 1]]))))
}

## 3. synthetic map + MAF filter + genome-wide thin ------------------------
o <- order(G$chr, G$pos); G <- G[o, , drop = FALSE]; M <- M[o, , drop = FALSE]
G$cm <- synthetic_map(G$chr, G$pos)
p   <- (rowMeans(M) + 1) / 2
maf <- pmin(p, 1 - p)
keep <- which(maf >= 0.05)
message(sprintf("SNPs with MAF >= 0.05 after imputation: %d / %d (chromosome labels: %d distinct).",
                length(keep), length(maf), length(unique(G$chr[keep]))))
ord <- keep[order(G$chr[keep], G$pos[keep])]
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]

SugarBeet_div <- data.frame(G[ord, 1:5, drop = FALSE],
                            as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                            check.names = FALSE, stringsAsFactors = FALSE)
rownames(SugarBeet_div) <- NULL

## 4. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(SugarBeet_div[, 6:ncol(SugarBeet_div)]) %in% c(-1L, 0L, 1L)))
pop <- as_population(SugarBeet_div,
                     individuals = seq_len(min(10L, ncol(SugarBeet_div) - 5L)))
gt  <- table(unlist(SugarBeet_div[, 6:ncol(SugarBeet_div)]))
message(sprintf("Built SugarBeet_div: %d SNPs x %d accessions (as_population OK: %s).",
                nrow(SugarBeet_div), ncol(SugarBeet_div) - 5L, inherits(pop, "Population")))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(SugarBeet_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
