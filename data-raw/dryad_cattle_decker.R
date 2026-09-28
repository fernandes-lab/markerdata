#!/usr/bin/env Rscript
# data-raw/dryad_cattle_decker.R
#
# Build the *cattle* marker dataset (Cattle_div) for the markerdata package,
# in the SAME numeric format as simplePHENOTYPES::SNP55K_maize282_maf04 (snp,
# allele, chr, pos, cm, then one -1/0/1 column per animal). markerdata does
# not ship this object; it publishes the packed result as a GitHub release
# asset that marker_data() downloads and caches on first use (see
# data-raw/release.R).
#
# Source: worldwide cattle diversity panel (134 breeds), Illumina BovineSNP50,
# deposited (CC0) at
#   Dryad doi:10.5061/dryad.th092  (WorldWideCattleDryad.zip -> PLINK .ped/.map)
# from
#   Decker J.E., McKay S.D., Rolf M.M., Kim J., Molina Alcala A., Sonstegard T.S.,
#     Hanotte O., ... Taylor J.F. (2014) "Worldwide Patterns of Ancestry,
#     Divergence, and Admixture in Domesticated Cattle", PLoS Genetics
#     10(3):e1004254, doi:10.1371/journal.pgen.1004254
# License: CC0 (public domain, Dryad); cite the paper above.
#
# PLINK .ped alleles are coded 1 (Illumina A) / 2 (Illumina B), 0 = missing; the
# .map gives chr (UMD3.1) + bp. This recodes allele dosage to -1/0/1 (count of the
# B allele, minus 1), drops the unassigned chromosome 0, adds a synthetic genetic
# map, MAF-filters (>=0.05) and thins evenly across the 29 autosomes to a curated
# 12,000-SNP subset -> data-raw/built/Cattle_div.rda. Cattle are OUTBRED, so
# heterozygosity is present (a dominance layer can be simulated, cf. EUChip60K_div).
#
# Usage:  Rscript data-raw/dryad_cattle_decker.R [WorldWideCattleDryad.zip] [n_snps] [out.rda]
#   defaults: data-raw/cattle/WorldWideCattleDryad.zip, 12000 SNPs ->
#             data-raw/built/Cattle_div.rda.

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else "data-raw/cattle/WorldWideCattleDryad.zip"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/Cattle_div.rda"
stopifnot(file.exists(infile))

tmp <- file.path(tempdir(), "cattle"); dir.create(tmp, showWarnings = FALSE)

## 1. unzip + locate .map / .ped -------------------------------------------
if (grepl("\\.zip$", infile)) utils::unzip(infile, exdir = tmp) else
  file.copy(infile, tmp)
mapf <- list.files(tmp, "\\.map$", full.names = TRUE, recursive = TRUE)[1]
pedf <- list.files(tmp, "\\.ped$", full.names = TRUE, recursive = TRUE)[1]
stopifnot(file.exists(mapf), file.exists(pedf))

map <- data.table::fread(mapf, header = FALSE, data.table = FALSE)  # chr snp cM bp
message(sprintf("Map: %d SNPs; chromosomes %s.", nrow(map),
                paste(range(as.integer(map$V1)), collapse = "-")))

## 2. read .ped (1/2/0 alleles) + recode to -1/0/1 dosage ------------------
ped <- data.table::fread(pedf, header = FALSE, sep = " ", data.table = FALSE)
iid <- as.character(ped[[2]])
a   <- as.matrix(ped[, 7:ncol(ped), drop = FALSE]); storage.mode(a) <- "integer"
rm(ped); invisible(gc())
stopifnot(ncol(a) == 2L * nrow(map))
odd <- seq(1L, ncol(a), by = 2L)
A1 <- a[, odd, drop = FALSE]; A2 <- a[, odd + 1L, drop = FALSE]  # animals x SNP
rm(a); invisible(gc())
B <- (A1 == 2L) + (A2 == 2L)                      # copies of the B allele: 0/1/2
B[A1 == 0L | A2 == 0L] <- NA                      # 0 = missing call
storage.mode(B) <- "integer"
rm(A1, A2); invisible(gc())
M <- t(B - 1L)                                    # SNP x animal, -1/0/1
colnames(M) <- iid
rm(B); invisible(gc())

## 3. assemble + drop unassigned chr 0 -------------------------------------
G <- data.frame(snp = map$V2, allele = "A/B", chr = as.integer(map$V1),
                pos = as.integer(map$V4), cm = NA_real_,
                check.names = FALSE, stringsAsFactors = FALSE)
ok <- G$chr >= 1 & G$chr <= 29 & !is.na(G$pos)
G <- G[ok, , drop = FALSE]; M <- M[ok, , drop = FALSE]

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
G$chr <- sprintf("%02d", G$chr)
G$cm  <- synthetic_map(G$chr, G$pos)
p   <- (rowMeans(M) + 1) / 2
maf <- pmin(p, 1 - p)
keep <- which(maf >= 0.05)
message(sprintf("SNPs with MAF >= 0.05: %d / %d.", length(keep), length(maf)))
ord <- keep[order(G$chr[keep], G$pos[keep])]
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]

Cattle_div <- data.frame(G[ord, 1:5, drop = FALSE],
                         as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                         check.names = FALSE, stringsAsFactors = FALSE)
rownames(Cattle_div) <- NULL

## 6. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(Cattle_div[, 6:ncol(Cattle_div)]) %in% c(-1L, 0L, 1L)))
pop <- as_population(Cattle_div, individuals = seq_len(min(10L, ncol(Cattle_div) - 5L)))
gt  <- table(unlist(Cattle_div[, 6:ncol(Cattle_div)]))
message(sprintf("Built Cattle_div: %d SNPs x %d animals (as_population OK: %s).",
                nrow(Cattle_div), ncol(Cattle_div) - 5L, inherits(pop, "Population")))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Cattle_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
