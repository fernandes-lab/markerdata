#!/usr/bin/env Rscript
# data-raw/dryad_rice_spindel.R
#
# Build the *rice* marker dataset (Rice_div) for the markerdata package, in
# the SAME numeric format as simplePHENOTYPES::SNP55K_maize282_maf04 (a data.frame
# whose first five columns are snp, allele, chr, pos, cm, followed by one -1/0/1
# column per accession). markerdata does not ship this object; it publishes the
# packed result as a GitHub release asset that marker_data() downloads and
# caches on first use (see data-raw/release.R).
#
# Source: the IRRI elite tropical rice genomic-selection breeding panel (indica /
# indica-admixed), genotyped by GBS, deposited (CC0) at
#   Dryad doi:10.5061/dryad.7369p  (file MET_crfilt_.75_allchrom.hmp.txt.zip)
# from
#   Spindel J., Begum H., Akdemir D., Virk P., Collard B., Redona E., Atlin G.,
#     Jannink J-L., McCouch S.R. (2015) "Genomic Selection and Association Mapping
#     in Rice (Oryza sativa): Effect of Trait Genetic Architecture, Training
#     Population Composition, Marker Number and Statistical Model on Accuracy of
#     Rice Genomic Selection in Elite, Tropical Rice Breeding Lines",
#     PLoS Genetics 11(2):e1004982, doi:10.1371/journal.pgen.1004982
# License: CC0 (public domain, Dryad); cite the paper above.
#
# The deposited HapMap (rs#/alleles/chrom/pos + IUPAC single-char genotypes) is
# converted to numeric -1/0/1 with simplePHENOTYPES::as_numeric(from="hapmap"),
# given a synthetic genetic map, MAF-filtered (>=0.05) and thinned evenly across
# the 12 rice chromosomes to a curated 12,000-SNP subset ->
# data-raw/built/Rice_div.rda.
#
# Usage:  Rscript data-raw/dryad_rice_spindel.R [hapmap.zip|.hmp.txt] [n_snps] [out.rda]
#   defaults: data-raw/rice/MET_crfilt_.75_allchrom.hmp.txt.zip, 12000 SNPs ->
#             data-raw/built/Rice_div.rda.

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else
  "data-raw/rice/MET_crfilt_.75_allchrom.hmp.txt.zip"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/Rice_div.rda"
stopifnot(file.exists(infile))

tmp <- tempdir()

## 1. get the .hmp.txt (unzip if needed) -----------------------------------
if (grepl("\\.zip$", infile)) {
  nm  <- utils::unzip(infile, list = TRUE)$Name
  hmp <- grep("\\.hmp\\.txt$", nm, value = TRUE)[1]
  utils::unzip(infile, files = hmp, exdir = tmp)
  hmp <- file.path(tmp, hmp)
} else hmp <- infile
message("HapMap: ", hmp, " (", round(file.size(hmp)/1e6), " MB)")

## 2. read HapMap + recode single-char IUPAC calls -> -1/0/1 ----------------
# These inbred rice lines are reported as single-character homozygous calls
# (e.g. "T" == TT), which as_numeric(from="hapmap") (expecting diploid two-char)
# cannot parse, so recode directly against the per-SNP `alleles` column.
h <- data.table::fread(hmp, header = TRUE, sep = "\t", data.table = FALSE,
                       colClasses = "character")
message(sprintf("HapMap: %d SNPs x %d accessions.", nrow(h), ncol(h) - 11L))
a1 <- substr(h$alleles, 1, 1)                    # e.g. "T/G" -> a1="T", a2="G"
a2 <- substr(h$alleles, 3, 3)
iupac_het <- function(x, y) {
  key <- paste0(pmin(x, y), pmax(x, y))
  unname(c(AC = "M", AG = "R", AT = "W", CG = "S", CT = "Y", GT = "K")[key])
}
het <- iupac_het(a1, a2)                          # ambiguity code for {a1,a2}
Gc  <- as.matrix(h[, 12:ncol(h), drop = FALSE])   # SNP x accession character calls
M   <- matrix(NA_integer_, nrow(Gc), ncol(Gc), dimnames = list(NULL, colnames(Gc)))
# a1/a2/het are length-nrow and recycle down each column (row-wise comparison):
M[Gc == a1]  <- -1L
M[Gc == a2]  <-  1L
M[Gc == het] <-  0L
G <- data.frame(snp = h[[1]], allele = h$alleles, chr = h$chrom,
                pos = suppressWarnings(as.integer(h$pos)), cm = NA_real_,
                check.names = FALSE, stringsAsFactors = FALSE)
ok <- !is.na(G$pos) & !is.na(a1) & !is.na(a2) & a1 != a2
G  <- G[ok, , drop = FALSE]; M <- M[ok, , drop = FALSE]

## 3. impute residual missing (rounded per-SNP mean) -----------------------
storage.mode(M) <- "double"
na_frac <- mean(is.na(M))
if (na_frac > 0) {
  message(sprintf("Imputing %.3f%% missing/uncalled to the rounded per-SNP mean.",
                  100 * na_frac))
  rmn <- rowMeans(M, na.rm = TRUE); rmn[is.na(rmn)] <- 0
  ix <- which(is.na(M), arr.ind = TRUE)
  M[ix] <- pmax(-1, pmin(1, round(rmn[ix[, 1]])))
}
storage.mode(M) <- "integer"

## 4. synthetic map + MAF filter + genome-wide thin ------------------------
G$cm <- synthetic_map(G$chr, G$pos)
p   <- (rowMeans(M) + 1) / 2
maf <- pmin(p, 1 - p)
keep <- which(maf >= 0.05)
message(sprintf("SNPs with MAF >= 0.05: %d / %d (chromosomes: %s).",
                length(keep), length(maf),
                paste(sort(unique(G$chr)), collapse = ",")))
ord <- keep[order(G$chr[keep], G$pos[keep])]
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]

Rice_div <- data.frame(G[ord, 1:5, drop = FALSE],
                       as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                       check.names = FALSE, stringsAsFactors = FALSE)
rownames(Rice_div) <- NULL
names(Rice_div)[1:5] <- c("snp", "allele", "chr", "pos", "cm")

## 5. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(Rice_div[, 6:ncol(Rice_div)]) %in% c(-1L, 0L, 1L)))
pop <- as_population(Rice_div, individuals = seq_len(min(10L, ncol(Rice_div) - 5L)))
gt  <- table(unlist(Rice_div[, 6:ncol(Rice_div)]))
message(sprintf("Built Rice_div: %d SNPs x %d accessions (as_population OK: %s).",
                nrow(Rice_div), ncol(Rice_div) - 5L, inherits(pop, "Population")))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Rice_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
