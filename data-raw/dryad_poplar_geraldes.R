#!/usr/bin/env Rscript
# data-raw/dryad_poplar_geraldes.R
#
# Build a *poplar* genotype dataset for the Breeding Designer / simplePHENOTYPES,
# in the SAME numeric format as simplePHENOTYPES::SNP55K_maize282_maf04 (snp,
# allele, chr, pos, cm, then one -1/0/1 column per accession).
#
# Source: Populus trichocarpa landscape-genomics panel genotyped on the
# PoplarConsortium-ORNL 34K SNP array, deposited (CC0) at
#   Dryad doi:10.5061/dryad.7s848  (FOR_DRYAD_TandB_fwd_june2013_434.txt)
# from
#   Geraldes A., Farzaneh N., Grassa C.J., McKown A.D., Guy R.D., Mansfield S.D.,
#     Douglas C.J., Cronk Q.C.B. (2014) "Landscape genomics of Populus
#     trichocarpa: the role of hybridization, limited gene flow, and natural
#     selection in shaping patterns of population structure", Evolution
#     68(11):3260-3280, doi:10.1111/evo.12497
#   Array: Geraldes A. et al. (2013) "A 34K SNP genotyping array for Populus
#     trichocarpa", Mol. Ecol. Resour. 13(2):306-323, doi:10.1111/1755-0998.12056
#
# The file is an Illumina GenomeStudio final report (CR line endings; an 8-line
# [Header], then [Data]: a marker x sample matrix). Each cell is "<GT>|<GenCall>"
# where <GT> is the forward-strand two-nucleotide genotype (e.g. "AG", "GG", "" =
# missing). Marker names encode position: scaffold_<chr>_<pos>. This strips the
# GenCall score, recodes each biallelic SNP to -1/0/1 (dosage of the alphabetically
# greater allele, minus 1), parses chr/pos from the name, keeps the 19 nuclear
# chromosomes, adds a synthetic map, MAF-filters (>=0.05) and thins evenly to a
# curated 12,000-SNP subset -> data/Poplar_div.rda. Poplar trees are OUTBRED, so a
# dominance layer can be simulated (cf. EUChip60K_div).
#
# NOTE: the raw file and the built .rda are NOT tracked in git (see .gitignore).
#
# Usage:  Rscript data-raw/dryad_poplar_geraldes.R [report.txt] [n_snps] [out.rda]

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else
  "data-raw/poplar/FOR_DRYAD_TandB_fwd_june2013_434.txt"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data/Poplar_div.rda"
stopifnot(file.exists(infile))

## 1. stream: CR->LF, drop the 9 header lines, strip "|<score>" ------------
cmd <- sprintf("tr '\\r' '\\n' < %s | sed '1,9d' | sed 's/|[0-9.]*//g'",
               shQuote(normalizePath(infile)))
dat <- data.table::fread(cmd = cmd, header = TRUE, sep = "\t", data.table = FALSE,
                         colClasses = "character")
marker <- dat[[1]]
Gn <- as.matrix(dat[, -1, drop = FALSE])          # marker x sample, "AG"/"GG"/""
message(sprintf("GenomeStudio matrix: %d markers x %d samples.",
                nrow(Gn), ncol(Gn)))

## 2. recode two-nucleotide genotypes -> -1/0/1 ----------------------------
c1 <- substr(Gn, 1, 1); c2 <- substr(Gn, 2, 2)
c1[c1 == ""] <- NA; c2[c2 == ""] <- NA
dim(c1) <- dim(c2) <- dim(Gn)
df1 <- as.data.frame(c1, stringsAsFactors = FALSE)
df2 <- as.data.frame(c2, stringsAsFactors = FALSE)
a2 <- do.call(pmax, c(df1, df2, list(na.rm = TRUE)))  # greater allele per marker
a1 <- do.call(pmin, c(df1, df2, list(na.rm = TRUE)))  # lesser allele per marker
M  <- (c1 == a2) + (c2 == a2)                     # copies of a2: 0/1/2 (a2 recycles)
M[is.na(c1) | is.na(c2)] <- NA
storage.mode(M) <- "integer"
M <- M - 1L                                       # -1/0/1
colnames(M) <- colnames(Gn)

## 3. parse chr/pos from scaffold_<chr>_<pos>; keep nuclear chromosomes -----
m   <- regmatches(marker, regexec("^scaffold_([0-9]+)_([0-9]+)$", marker))
chr <- vapply(m, function(x) if (length(x) == 3L) as.integer(x[[2]]) else NA_integer_, 0L)
pos <- vapply(m, function(x) if (length(x) == 3L) as.integer(x[[3]]) else NA_integer_, 0L)
biallelic <- !is.na(a1) & !is.na(a2) & a1 != a2
ok  <- !is.na(chr) & chr >= 1 & chr <= 19 & !is.na(pos) & biallelic
G <- data.frame(snp = marker[ok], allele = paste(a1[ok], a2[ok], sep = "/"),
                chr = sprintf("%02d", chr[ok]), pos = pos[ok], cm = NA_real_,
                check.names = FALSE, stringsAsFactors = FALSE)
M <- M[ok, , drop = FALSE]
o <- order(G$chr, G$pos)                          # array order != position order
G <- G[o, , drop = FALSE]; M <- M[o, , drop = FALSE]
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

Poplar_div <- data.frame(G[ord, 1:5, drop = FALSE],
                         as.data.frame(M[ord, , drop = FALSE], check.names = FALSE),
                         check.names = FALSE, stringsAsFactors = FALSE)
rownames(Poplar_div) <- NULL

## 6. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(Poplar_div[, 6:ncol(Poplar_div)]) %in% c(-1L, 0L, 1L)))
pop <- as_population(Poplar_div, individuals = seq_len(min(10L, ncol(Poplar_div) - 5L)))
gt  <- table(unlist(Poplar_div[, 6:ncol(Poplar_div)]))
message(sprintf("Built Poplar_div: %d SNPs x %d accessions (as_population OK: %s).",
                nrow(Poplar_div), ncol(Poplar_div) - 5L, inherits(pop, "Population")))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Poplar_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
