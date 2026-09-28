#!/usr/bin/env Rscript
# data-raw/dryad_euchip60k.R
#
# Build the *eucalyptus* marker dataset (EUChip60K_div) for the markerdata
# package, in the SAME numeric format as
# simplePHENOTYPES::SNP55K_maize282_maf04 (a data.frame whose first five
# columns are snp, allele, chr, pos, cm, followed by one -1/0/1 column per
# individual). markerdata does not ship this object; it publishes the packed
# result as a GitHub release asset that marker_data() downloads and caches on
# first use (see data-raw/release.R).
#
# Source: Eucalyptus grandis x E. urophylla synthetic breeding population,
# genotyped on the EUChip60K Illumina Infinium array (Silva-Junior et al. 2015,
# New Phytologist 206:1527), 24,806 QC SNPs x 768 OUTBRED trees. Deposited with:
#   Resende R.T., Resende M.D.V., Silva F.F. et al. (2017)
#     "Assessing the expected response to genomic selection of individuals and
#      families in Eucalyptus breeding with an additive-dominant model",
#     Heredity 119:245-255, doi:10.1038/hdy.2017.37
#   Data: Dryad Digital Repository, doi:10.5061/dryad.ms580  (license CC0-1.0,
#         public domain -- redistributable). File: SNP_Data.txt.
# License: CC0-1.0 (public domain, Dryad); cite the papers above.
# Unlike inbred soybean, these trees are OUTBRED (heterozygous loci present), so
# a dominance layer can be simulated -- useful for schemes that exploit
# heterosis / combining ability, such as reciprocal recurrent selection.
#
# SNP_Data.txt layout: an individuals x SNPs matrix, space-separated, header row
# = SNP IDs "EuBR{chr}s{pos}A..B..C..D..E.." (chr 01-11, pos in bp); the first
# field of each data row is the tree ID; genotypes coded 0/1/2 (allele dosage).
#
# This script transposes to SNPs x individuals, recodes 0/1/2 -> -1/0/1, parses
# chr/pos from the SNP IDs, adds a synthetic genetic map (the array carries no
# cM), sanity-checks with as_population(), and saves a curated MAF-filtered
# subset as data-raw/built/EUChip60K_div.rda (mirroring
# data-raw/soybase_soysnp50k.R).
#
# Usage:  Rscript data-raw/dryad_euchip60k.R \
#           [data-raw/eucalyptus/SNP_Data.txt] [n_snps] [out.rda]
#         (defaults: 12000 SNPs -> data-raw/built/EUChip60K_div.rda, all 768 trees)

suppressMessages(library(simplePHENOTYPES))
args  <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else "data-raw/eucalyptus/SNP_Data.txt"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data-raw/built/EUChip60K_div.rda"

stopifnot(file.exists(infile))

## 1. read -----------------------------------------------------------------
# Header = SNP IDs (24,806); each data row = tree_id + 24,806 dosages.
snp_ids <- scan(infile, what = character(), nlines = 1L, quiet = TRUE)
have_dt <- requireNamespace("data.table", quietly = TRUE)
if (have_dt) {
  message("Reading with data.table::fread ...")
  raw <- data.table::fread(infile, skip = 1L, header = FALSE,
                           sep = " ", data.table = FALSE)
} else {
  message("Reading with read.table (data.table not installed) ...")
  raw <- utils::read.table(infile, skip = 1L, header = FALSE,
                           sep = " ", check.names = FALSE,
                           stringsAsFactors = FALSE)
}
tree_ids <- as.character(raw[[1]])
G <- as.matrix(raw[, -1, drop = FALSE])          # individuals x SNPs, 0/1/2
storage.mode(G) <- "integer"
stopifnot(ncol(G) == length(snp_ids))
message(sprintf("Read %d trees x %d SNPs (codes: %s).",
                nrow(G), ncol(G), paste(sort(unique(as.vector(G[1:5, 1:50]))),
                                        collapse = "/")))

## 2. transpose + recode 0/1/2 -> -1/0/1 -----------------------------------
Gt <- t(G) - 1L                                  # SNPs x individuals, -1/0/1
colnames(Gt) <- make.unique(tree_ids)
rm(G, raw); invisible(gc())

## 3. parse chr / pos from the SNP IDs -------------------------------------
m   <- regmatches(snp_ids, regexec("^EuBR([0-9]+)s([0-9]+)", snp_ids))
chr <- vapply(m, function(x) if (length(x) == 3L) x[[2]] else NA_character_, "")
pos <- vapply(m, function(x) if (length(x) == 3L) as.integer(x[[3]]) else NA_integer_,
              integer(1))
ok  <- !is.na(chr) & !is.na(pos)
if (!all(ok)) message(sprintf("Dropping %d SNP(s) with unparseable IDs.", sum(!ok)))
snp_ids <- snp_ids[ok]; chr <- chr[ok]; pos <- pos[ok]; Gt <- Gt[ok, , drop = FALSE]
chr <- sprintf("EuBR%02d", as.integer(chr))       # canonical linkage-group label

## 4. MAF filter + genome-wide thin to n_snps ------------------------------
p   <- (rowMeans(Gt) + 1) / 2                      # allele freq from -1/0/1
maf <- pmin(p, 1 - p)
keep <- which(maf >= 0.05)
message(sprintf("SNPs with MAF >= 0.05: %d / %d.", length(keep), length(maf)))
# order by (chr, pos) then thin evenly so the subset spans all 11 chromosomes
ord  <- keep[order(chr[keep], pos[keep])]
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]
sel  <- ord

EUChip60K_div <- data.frame(
  snp    = snp_ids[sel],
  allele = "A/G",                                 # array alleles not in the file
  chr    = chr[sel],
  pos    = pos[sel],
  cm     = NA_real_,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
EUChip60K_div$cm <- synthetic_map(EUChip60K_div$chr, EUChip60K_div$pos)
EUChip60K_div <- cbind(EUChip60K_div, as.data.frame(Gt[sel, , drop = FALSE],
                                                    check.names = FALSE))

## 5. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(EUChip60K_div[, 6:ncol(EUChip60K_div)]) %in% c(-1L, 0L, 1L)))
pop <- as_population(EUChip60K_div,
                     individuals = seq_len(min(10L, ncol(EUChip60K_div) - 5L)))
gt  <- table(unlist(EUChip60K_div[, 6:ncol(EUChip60K_div)]))
message(sprintf("Built EUChip60K_div: %d SNPs x %d trees (as_population OK: %s).",
                nrow(EUChip60K_div), ncol(EUChip60K_div) - 5L,
                inherits(pop, "Population")))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)",
                gt[["0"]] / sum(gt)))
save(EUChip60K_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out),
                file.size(out) / 1024))
