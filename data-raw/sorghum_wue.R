#!/usr/bin/env Rscript
# data-raw/sorghum_wue.R
#
# Build a *sorghum* genotype dataset for the Breeding Designer / simplePHENOTYPES,
# in the SAME numeric format as simplePHENOTYPES::SNP55K_maize282_maf04 and
# breedingDesigner::SoySNP50K_div / EUChip60K_div (a data.frame whose first five
# columns are snp, allele, chr, pos, cm, followed by one -1/0/1 column per
# individual).
#
# Source: the sorghum WEST Diversity Panel genotyped by GBS, deposited (CC-BY-4.0,
# redistributable with attribution) at Zenodo doi:10.5281/zenodo.5019227
# ("WEST Diversity Panel GBS information", rawGBS_imputed_850x5mi.vcf.gz: 5.5M SNPs
# x 850 [842 unique] individuals). Used in the water-use-efficiency GWAS/TWAS study
#   Ferguson J.N., Fernandes S.B., Monier B., Miller N.D., Allen D., Dmitrieva A.,
#     Schmuker P., Lozano R., Valluru R., Buckler E.S., Gore M.A., Brown P.J.,
#     Spalding E.P., Leakey A.D.B. (2021) "Machine learning-enabled phenotyping for
#     GWAS and TWAS of WUE traits in 869 field-grown sorghum accessions",
#     Plant Physiology 187(3):1481-1500, doi:10.1093/plphys/kiab346
# The input markers.txt below is a recoded -1/0/1 genotype matrix derived from that
# deposit (842 accessions x ~499k SNPs).
#
# markers.txt layout: a SNPs x (meta + individuals) matrix, TAB-separated. Header:
#   snp  allele  chr  pos  cm  <indiv_1> <indiv_2> ...   (snp and cm are EMPTY;
#   allele like "A/T"; chr 1-10; pos in bp; individual columns named NAME_NAME).
# Genotypes are already coded -1/0/1 (ref-hom / het / alt-hom) with no missing.
#
# This script (1) streams the big file with awk to compute per-SNP MAF without
# loading it into R, (2) thins the MAF>=0.05 SNPs evenly across the 10 sorghum
# chromosomes to n_snps, (3) extracts just those rows, (4) synthesizes SNP IDs
# (S<chr>_<pos>) and a synthetic genetic map, cleans the doubled accession names,
# and saves data/Sorghum_div.rda.
#
# NOTE: the raw markers.txt and the built .rda are NOT tracked in git yet (see
# .gitignore) -- the example marker datasets are slated to move to a separate
# companion data package (breedingDesignerData) so breedingDesigner stays small.
#
# Usage:  Rscript data-raw/sorghum_wue.R [markers.txt] [n_snps] [out.rda]
#   defaults: the collaborator's markers.txt, 12000 SNPs -> data/Sorghum_div.rda,
#             all accessions.

suppressMessages(library(simplePHENOTYPES))
args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else
  "/Users/samuelbf/Library/CloudStorage/OneDrive-UniversityofArkansas/UARK/collaboration/Danilo/ML_haplotype_gwas/markers.txt"
n_snps <- if (length(args) >= 2) as.integer(args[[2]]) else 12000L
out    <- if (length(args) >= 3) args[[3]] else "data/Sorghum_div.rda"
stopifnot(file.exists(infile))

tmp     <- tempdir()
maf_tsv <- file.path(tmp, "sorghum_maf.tsv")     # NR chr pos allele maf (MAF>=0.05)
keep_tx <- file.path(tmp, "sorghum_keep.txt")    # selected markers.txt line numbers
sub_tsv <- file.path(tmp, "sorghum_subset.tsv")  # header + selected rows only

## 1. streaming MAF pass (awk) ---------------------------------------------
# For each data row: p = (mean(g)+1)/2 from the -1/0/1 dosages in cols 6..NF;
# emit line number + locus info when MAF >= 0.05. The awk program is written to a
# file (-f) so tab/newline escapes reach awk verbatim, not the shell or R.
message("Pass A: streaming MAF over ", infile, " ...")
progA <- file.path(tmp, "mafA.awk")
writeLines(c(
  'BEGIN { FS="\\t" }',
  'NR==1 { next }',
  '{',
  '  s=0; n=0;',
  '  for (i=6; i<=NF; i++) { s+=$i; n++; }',
  '  if (n==0) next;',
  '  p=(s/n + 1)/2; maf=(p<1-p)?p:1-p;',
  '  if (maf>=0.05) printf "%d\\t%s\\t%s\\t%s\\t%.4f\\n", NR, $3, $4, $2, maf;',
  '}'
), progA)
system2("awk", c("-f", shQuote(progA), shQuote(infile)), stdout = maf_tsv)
maf <- data.table::fread(maf_tsv, header = FALSE, sep = "\t", data.table = FALSE,
                         col.names = c("nr", "chr", "pos", "allele", "maf"))
message(sprintf("  SNPs with MAF >= 0.05: %d (of %d chromosomes).",
                nrow(maf), length(unique(maf$chr))))

## 2. thin evenly across the genome to n_snps ------------------------------
ord <- order(maf$chr, maf$pos)
if (length(ord) > n_snps) ord <- ord[round(seq(1, length(ord), length.out = n_snps))]
sel <- maf[ord, , drop = FALSE]
writeLines(as.character(sel$nr), keep_tx)
message(sprintf("  Selected %d SNPs, spanning chromosomes: %s.",
                nrow(sel), paste(sort(unique(sel$chr)), collapse = ",")))

## 3. extract just the selected rows (awk) ---------------------------------
message("Pass B: extracting selected rows ...")
progB <- file.path(tmp, "subB.awk")
writeLines(c(
  'BEGIN { FS="\\t" }',
  'NR==FNR { keep[$1]=1; next }',   # first file: the selected line numbers
  'FNR==1 { print; next }',          # markers.txt header
  '(FNR in keep) { print }'          # selected SNP rows ($0 kept verbatim)
), progB)
system2("awk", c("-f", shQuote(progB), shQuote(keep_tx), shQuote(infile)),
        stdout = sub_tsv)
sub <- data.table::fread(sub_tsv, header = TRUE, sep = "\t", data.table = FALSE,
                         check.names = FALSE)
message(sprintf("  Read subset: %d SNPs x %d columns.", nrow(sub), ncol(sub)))

## 4. assemble the standard 5-meta + individuals data.frame ----------------
acc <- sub(pattern = "_.*$", replacement = "", x = colnames(sub)[6:ncol(sub)])
if (anyDuplicated(acc)) acc <- make.unique(acc)     # guard doubled-name collisions
chr <- as.integer(sub$chr)
pos <- as.integer(sub$pos)
Sorghum_div <- data.frame(
  snp    = sprintf("S%d_%d", chr, pos),             # the file carries no SNP IDs
  allele = as.character(sub$allele),
  chr    = sprintf("%02d", chr),                    # canonical chromosome label
  pos    = pos,
  cm     = NA_real_,
  stringsAsFactors = FALSE, check.names = FALSE
)
if (anyDuplicated(Sorghum_div$snp))
  Sorghum_div$snp <- make.unique(Sorghum_div$snp)
Sorghum_div$cm <- synthetic_map(Sorghum_div$chr, Sorghum_div$pos)
G <- as.matrix(sub[, 6:ncol(sub), drop = FALSE]); storage.mode(G) <- "integer"
colnames(G) <- acc
Sorghum_div <- cbind(Sorghum_div, as.data.frame(G, check.names = FALSE))

## 5. sanity check + save ---------------------------------------------------
stopifnot(all(unlist(Sorghum_div[, 6:ncol(Sorghum_div)]) %in% c(-1L, 0L, 1L)))
pop <- as_population(Sorghum_div,
                     individuals = seq_len(min(10L, ncol(Sorghum_div) - 5L)))
gt  <- table(unlist(Sorghum_div[, 6:ncol(Sorghum_div)]))
message(sprintf("Built Sorghum_div: %d SNPs x %d accessions (as_population OK: %s).",
                nrow(Sorghum_div), ncol(Sorghum_div) - 5L, inherits(pop, "Population")))
message("Genotype codes: ", paste(names(gt), gt, sep = "=", collapse = "  "),
        sprintf("  (het fraction %.3f)", gt[["0"]] / sum(gt)))
save(Sorghum_div, file = out, compress = "xz")
message(sprintf("Saved -> %s (%.0f KB).", normalizePath(out), file.size(out) / 1024))
