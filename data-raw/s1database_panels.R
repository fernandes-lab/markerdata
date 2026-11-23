#!/usr/bin/env Rscript
# data-raw/s1database_panels.R
#
# Build the two sample-only data sets, Arabidopsis_div and Human_div, from the
# already-processed pilot panels of the s1database project (6,000 markers
# each). markerdata does NOT host a whole panel for them: only the 2,000 SNP
# x all-individuals sample ships in inst/extdata/samples/<name>.rds, and the
# catalog's `source_url` points to the original full download.
#
#   Arabidopsis_div  1001 Genomes (1,135 accessions; 1,097 in the pilot panel),
#                    TAIR10, https://1001genomes.org/data/GMI-MPI/releases/v3.1/
#   Human_div        1000 Genomes Project phase 3 (2,504 people), GRCh37,
#                    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#
# Inputs (read-only; per panel <dir>/pilot-v1/): variants.tsv,
# dosages.simulation.tsv.gz (rows = samples, columns = markers, -1/0/1; the
# imputed file), manifest.json. The raw VCFs are never touched.
#
# Steps: read and check -> build the marker data frame (snp = marker_id,
# allele = "REF/ALT", chr, pos, synthetic cm) sorted by chromosome and
# position -> thin with the rule in data-raw/thin.R (same as samples.R) ->
# .check_marker() -> .pack() -> inst/extdata/samples/<name>.rds.
#
# Run from the package root:  Rscript data-raw/s1database_panels.R [panels_dir]
# (default panels_dir: the s1database/data-panels folder next to markerdata)

args <- commandArgs(trailingOnly = TRUE)
panels <- if (length(args)) args[1L] else
  normalizePath("../s1database/data-panels", mustWork = TRUE)

suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(simplePHENOTYPES)) # synthetic_map(), as in the other build scripts
source("data-raw/thin.R") # N_SNP, N_IND, thin_sample()

build <- function(name, dir) {
  f <- file.path(panels, dir, "pilot-v1")
  v <- utils::read.delim(file.path(f, "variants.tsv"), colClasses = "character")
  gz <- file.path(f, "dosages.simulation.tsv.gz")
  hdr <- names(utils::read.delim(gzfile(gz), nrows = 1L, check.names = FALSE))
  d <- utils::read.table(gzfile(gz), header = TRUE, sep = "\t", quote = "",
                         comment.char = "", check.names = FALSE,
                         colClasses = c("character", rep("integer", length(hdr) - 1L)))
  ids <- d[[1L]]
  M <- t(as.matrix(d[-1L])) # markers x individuals
  stopifnot(identical(rownames(M), v$marker_id), !anyDuplicated(ids), all(nzchar(ids)))
  if (anyNA(M)) {
    # shortcut: not triggered for the pilot-v1 files (checked: 0 NA); impute to the major homozygote
    message(name, ": ", sum(is.na(M)), " NA in the simulation file; imputing to the major homozygote.")
    maj <- ifelse(rowSums(M == 1L, na.rm = TRUE) > rowSums(M == -1L, na.rm = TRUE), 1L, -1L)
    M[is.na(M)] <- maj[row(M)[is.na(M)]]
  }
  storage.mode(M) <- "integer"
  stopifnot(all(M %in% -1:1))

  x <- data.frame(snp = v$marker_id, allele = paste0(v$ref, "/", v$alt),
                  chr = v$chromosome, pos = as.integer(v$position),
                  cm = NA_real_, stringsAsFactors = FALSE)
  x <- cbind(x, as.data.frame(M, check.names = FALSE))
  colnames(x)[-(1:5)] <- ids # ids stay character names even when numeric
  x <- x[order(suppressWarnings(as.integer(x$chr)), x$chr, x$pos), ]
  x$cm <- synthetic_map(x$chr, x$pos)
  rownames(x) <- NULL

  s <- thin_sample(x)
  .check_marker(s, name)
  p <- .pack(s, name)
  stopifnot(identical(.unpack(p), s))
  out <- file.path("inst", "extdata", "samples", paste0(name, ".rds"))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(p, out, compress = "xz")
  message(sprintf("%s: %d SNPs x %d individuals (from %d x %d) -> %s (%.0f KB)",
                  name, nrow(s), ncol(s) - 5L, nrow(x), ncol(x) - 5L, out, file.size(out) / 1024))
}

build("Arabidopsis_div", "arabidopsis_1001g")
build("Human_div", "human_1kg_phase3")
