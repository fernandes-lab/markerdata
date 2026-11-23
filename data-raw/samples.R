#!/usr/bin/env Rscript
# data-raw/samples.R
#
# Build the small sample of every remote panel that ships inside the package
# (inst/extdata/samples/<name>.rds, same packed format as the release
# assets). marker_data("<name>") returns this sample; marker_data("<name>",
# full = TRUE) downloads the whole panel.
#
# A sample is a deterministic thinning of the built object in
# data-raw/built/<name>.rda (no random numbers, so it is reproducible on any
# R version):
#   1. keep all individuals (N_IND in thin.R; Inf = all);
#   2. drop SNPs that are monomorphic among those individuals;
#   3. keep N_SNP of the remaining SNPs, evenly spaced over the rows, so every
#      chromosome keeps its share and `cm` stays sorted.
#
# The rule itself lives in data-raw/thin.R. Catalog rows with no built .rda
# (Arabidopsis_div, Human_div) are skipped here: their samples come from
# data-raw/s1database_panels.R.
#
# Run from the package root:  Rscript data-raw/samples.R

suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("data-raw/thin.R") # N_SNP, N_IND, thin_sample()

catalog <- utils::read.csv("inst/extdata/catalog.csv", stringsAsFactors = FALSE,
                           na.strings = "")
out_dir <- "inst/extdata/samples"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (name in catalog$name[!catalog$bundled]) {
  rda <- file.path("data-raw", "built", paste0(name, ".rda"))
  if (!file.exists(rda)) {
    message("samples.R: ", rda, " not found; skipping ", name, ".")
    next
  }
  e <- new.env()
  x <- get(load(rda, envir = e)[1L], envir = e)

  s <- thin_sample(x, name)

  .check_marker(s, name)
  p <- .pack(s, name)
  stopifnot(identical(.unpack(p), s))
  out <- file.path(out_dir, paste0(name, ".rds"))
  saveRDS(p, out, compress = "xz")
  message(sprintf("%s: %d SNPs x %d individuals -> %s (%.0f KB)",
                  name, nrow(s), ncol(s) - 5L, out, file.size(out) / 1024))
}
