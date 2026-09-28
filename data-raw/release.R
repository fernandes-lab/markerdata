#!/usr/bin/env Rscript
# data-raw/release.R
#
# Pack the built marker objects (data-raw/built/<name>.rda) into markerdata's
# packed release format (see data-raw/README.md for the format), stage them
# under data-raw/release/ for a GitHub release, and update
# inst/extdata/catalog.csv with the resulting file name, release tag, size
# and checksum for each dataset that was packed.
#
# Run from the package root, after building the datasets you want to release
# (see the data-raw/*.R build scripts) and once the package's internal
# .pack() / .unpack() (R/pack.R) and .check_marker() helpers are available:
#
#   Rscript data-raw/release.R
#
# This script does NOT upload anything itself -- it prints the exact `gh`
# commands to run afterward (see the end of this script / its console output).

RELEASE <- "data-v1"

suppressMessages(pkgload::load_all(".", quiet = TRUE))

catalog_path <- "inst/extdata/catalog.csv"
stopifnot(file.exists(catalog_path))
catalog <- utils::read.csv(catalog_path, stringsAsFactors = FALSE, na.strings = "")
catalog$file    <- as.character(catalog$file)
catalog$release <- as.character(catalog$release)
catalog$bytes   <- as.character(catalog$bytes)
catalog$md5     <- as.character(catalog$md5)

dir.create("data-raw/release", recursive = TRUE, showWarnings = FALSE)

packed   <- data.frame(name = character(), bytes = double(), md5 = character(),
                       stringsAsFactors = FALSE)
missing_ <- character()

for (i in seq_len(nrow(catalog))) {
  if (isTRUE(catalog$bundled[i])) next  # bundled datasets are not released as assets

  name <- catalog$name[i]
  rda  <- file.path("data-raw", "built", paste0(name, ".rda"))
  if (!file.exists(rda)) {
    missing_ <- c(missing_, name)
    message("release.R: ", rda, " not found; leaving ", name, " unchanged in the catalog.")
    next
  }

  e  <- new.env()
  nm <- load(rda, envir = e)
  stopifnot(identical(nm[1L], name))
  x <- get(nm[1L], envir = e)

  .check_marker(x)
  p <- .pack(x, name)
  stopifnot(identical(.unpack(p), x))

  out_rds <- file.path("data-raw", "release", paste0(name, ".rds"))
  saveRDS(p, file = out_rds, compress = "xz")

  bytes <- file.size(out_rds)
  md5   <- unname(tools::md5sum(out_rds))

  catalog$file[i]    <- paste0(name, ".rds")
  catalog$release[i] <- RELEASE
  catalog$bytes[i]   <- as.character(bytes)
  catalog$md5[i]     <- md5

  packed <- rbind(packed, data.frame(name = name, bytes = bytes, md5 = md5,
                                     stringsAsFactors = FALSE))
  message(sprintf("Packed %s -> %s (%.0f KB, md5 %s).",
                  name, out_rds, bytes / 1024, md5))
}

utils::write.csv(catalog, file = catalog_path, row.names = FALSE, na = "")
message("Updated ", catalog_path, ".")

cat("\n")
if (nrow(packed)) {
  base::print(packed, row.names = FALSE)
  cat(sprintf("\nTotal size: %.2f MB (%d dataset%s packed)\n",
             sum(packed$bytes) / 1024^2, nrow(packed),
             if (nrow(packed) != 1L) "s" else ""))
} else {
  cat("Nothing packed (no data-raw/built/<name>.rda files were found for the\n",
     "non-bundled catalog rows).\n", sep = "")
}
if (length(missing_)) {
  cat("\nSkipped (no data-raw/built/<name>.rda found; catalog left unchanged):\n  ",
     paste(missing_, collapse = ", "), "\n", sep = "")
}

cat("\nTo publish this release (first time):\n")
cat(sprintf(
  '  gh release create %s --title "markerdata data v1" --notes "Packed marker panels for markerdata 0.1.0" data-raw/release/*.rds\n',
  RELEASE))
cat("\nTo update an existing release with newly packed assets:\n")
cat(sprintf("  gh release upload %s data-raw/release/*.rds --clobber\n", RELEASE))
