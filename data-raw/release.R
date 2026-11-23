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

# Each hosted catalog row names its own release tag (the `release` column,
# set in data-raw/catalog.R). A row that already has a file/release/bytes/md5
# is left as it is -- never re-packed and never re-tagged -- so the published
# assets of an earlier tag stay valid. Only rows still missing a checksum are
# packed, into data-raw/release/<name>.rds. Set FORCE <- TRUE to re-pack all.
FORCE <- FALSE
DEFAULT_RELEASE <- "data-v1" # used only for a row whose release is empty

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
  if (is.na(catalog$file[i])) next      # sample-only datasets (no hosted file): nothing to release

  name <- catalog$name[i]
  if (!FORCE && !is.na(catalog$md5[i]) && nzchar(catalog$md5[i])) next # already released
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
  if (is.na(catalog$release[i]) || !nzchar(catalog$release[i])) catalog$release[i] <- DEFAULT_RELEASE
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
  cat("Nothing packed (every hosted row already has a checksum, or no\n",
     "data-raw/built/<name>.rda was found).\n", sep = "")
}
if (length(missing_)) {
  cat("\nSkipped (no data-raw/built/<name>.rda found; catalog left unchanged):\n  ",
     paste(missing_, collapse = ", "), "\n", sep = "")
}

# One `gh` command per release tag, listing only that tag's packed assets.
for (tag in unique(catalog$release[catalog$name %in% packed$name])) {
  files <- catalog$file[catalog$name %in% packed$name & catalog$release == tag]
  cat(sprintf("\nTo publish release %s (first time):\n", tag))
  cat(sprintf('  gh release create %s --title "markerdata %s" --notes "Packed marker panels for markerdata" %s\n',
              tag, sub("^data-v", "data v", tag),
              paste(file.path("data-raw/release", files), collapse = " ")))
  cat(sprintf("To add or replace assets in an existing release %s:\n", tag))
  cat(sprintf("  gh release upload %s %s --clobber\n", tag,
              paste(file.path("data-raw/release", files), collapse = " ")))
}
