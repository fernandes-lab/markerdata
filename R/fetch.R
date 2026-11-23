#' Build the remote URL for a catalog entry
#'
#' @param entry One row of the catalog (as returned by `.catalog()` or
#'   [marker_catalog()]), coerced to a list-like with `$release` and
#'   `$file`.
#' @return `character(1)` URL.
#' @noRd
.remote_url <- function(entry) {
  mirror <- getOption("markerdata.mirror")
  if (!is.null(mirror)) {
    base <- mirror
    if (!grepl("/$", base)) base <- paste0(base, "/")
    return(paste0(base, entry$file))
  }
  base <- "https://github.com/fernandes-lab/markerdata/releases/download/"
  paste0(base, entry$release, "/", entry$file)
}

#' Download a file with graceful failure
#'
#' Downloads `url` to `dest`, verifying size and/or md5 when known. Never
#' raises an error: on any failure it emits a `warning()` and returns
#' `FALSE`, per CRAN policy for packages that fetch data over the network.
#'
#' @param url `character(1)` source URL (`file://` URLs are supported via
#'   libcurl).
#' @param dest `character(1)` destination file path.
#' @param md5 `character(1)` expected MD5 checksum, or `NA` if unknown.
#' @param bytes `integer(1)` expected file size in bytes, or `NA` if
#'   unknown.
#' @param quiet `logical(1)` passed through to
#'   [utils::download.file()].
#' @return `logical(1)`: `TRUE` on success, `FALSE` on any failure.
#' @noRd
.download <- function(url, dest, md5 = NA, bytes = NA, quiet = FALSE) {
  old_timeout <- getOption("timeout")
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  tmp <- tempfile()
  status <- tryCatch(
    utils::download.file(url, tmp, mode = "wb", quiet = quiet, method = "libcurl"),
    error = function(e) {
      warning(sprintf(
        "Download failed for %s: %s. Is the machine online? If the problem persists, check that the URL is reachable in a browser.",
        url, conditionMessage(e)
      ), call. = FALSE)
      NA_integer_
    },
    warning = function(w) {
      warning(sprintf(
        "Download failed for %s: %s. Is the machine online? If the problem persists, check that the URL is reachable in a browser.",
        url, conditionMessage(w)
      ), call. = FALSE)
      NA_integer_
    }
  )
  if (is.na(status) || status != 0) {
    if (file.exists(tmp)) unlink(tmp)
    return(FALSE)
  }

  if (!is.na(bytes)) {
    got_bytes <- file.size(tmp)
    if (is.na(got_bytes) || got_bytes != bytes) {
      unlink(tmp)
      warning(sprintf(
        "Downloaded file size (%s) does not match expected size (%s) for %s.",
        got_bytes, bytes, url
      ), call. = FALSE)
      return(FALSE)
    }
  }
  if (!is.na(md5) && nzchar(md5)) {
    got_md5 <- unname(tools::md5sum(tmp))
    if (is.na(got_md5) || !identical(got_md5, md5)) {
      unlink(tmp)
      warning(sprintf(
        "Downloaded file checksum (%s) does not match expected checksum (%s) for %s.",
        got_md5, md5, url
      ), call. = FALSE)
      return(FALSE)
    }
  }

  ok <- suppressWarnings(file.rename(tmp, dest))
  if (!ok) {
    ok2 <- file.copy(tmp, dest, overwrite = TRUE)
    unlink(tmp)
    if (!ok2) {
      warning(sprintf("Could not move downloaded file into place at %s.", dest),
        call. = FALSE
      )
      return(FALSE)
    }
  }
  TRUE
}

#' Get a marker data set
#'
#' Returns a marker data set as a `data.frame`. Every data set is available
#' offline: the maize panel (`SNP55K_maize282_maf04`) is bundled whole, and
#' each of the other data sets ships a sample (2,000 SNPs, every individual)
#' inside the package. Ask for the whole panel with `full = TRUE`; it is
#' downloaded from the package's GitHub releases on first use and cached
#' locally (see [marker_cache_dir()]) so later calls are instant.
#'
#' @section The sample:
#' A sample is a fixed thinning of the full panel, built by
#' `data-raw/samples.R` without random numbers: every individual of the
#' panel, and 2,000 SNPs (1,000 for the four largest panels) evenly spaced over those that are polymorphic among
#' them. It has the same columns and coding as the
#' full panel and is meant for examples, tests and quick simulations; use
#' `full = TRUE` for real studies.
#'
#' @param name `character(1)` data set name; must be one of
#'   `marker_catalog()$name`.
#' @param full `logical(1)`: if `FALSE` (the default), return the sample
#'   shipped with the package; if `TRUE`, return the whole panel,
#'   downloading it when it is not cached yet. Ignored for
#'   `SNP55K_maize282_maf04`, which is bundled whole. For the sample-only
#'   data sets (`Arabidopsis_div`, `Human_div`), markerdata hosts no whole
#'   panel: `full = TRUE` warns, gives the `source_url` of the original
#'   data and returns `NULL`; nothing is downloaded. Use
#'   [marker_download()] to fetch every full panel at once.
#' @param cache `character(1)` cache directory for full panels.
#'   Defaults to [marker_cache_dir()].
#' @param download `logical(1)`: if `FALSE`, never download -- with
#'   `full = TRUE`, only return the panel if it is already cached.
#' @param quiet `logical(1)`: if `FALSE` (the default), print a message
#'   when a sample is returned or a download happens.
#' @return The requested marker `data.frame`, or (invisibly) `NULL` with a
#'   `warning()` if a full panel could not be obtained (for example, no
#'   network connection, or a sample-only data set with `full = TRUE`).
#' @examples
#' head(marker_data("SNP55K_maize282_maf04")[, 1:8])
#'
#' # the sample shipped with the package: no download
#' dim(marker_data("Rice_div"))
#'
#' \donttest{
#' # the whole panel: downloaded on first use, then cached
#' geno <- marker_data("Rice_div", full = TRUE)
#' if (!is.null(geno)) dim(geno)
#' }
#' @export
marker_data <- function(name, full = FALSE, cache = marker_cache_dir(),
                         download = TRUE, quiet = FALSE) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    stop("`name` must be a single string.", call. = FALSE)
  }
  cat <- marker_catalog(cache = cache)
  if (!name %in% cat$name) {
    stop(sprintf(
      "Unknown data set '%s'. Valid names: %s",
      name, paste(cat$name, collapse = ", ")
    ), call. = FALSE)
  }
  entry <- cat[cat$name == name, , drop = FALSE][1, ]

  if (isTRUE(entry$bundled)) {
    e <- new.env()
    utils::data(list = name, package = "markerdata", envir = e)
    return(get(name, envir = e))
  }

  hosted <- .is_hosted(entry)
  # A catalog entry without a shipped sample falls through to the full panel.
  smp <- system.file("extdata", "samples", paste0(name, ".rds"),
    package = "markerdata"
  )
  if (!isTRUE(full) && nzchar(smp)) {
    x <- .unpack(readRDS(smp))
    if (!quiet) {
      message(if (hosted) {
        sprintf(
          "'%s': returning the bundled sample (%d SNPs x %d individuals of %d x %d). Use full = TRUE for the whole panel.",
          name, nrow(x), ncol(x) - 5L, entry$n_snp, entry$n_ind
        )
      } else {
        sprintf(
          "'%s': returning the bundled sample (%d SNPs x %d individuals); the whole data set is at %s.",
          name, nrow(x), ncol(x) - 5L, entry$source_url
        )
      })
    }
    return(x)
  }

  if (!hosted) {
    .warn_not_hosted(entry)
    return(invisible(NULL))
  }

  dest <- file.path(cache, entry$file)
  if (file.exists(dest)) {
    md5 <- entry$md5
    if (!is.na(md5) && nzchar(md5)) {
      got <- unname(tools::md5sum(dest))
      if (!identical(got, md5)) {
        warning(sprintf(
          "Cached file for '%s' is corrupt; re-downloading.", name
        ), call. = FALSE)
        unlink(dest)
      }
    }
  }

  if (!file.exists(dest)) {
    if (!download) {
      warning(sprintf(
        "'%s' is not cached at '%s'. Use marker_download('%s') or marker_data('%s', full = TRUE, download = TRUE) to fetch it.",
        name, cache, name, name
      ), call. = FALSE)
      return(invisible(NULL))
    }
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    url <- .remote_url(entry)
    ok <- .download(url, dest,
      md5 = entry$md5, bytes = entry$bytes,
      quiet = quiet
    )
    if (!ok) {
      return(invisible(NULL))
    }
    if (!quiet) {
      sz <- file.size(dest)
      message(sprintf(
        "Downloaded '%s' (%s, %.1f MB) to cache: %s",
        name, entry$file, sz / 1e6, cache
      ))
    }
  }

  p <- readRDS(dest)
  .unpack(p)
}

#' Download marker data sets into the local cache
#'
#' Downloads one or more full marker panels so subsequent
#' `marker_data(name, full = TRUE)` calls are served from the local cache
#' without network access. Called with no arguments, it downloads every
#' full panel (roughly 19 MB in total; the `bytes` column of
#' [marker_catalog()] gives each size).
#'
#' @param names `character()` vector of data set names to download, or
#'   `NULL` (the default) for every data set whose whole panel markerdata
#'   hosts. The sample-only data sets (`Arabidopsis_div`, `Human_div`) are
#'   skipped silently then; asking for one by name gives a warning with
#'   its `source_url`, and its entry in the result is `FALSE`.
#' @param cache `character(1)` cache directory. Defaults to
#'   [marker_cache_dir()].
#' @param overwrite `logical(1)`: if `TRUE`, re-download files already
#'   present in the cache. Default `FALSE` skips them.
#' @param quiet `logical(1)`: if `FALSE` (the default), print a message
#'   for each download.
#' @return Invisibly, a named `logical` vector (one entry per requested
#'   data set) indicating whether the file is present in the cache after
#'   the call.
#' @examples
#' \donttest{
#' marker_download("Rice_div")
#' }
#' @export
marker_download <- function(names = NULL, cache = marker_cache_dir(),
                             overwrite = FALSE, quiet = FALSE) {
  cat <- marker_catalog(cache = cache)
  remote <- cat[.is_hosted(cat), , drop = FALSE]
  sample_only <- cat[!cat$bundled & !.is_hosted(cat), , drop = FALSE]
  if (is.null(names)) {
    names <- remote$name
  } else {
    unknown <- setdiff(names, c(remote$name, sample_only$name))
    if (length(unknown)) {
      stop(sprintf(
        "Unknown or non-remote data set(s): %s. Valid remote names: %s",
        paste(unknown, collapse = ", "),
        paste(remote$name, collapse = ", ")
      ), call. = FALSE)
    }
  }

  result <- logical(length(names))
  names(result) <- names
  for (nm in names) {
    if (nm %in% sample_only$name) { # no whole panel is hosted: nothing to fetch
      .warn_not_hosted(sample_only[sample_only$name == nm, , drop = FALSE][1, ])
      next
    }
    entry <- remote[remote$name == nm, , drop = FALSE][1, ]
    dest <- file.path(cache, entry$file)
    if (file.exists(dest) && !overwrite) {
      result[[nm]] <- TRUE
      next
    }
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    url <- .remote_url(entry)
    ok <- .download(url, dest, md5 = entry$md5, bytes = entry$bytes, quiet = quiet)
    if (ok && !quiet) {
      sz <- file.size(dest)
      message(sprintf(
        "Downloaded '%s' (%s, %.1f MB) to cache: %s",
        nm, entry$file, sz / 1e6, cache
      ))
    }
    result[[nm]] <- ok
  }
  invisible(result)
}
