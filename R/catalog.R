#' Read the raw catalog CSV with coerced column types
#'
#' @return A `data.frame` with one row per data set, column types coerced as
#'   documented in [marker_catalog()].
#' @noRd
.catalog <- function() {
  path <- getOption("markerdata.catalog")
  if (is.null(path)) {
    path <- system.file("extdata", "catalog.csv", package = "markerdata")
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop(
      "Could not find the markerdata catalog file. ",
      "Expected it at inst/extdata/catalog.csv (or the path set in ",
      "options(markerdata.catalog = ...)).",
      call. = FALSE
    )
  }
  cat <- utils::read.csv(path,
    stringsAsFactors = FALSE, na.strings = "",
    check.names = FALSE
  )
  int_cols <- c("n_snp", "n_ind", "bytes")
  for (col in int_cols) {
    if (col %in% names(cat)) cat[[col]] <- as.integer(cat[[col]])
  }
  if ("het" %in% names(cat)) cat$het <- as.numeric(cat$het)
  if ("bundled" %in% names(cat)) cat$bundled <- as.logical(cat$bundled)
  char_cols <- setdiff(names(cat), c(int_cols, "het", "bundled"))
  for (col in char_cols) {
    if (!is.character(cat[[col]])) cat[[col]] <- as.character(cat[[col]])
  }
  cat
}

#' Which catalog rows have a whole panel hosted by markerdata
#'
#' Hosted = not bundled and with a release `file`. Rows with neither are
#' sample-only (e.g. `Arabidopsis_div`, `Human_div`): only the sample ships,
#' and `source_url` points to the original full download.
#'
#' @param cat A catalog `data.frame`.
#' @return `logical` vector, one per row.
#' @noRd
.is_hosted <- function(cat) {
  !cat$bundled & !is.na(cat$file) & nzchar(cat$file)
}

#' Warn that a data set has no hosted whole panel, giving its source URL
#' @noRd
.warn_not_hosted <- function(entry, what = "The whole data set") {
  url <- if ("source_url" %in% names(entry)) entry$source_url else NA_character_
  where <- if (is.na(url) || !nzchar(url)) {
    "see the original source (no URL is recorded in the catalog)"
  } else {
    paste0("download it from the original source: ", url)
  }
  warning(sprintf(
    "'%s' is a sample-only data set. %s is not hosted by markerdata; %s",
    entry$name, what, where
  ), call. = FALSE)
}

#' List available marker data sets
#'
#' Returns the package's catalog of marker data sets: one small panel
#' bundled whole plus a set of larger panels, each shipped as a sample and
#' downloaded whole on demand (see [marker_data()]). Every row describes
#' the full panel, except the rows of the sample-only data sets
#' (`Arabidopsis_div`, `Human_div`), where `n_snp`, `n_ind` and `het`
#' describe the sample, because markerdata does not host the whole data
#' set: its `file`, `release`, `bytes` and `md5` are empty and `source_url`
#' says where to download the original. The catalog reports which full
#' panels are already available on this machine.
#'
#' @param cache `character(1)` cache directory to check for already
#'   downloaded remote data sets. Defaults to [marker_cache_dir()].
#' @return A `data.frame`, one row per data set, with the catalog columns
#'   (`name`, `species`, `common_name`, `panel`, `ploidy`, `mating`,
#'   `n_snp`, `n_ind`, `het`, `chr_labels`, `license`, `source`,
#'   `source_doi`, `source_url`, `reference`, `doi`, `bundled`, `file`, `release`,
#'   `bytes`, `md5`) plus two logical columns: `cached` (the full panel's
#'   file is already present in `cache`) and `available`
#'   (`bundled | cached`, i.e. the full panel can be loaded without a
#'   download; the sample always can). `source_url` is the page or
#'   directory where the original full data can be downloaded (empty when
#'   none is recorded).
#' @examples
#' cat <- marker_catalog()
#' cat[, c("name", "bundled", "available")]
#' @export
marker_catalog <- function(cache = marker_cache_dir()) {
  cat <- .catalog()
  cached <- rep(FALSE, nrow(cat))
  has_file <- !is.na(cat$file) & nzchar(cat$file)
  if (any(has_file)) {
    cached[has_file] <- file.exists(file.path(cache, cat$file[has_file]))
  }
  cat$cached <- cached
  cat$available <- cat$bundled | cat$cached
  rownames(cat) <- NULL
  cat
}

#' Bibliographic references for marker data sets
#'
#' Returns the citation(s) for one or more marker data sets, drawn from the
#' package's internal references database.
#'
#' @param name `character()` vector of data set name(s), or `NULL` (the
#'   default) for every data set in the catalog.
#' @return A single `bibentry` object (the requested entries concatenated
#'   with [c()] when more than one is requested).
#' @examples
#' refs <- marker_references("SNP55K_maize282_maf04")
#' print(refs)
#' @export
marker_references <- function(name = NULL) {
  if (!exists(".marker_refs", envir = asNamespace("markerdata"), inherits = FALSE)) {
    stop("References database missing.", call. = FALSE)
  }
  refs <- get(".marker_refs", envir = asNamespace("markerdata"))
  valid <- names(refs)
  if (is.null(name)) {
    name <- valid
  }
  unknown <- setdiff(name, valid)
  if (length(unknown)) {
    stop(sprintf(
      "Unknown data set(s): %s. Valid names: %s",
      paste(unknown, collapse = ", "),
      paste(valid, collapse = ", ")
    ), call. = FALSE)
  }
  do.call(c, unname(refs[name]))
}
