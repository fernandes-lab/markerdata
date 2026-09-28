#' Local cache directory for downloaded marker data sets
#'
#' Resolves the directory where remote (non-bundled) marker data sets are
#' cached after download. The directory is not created by this function;
#' see [marker_data()] and [marker_download()], which create it as needed.
#'
#' Resolution order:
#' 1. `getOption("markerdata.cache")`, if set.
#' 2. The environment variable `MARKERDATA_CACHE`, if set to a non-empty
#'    string.
#' 3. `tools::R_user_dir("markerdata", "cache")`.
#'
#' @return `character(1)`, the path to the cache directory.
#' @examples
#' marker_cache_dir()
#' @export
marker_cache_dir <- function() {
  opt <- getOption("markerdata.cache")
  if (!is.null(opt)) {
    return(opt)
  }
  env <- Sys.getenv("MARKERDATA_CACHE", unset = "")
  if (nzchar(env)) {
    return(env)
  }
  tools::R_user_dir("markerdata", "cache")
}

#' Remove cached marker data files
#'
#' Deletes the cached file(s) for one or more remote marker data sets from
#' the local cache directory.
#'
#' @param names `character()` vector of data set names to clear, or `NULL`
#'   (the default) to clear every remote data set present in the catalog.
#' @param cache `character(1)` cache directory to clear from. Defaults to
#'   [marker_cache_dir()].
#' @return Invisibly, a `character()` vector of the file paths that were
#'   removed.
#' @examples
#' \donttest{
#' marker_cache_clear("Rice_div")
#' }
#' @export
marker_cache_clear <- function(names = NULL, cache = marker_cache_dir()) {
  cat <- .catalog()
  remote <- cat[!cat$bundled, , drop = FALSE]
  if (!is.null(names)) {
    unknown <- setdiff(names, remote$name)
    if (length(unknown)) {
      stop(sprintf(
        "Unknown or non-remote data set(s): %s. Valid remote names: %s",
        paste(unknown, collapse = ", "),
        paste(remote$name, collapse = ", ")
      ), call. = FALSE)
    }
    remote <- remote[remote$name %in% names, , drop = FALSE]
  }
  files <- file.path(cache, remote$file)
  present <- files[nzchar(remote$file) & file.exists(files)]
  if (length(present)) {
    unlink(present)
  }
  invisible(present)
}
