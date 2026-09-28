#' markerdata: Marker Data Sets for Breeding Simulation
#'
#' Curated SNP marker panels for breeding-program simulation. One small
#' maize panel is bundled with the package; the larger panels (cattle,
#' eucalyptus, poplar, rice, sorghum, soybean, wheat) are downloaded on
#' demand from the package's GitHub release and cached locally, so they do
#' not bloat the installed package. See [marker_data()] to load a data
#' set (downloading and caching it if needed), [marker_catalog()] to list
#' every available data set and its availability, and
#' [marker_references()] for the literature citations behind each one.
#'
#' @keywords internal
#' @importFrom utils download.file read.csv data
#' @importFrom tools md5sum R_user_dir
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
