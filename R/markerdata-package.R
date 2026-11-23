#' markerdata: Marker Data Sets for Breeding Simulation
#'
#' Curated SNP marker panels for breeding-program simulation, across
#' crops, trees, livestock and fish. One small maize panel is bundled whole;
#' every other panel ships a 2,000-SNP sample (every individual), and
#' the whole panel is downloaded on demand from the package's GitHub
#' releases and cached locally, so it does not bloat the installed package.
#' Arabidopsis and human (1001 Genomes, 1000 Genomes Project) are sample
#' only, with the original download in the catalog's `source_url`. See [marker_data()] to load a data set (the
#' sample, or the whole panel with `full = TRUE`), [marker_catalog()] to list
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
