# Test fixtures shared across test files: a small synthetic marker
# data.frame, a packed .rds for it, and a fixture catalog.csv pointing at a
# flat file:// "mirror" directory plus a separate empty cache directory.
# Everything is offline and self-contained under withr-managed temp dirs.

#' Build a small synthetic marker data.frame
#'
#' 40 SNPs across 3 chromosomes (character `chr`), 12 individuals, integer
#' genotypes in `{-1, 0, 1}` with some heterozygotes, automatic row names.
#'
#' @return A marker `data.frame` satisfying the package's marker contract.
make_synthetic_marker <- function() {
  set.seed(42)
  n_snp <- 40L
  n_ind <- 12L
  chr <- rep(c("1", "2", "3"), c(14, 13, 13))
  pos <- unlist(lapply(c(14, 13, 13), function(n) sort(sample(1:1e6, n))))
  cm <- unlist(lapply(c(14, 13, 13), function(n) sort(runif(n, 0, 100))))
  snp <- sprintf("S%s_%d", chr, pos)
  allele <- sample(c("A/G", "C/T", "A/T", "C/G"), n_snp, replace = TRUE)

  geno <- matrix(
    sample(c(-1L, 0L, 1L), n_snp * n_ind,
      replace = TRUE, prob = c(0.4, 0.3, 0.3)
    ),
    nrow = n_snp, ncol = n_ind
  )
  colnames(geno) <- sprintf("Ind%02d", seq_len(n_ind))

  data.frame(
    snp = snp, allele = allele, chr = chr, pos = as.integer(pos), cm = cm,
    as.data.frame(geno, check.names = FALSE),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

#' Set up an offline fixture catalog + mirror + cache for a test
#'
#' Builds the synthetic marker data.frame, packs and writes it to
#' `<mirror_dir>/Synthetic_div.rds`, writes a fixture `catalog.csv`
#' listing the bundled `SNP55K_maize282_maf04` and the remote
#' `Synthetic_div`, and sets `markerdata.catalog`, `markerdata.mirror`
#' (`file://<mirror_dir>`) and `markerdata.cache` options scoped to the
#' calling test via [withr::local_options()].
#'
#' @param env Environment the option/scoping should be tied to; defaults
#'   to the caller's frame so options are restored when that test exits.
#' @return A list with elements `synthetic` (the synthetic data.frame),
#'   `catalog_path`, `mirror_dir`, `cache_dir`, `file_name`, `md5`,
#'   `bytes`.
local_marker_fixture <- function(env = parent.frame()) {
  mirror_dir <- withr::local_tempdir(.local_envir = env)
  cache_dir <- withr::local_tempdir(.local_envir = env)

  synthetic <- make_synthetic_marker()
  packed <- .pack(synthetic, "Synthetic_div")
  file_name <- "Synthetic_div.rds"
  dest <- file.path(mirror_dir, file_name)
  saveRDS(packed, dest, compress = "xz")

  bytes <- file.size(dest)
  md5 <- unname(tools::md5sum(dest))

  catalog_path <- file.path(mirror_dir, "catalog.csv")
  cat_df <- data.frame(
    name = c("SNP55K_maize282_maf04", "Synthetic_div"),
    species = c("Zea mays", "Test species"),
    common_name = c("maize", "synthetic"),
    panel = c("maize282", "synthetic-panel"),
    ploidy = c("2", "2"),
    mating = c("inbred", "mixed"),
    n_snp = c(10650L, nrow(synthetic)),
    n_ind = c(280L, ncol(synthetic) - 5L),
    het = c(0.01, 0.3),
    chr_labels = c("1-10", "1-3"),
    license = c("CC0", "CC0"),
    source = c("Panzea", "unit-test fixture"),
    source_doi = c("", ""),
    reference = c("", ""),
    doi = c("", ""),
    bundled = c(TRUE, FALSE),
    file = c("", file_name),
    release = c("", "data-v1"),
    bytes = c("", as.character(bytes)),
    md5 = c("", md5),
    stringsAsFactors = FALSE
  )
  utils::write.csv(cat_df, catalog_path, row.names = FALSE, na = "")

  withr::local_options(
    list(
      markerdata.catalog = catalog_path,
      markerdata.mirror = paste0("file://", mirror_dir),
      markerdata.cache = cache_dir
    ),
    .local_envir = env
  )

  list(
    synthetic = synthetic,
    catalog_path = catalog_path,
    mirror_dir = mirror_dir,
    cache_dir = cache_dir,
    file_name = file_name,
    md5 = md5,
    bytes = bytes
  )
}
