# data-raw/thin.R -- the one thinning rule used for every package sample.
# Sourced by data-raw/samples.R and data-raw/s1database_panels.R.
#   1. keep all individuals (N_IND = Inf; a finite N_IND would keep that many,
#      evenly spaced over the columns);
#   2. drop SNPs that are monomorphic among those individuals;
#   3. keep N_SNP of the remaining SNPs, evenly spaced over the rows, so every
#      chromosome keeps its share and `cm` stays sorted.
N_SNP <- 2000L
N_IND <- Inf
N_SNP_OVERRIDE <- c(SoySNP50K_div = 1000L, Cattle_div = 1000L, Wheat_div = 1000L, SugarBeet_div = 1000L) # keeps the package tarball under CRAN's 5 MB

even <- function(n, k) unique(round(seq(1, n, length.out = min(k, n))))

thin_sample <- function(x, name = NA_character_) {
  n_snp <- if (name %in% names(N_SNP_OVERRIDE)) N_SNP_OVERRIDE[[name]] else N_SNP
  ind <- 5L + even(ncol(x) - 5L, N_IND)
  g <- as.matrix(x[, ind])
  poly <- which(rowSums(g != g[, 1L]) > 0L)
  rows <- poly[even(length(poly), n_snp)]
  s <- x[rows, c(1:5, ind)]
  rownames(s) <- NULL
  s
}
