# .pack() / .unpack() round-trip and .check_marker() contract checks.

make_marker <- function(n_snp, n_ind, chr_type = c("character", "integer"),
                         geno_type = c("integer", "double"),
                         row_names = c("automatic", "character")) {
  chr_type <- match.arg(chr_type)
  geno_type <- match.arg(geno_type)
  row_names <- match.arg(row_names)

  n_chr <- 3L
  chr_id <- rep(seq_len(n_chr), length.out = n_snp)
  chr_id <- sort(chr_id)
  pos <- seq_len(n_snp) * 100L
  cm <- as.numeric(seq_len(n_snp))
  chr <- if (chr_type == "character") as.character(chr_id) else as.integer(chr_id)

  geno <- matrix(sample(c(-1L, 0L, 1L), n_snp * n_ind, replace = TRUE),
    nrow = n_snp, ncol = n_ind
  )
  if (geno_type == "double") storage.mode(geno) <- "double"
  colnames(geno) <- sprintf("I%d", seq_len(n_ind))

  df <- data.frame(
    snp = sprintf("s%d", seq_len(n_snp)),
    allele = rep("A/T", n_snp),
    chr = chr,
    pos = pos,
    cm = cm,
    as.data.frame(geno, check.names = FALSE),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (row_names == "character") {
    rownames(df) <- sprintf("row_%d", seq_len(n_snp))
  }
  df
}

test_that("round trip is identical: integer genotypes + character chr", {
  x <- make_marker(20, 6, "character", "integer", "automatic")
  p <- .pack(x, "t1")
  expect_true(.row_names_info(x) < 0)
  expect_identical(.unpack(p), x)
})

test_that("round trip is identical: double genotypes + integer chr", {
  x <- make_marker(20, 6, "integer", "double", "automatic")
  p <- .pack(x, "t2")
  expect_identical(.unpack(p), x)
})

test_that("round trip is identical: character row names", {
  x <- make_marker(15, 5, "character", "integer", "character")
  p <- .pack(x, "t3")
  expect_false(.row_names_info(x) < 0)
  expect_identical(.unpack(p), x)
})

test_that("round trip is identical: n_snp * n_ind not divisible by 4", {
  x <- make_marker(7, 3, "character", "integer", "automatic") # 21 calls
  expect_true((nrow(x) * (ncol(x) - 5L)) %% 4L != 0L)
  p <- .pack(x, "t4")
  expect_identical(.unpack(p), x)
})

test_that("layout selection returns 'snp' or 'ind'", {
  x <- make_marker(30, 8, "character", "integer", "automatic")
  p <- .pack(x, "t5")
  expect_true(p$layout %in% c("snp", "ind"))
})

test_that(".unpack() rejects format != 1", {
  x <- make_marker(10, 4, "character", "integer", "automatic")
  p <- .pack(x, "t6")
  p_bad <- p
  p_bad$format <- 2L
  expect_error(.unpack(p_bad), "format")

  p_null <- p
  p_null$format <- NULL
  expect_error(.unpack(p_null), "format")
})

test_that(".check_marker() accepts a valid marker data.frame", {
  x <- make_marker(10, 4, "character", "integer", "automatic")
  expect_true(.check_marker(x))
})

test_that(".check_marker() rejects wrong column names", {
  x <- make_marker(10, 4)
  names(x)[1] <- "SNP"
  expect_error(.check_marker(x), "snp")
})

test_that(".check_marker() rejects NA in cm", {
  x <- make_marker(10, 4)
  x$cm[3] <- NA
  expect_error(.check_marker(x), "cm")
})

test_that(".check_marker() rejects genotype value 2", {
  x <- make_marker(10, 4)
  x[[6]][1] <- 2
  expect_error(.check_marker(x), "-1/0/1")
})

test_that(".check_marker() rejects duplicated individual names", {
  x <- make_marker(10, 4)
  names(x)[7] <- names(x)[6]
  expect_error(.check_marker(x), "unique")
})

test_that("built marker data sets round-trip through pack/unpack", {
  built_dir <- testthat::test_path("..", "..", "data-raw", "built")
  skip_if_not(dir.exists(built_dir), "data-raw/built not present (other agent's output)")

  files <- list.files(built_dir, pattern = "\\.rda$", full.names = TRUE)
  skip_if(length(files) == 0, "no built .rda files found")

  for (f in files) {
    e <- new.env()
    load(f, envir = e)
    nm <- ls(e)
    expect_length(nm, 1L)
    x <- get(nm, envir = e)
    expect_true(.check_marker(x, nm))
    p <- .pack(x, nm)
    expect_identical(.unpack(p), x, info = nm)
  }
})
