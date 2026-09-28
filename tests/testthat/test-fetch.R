# marker_data() / marker_download() / marker_cache_clear() over an offline
# fixture (a file:// mirror), covering the download/caching/corruption
# state machine without any real network access.

test_that("first call downloads (with a message) and matches the synthetic data", {
  fx <- local_marker_fixture()

  expect_message(
    geno <- marker_data("Synthetic_div", cache = fx$cache_dir),
    "Downloaded"
  )
  expect_identical(geno, fx$synthetic)
  expect_true(file.exists(file.path(fx$cache_dir, fx$file_name)))
})

test_that("a cached file is reused without re-downloading", {
  fx <- local_marker_fixture()
  geno1 <- marker_data("Synthetic_div", cache = fx$cache_dir, quiet = TRUE)
  expect_identical(geno1, fx$synthetic)

  # Point the mirror at a nonexistent location: if marker_data() tried to
  # download again, it would fail/warn. A silent, successful call proves
  # the cached file was used instead.
  withr::local_options(markerdata.mirror = "file:///no/such/mirror/dir")
  expect_no_condition({
    geno2 <- marker_data("Synthetic_div", cache = fx$cache_dir)
  })
  expect_identical(geno2, fx$synthetic)
})

test_that("a corrupted cached file is detected, re-downloaded and fixed", {
  fx <- local_marker_fixture()
  marker_data("Synthetic_div", cache = fx$cache_dir, quiet = TRUE)

  dest <- file.path(fx$cache_dir, fx$file_name)
  writeBin(as.raw(c(1, 2, 3, 4, 5)), dest)

  expect_warning(
    geno <- marker_data("Synthetic_div", cache = fx$cache_dir, quiet = TRUE),
    "corrupt"
  )
  expect_identical(geno, fx$synthetic)
})

test_that("md5 mismatch in the catalog: warning, NULL, and no bad file left cached", {
  fx <- local_marker_fixture()

  cat_df <- utils::read.csv(fx$catalog_path, stringsAsFactors = FALSE, na.strings = "")
  cat_df$md5[cat_df$name == "Synthetic_div"] <- "0000000000000000000000000000bad"
  utils::write.csv(cat_df, fx$catalog_path, row.names = FALSE, na = "")

  expect_warning(
    geno <- marker_data("Synthetic_div", cache = fx$cache_dir),
    "checksum"
  )
  expect_null(geno)
  expect_false(file.exists(file.path(fx$cache_dir, fx$file_name)))
})

test_that("download = FALSE on a cache miss warns and returns NULL", {
  fx <- local_marker_fixture()
  expect_warning(
    geno <- marker_data("Synthetic_div", cache = fx$cache_dir, download = FALSE),
    "not cached"
  )
  expect_null(geno)
  expect_false(file.exists(file.path(fx$cache_dir, fx$file_name)))
})

test_that("an unreachable mirror warns and returns NULL (never errors)", {
  fx <- local_marker_fixture()
  withr::local_options(markerdata.mirror = "file:///no/such/mirror/dir")

  expect_warning(
    geno <- marker_data("Synthetic_div", cache = fx$cache_dir)
  )
  expect_null(geno)
  expect_false(file.exists(file.path(fx$cache_dir, fx$file_name)))
})

test_that("marker_download() returns a named logical and populates the cache", {
  fx <- local_marker_fixture()
  res <- marker_download("Synthetic_div", cache = fx$cache_dir, quiet = TRUE)
  expect_type(res, "logical")
  expect_named(res, "Synthetic_div")
  expect_true(res[["Synthetic_div"]])
  expect_true(file.exists(file.path(fx$cache_dir, fx$file_name)))
})

test_that("marker_download() skips already-cached files unless overwrite = TRUE", {
  fx <- local_marker_fixture()
  marker_download("Synthetic_div", cache = fx$cache_dir, quiet = TRUE)
  dest <- file.path(fx$cache_dir, fx$file_name)
  mtime1 <- file.mtime(dest)
  Sys.sleep(1.1)
  res <- marker_download("Synthetic_div", cache = fx$cache_dir, quiet = TRUE)
  expect_true(res[["Synthetic_div"]])
  expect_equal(file.mtime(dest), mtime1)
})

test_that("marker_cache_clear() removes the cached file", {
  fx <- local_marker_fixture()
  marker_download("Synthetic_div", cache = fx$cache_dir, quiet = TRUE)
  dest <- file.path(fx$cache_dir, fx$file_name)
  expect_true(file.exists(dest))

  removed <- marker_cache_clear("Synthetic_div", cache = fx$cache_dir)
  expect_true(dest %in% removed)
  expect_false(file.exists(dest))
})

test_that("real network: fetch the smallest remote data set", {
  skip_on_cran()
  skip_if_offline()
  path <- system.file("extdata", "catalog.csv", package = "markerdata")
  skip_if(!nzchar(path) || !file.exists(path), "real catalog.csv not present yet")

  cat <- marker_catalog()
  remote <- cat[!cat$bundled & !is.na(cat$bytes), , drop = FALSE]
  skip_if(nrow(remote) == 0, "no remote data set with known size yet")
  remote <- remote[order(remote$bytes), ]
  row <- remote[1, ]
  skip_if(is.na(row$md5) || !nzchar(row$md5), "no md5 recorded for the smallest remote data set yet")

  tmp_cache <- withr::local_tempdir()
  geno <- suppressWarnings(marker_data(row$name, cache = tmp_cache))
  if (is.null(geno)) {
    skip("release asset not reachable")
  }
  expect_true(is.data.frame(geno))
})
