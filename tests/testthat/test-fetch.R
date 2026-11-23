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
  geno <- suppressWarnings(marker_data(row$name, full = TRUE, cache = tmp_cache))
  if (is.null(geno)) {
    skip("release asset not reachable")
  }
  expect_true(is.data.frame(geno))
})

test_that("every remote panel ships a sample that loads offline", {
  cat <- .catalog()
  remote <- cat$name[!cat$bundled]
  expect_length(remote, 16L) # 17 data sets in all, one of them bundled whole
  withr::local_options(markerdata.mirror = "file:///no/such/mirror/dir")
  cache <- withr::local_tempdir()
  for (nm in remote) {
    expect_message(s <- marker_data(nm, cache = cache), "bundled sample")
    expect_true(.check_marker(s, nm))
    # 2,000 SNPs (1,000 for soybean, cattle, wheat and sugar beet) and every individual of the panel
    n_ind <- cat$n_ind[cat$name == nm]
    expect_identical(dim(s), c(if (nm %in% c("SoySNP50K_div", "Cattle_div", "Wheat_div", "SugarBeet_div")) 1000L else 2000L, 5L + n_ind), info = nm)
    # the sample is polymorphic at every SNP
    g <- as.matrix(s[, -(1:5)])
    expect_true(all(rowSums(g != g[, 1L]) > 0L))
  }
  expect_length(list.files(cache), 0L) # nothing was downloaded
})

test_that("every hosted row names its release tag and asset; the data-v2 tag is the newer seven", {
  cat <- .catalog()
  hosted <- cat[.is_hosted(cat), , drop = FALSE]
  expect_length(hosted$name, 14L)
  expect_false(anyNA(hosted$bytes))
  expect_true(all(grepl("^[0-9a-f]{32}$", hosted$md5)))
  expect_identical(hosted$file, paste0(hosted$name, ".rds"))
  expect_setequal(hosted$release, c("data-v1", "data-v2"))
  expect_identical(sum(hosted$release == "data-v1"), 7L)
  expect_setequal(hosted$name[hosted$release == "data-v2"], c(
    "SugarBeet_div", "Tomato_div", "WheatHex_div", "Salmon_div", "Sheep_div",
    "Potato_div", "Sugarcane_div"
  ))
  # the URL of an asset contains its own release tag
  v2 <- hosted[hosted$release == "data-v2", , drop = FALSE][1, ]
  expect_match(.remote_url(v2), paste0("/", v2$release, "/", v2$file), fixed = TRUE)
})

test_that("full = TRUE on every hosted panel fails gracefully when its asset is unreachable", {
  # e.g. the data-v2 release is not published yet: warning, NULL, no error
  withr::local_options(markerdata.mirror = "file:///no/such/mirror/dir")
  cache <- withr::local_tempdir()
  cat <- .catalog()
  for (nm in cat$name[.is_hosted(cat)]) {
    expect_warning(geno <- marker_data(nm, full = TRUE, cache = cache, quiet = TRUE),
      "Download failed", info = nm)
    expect_null(geno)
  }
  expect_length(list.files(cache), 0L)
})

test_that("full = TRUE skips the sample and goes to the cache/download path", {
  withr::local_options(markerdata.mirror = "file:///no/such/mirror/dir")
  cache <- withr::local_tempdir()
  expect_warning(
    geno <- marker_data("Rice_div", full = TRUE, cache = cache, download = FALSE),
    "not cached"
  )
  expect_null(geno)
  expect_silent(marker_data("Rice_div", cache = cache, quiet = TRUE))
})

test_that("full = TRUE on a sample-only data set warns with the source URL, returns NULL, never downloads", {
  # a mirror that cannot be reached: any download attempt would add a second warning
  withr::local_options(markerdata.mirror = "file:///no/such/mirror/dir")
  cache <- withr::local_tempdir()
  cat <- .catalog()
  only <- cat[!cat$bundled & (is.na(cat$file) | !nzchar(cat$file)), ]
  expect_setequal(only$name, c("Arabidopsis_div", "Human_div"))
  expect_true(all(is.na(only$bytes) & is.na(only$md5) & is.na(only$release)))
  for (nm in only$name) {
    url <- only$source_url[only$name == nm]
    expect_match(url, "^(https|ftp)://")
    w <- testthat::capture_warnings(
      geno <- marker_data(nm, full = TRUE, cache = cache)
    )
    expect_length(w, 1L)
    expect_match(w, "not hosted by markerdata", fixed = TRUE)
    expect_match(w, url, fixed = TRUE)
    expect_null(geno)
    # also with download = FALSE
    expect_warning(marker_data(nm, full = TRUE, cache = cache, download = FALSE), url, fixed = TRUE)
  }
  expect_length(list.files(cache), 0L)
  # the sample itself is fine and points at the original data
  expect_message(s <- marker_data("Human_div", cache = cache), "ftp.1000genomes")
  expect_identical(dim(s), c(2000L, 5L + 2504L))
})

test_that("marker_download() and marker_cache_clear() skip sample-only rows", {
  fx <- local_marker_fixture()
  cat_df <- utils::read.csv(fx$catalog_path, stringsAsFactors = FALSE, na.strings = "")
  cat_df$source_url <- NA_character_
  extra <- cat_df[cat_df$name == "SNP55K_maize282_maf04", ]
  extra$name <- "Sample_only_div"
  extra$bundled <- FALSE
  extra$source_url <- "https://example.org/whole-data/"
  utils::write.csv(rbind(cat_df, extra), fx$catalog_path, row.names = FALSE, na = "")

  # names = NULL: downloads the hosted panel only, no error, no warning
  expect_no_warning(res <- marker_download(cache = fx$cache_dir, quiet = TRUE))
  expect_named(res, "Synthetic_div")
  expect_true(res[["Synthetic_div"]])
  expect_no_error(removed <- marker_cache_clear(cache = fx$cache_dir))
  expect_length(removed, 1L)

  # explicit name: warning that gives the URL, FALSE in the result, nothing fetched
  expect_warning(
    res <- marker_download(c("Sample_only_div", "Synthetic_div"), cache = fx$cache_dir, quiet = TRUE),
    "https://example.org/whole-data/",
    fixed = TRUE
  )
  expect_identical(res, c(Sample_only_div = FALSE, Synthetic_div = TRUE))
  expect_identical(list.files(fx$cache_dir), fx$file_name)
  expect_warning(marker_cache_clear("Sample_only_div", cache = fx$cache_dir), "nothing to clear")
  # truly unknown names still error
  expect_error(marker_download("NoSuchDataset", cache = fx$cache_dir), "Unknown")
})
