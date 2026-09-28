# marker_cache_dir() resolution order: option > env var > R_user_dir().

test_that("option takes precedence over env var and R_user_dir()", {
  withr::local_envvar(MARKERDATA_CACHE = "/env/cache/dir")
  withr::local_options(markerdata.cache = "/option/cache/dir")
  expect_identical(marker_cache_dir(), "/option/cache/dir")
})

test_that("env var takes precedence over R_user_dir() when option unset", {
  withr::local_options(markerdata.cache = NULL)
  withr::local_envvar(MARKERDATA_CACHE = "/env/cache/dir")
  expect_identical(marker_cache_dir(), "/env/cache/dir")
})

test_that("an empty env var falls through to R_user_dir()", {
  withr::local_options(markerdata.cache = NULL)
  withr::local_envvar(MARKERDATA_CACHE = "")
  expect_identical(marker_cache_dir(), tools::R_user_dir("markerdata", "cache"))
})

test_that("falls back to R_user_dir() when neither option nor env var is set", {
  withr::local_options(markerdata.cache = NULL)
  withr::local_envvar(MARKERDATA_CACHE = NA)
  expect_identical(marker_cache_dir(), tools::R_user_dir("markerdata", "cache"))
})
