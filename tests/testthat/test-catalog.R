# marker_catalog() / marker_references() / .catalog() behavior, using the
# offline fixture catalog, plus a few checks against the real catalog when
# it has been built by the other agent.

test_that("marker_catalog() returns correct types with the fixture catalog", {
  fx <- local_marker_fixture()

  cat <- marker_catalog(cache = fx$cache_dir)

  expect_true(is.data.frame(cat))
  expect_setequal(cat$name, c("SNP55K_maize282_maf04", "Synthetic_div"))
  expect_type(cat$n_snp, "integer")
  expect_type(cat$n_ind, "integer")
  expect_type(cat$bytes, "integer")
  expect_type(cat$het, "double")
  expect_type(cat$bundled, "logical")
  for (col in c(
    "name", "species", "common_name", "panel", "ploidy", "mating",
    "chr_labels", "license", "source", "source_doi", "reference", "doi",
    "file", "release", "md5"
  )) {
    expect_type(cat[[col]], "character")
  }
  expect_identical(rownames(cat), as.character(seq_len(nrow(cat))))
})

test_that("cached/available columns behave with the fixture catalog", {
  fx <- local_marker_fixture()

  cat0 <- marker_catalog(cache = fx$cache_dir)
  syn0 <- cat0[cat0$name == "Synthetic_div", ]
  bundled0 <- cat0[cat0$name == "SNP55K_maize282_maf04", ]

  expect_false(syn0$cached)
  expect_false(syn0$available)
  expect_true(bundled0$bundled)
  expect_false(bundled0$cached) # no remote file to be "cached"; it's bundled instead
  expect_true(bundled0$available) # available = bundled | cached

  # Now place the file in the cache directly and check cached/available flip.
  file.copy(file.path(fx$mirror_dir, fx$file_name), file.path(fx$cache_dir, fx$file_name))
  cat1 <- marker_catalog(cache = fx$cache_dir)
  syn1 <- cat1[cat1$name == "Synthetic_div", ]
  expect_true(syn1$cached)
  expect_true(syn1$available)
})

test_that("unknown name errors in marker_data() list valid names", {
  fx <- local_marker_fixture()
  expect_error(
    marker_data("NoSuchDataset", cache = fx$cache_dir),
    "Unknown data set.*SNP55K_maize282_maf04"
  )
})

test_that("bundled row is not marked 'cached' incorrectly when file is empty", {
  fx <- local_marker_fixture()
  cat <- marker_catalog(cache = fx$cache_dir)
  bundled <- cat[cat$bundled, , drop = FALSE]
  expect_true(all(bundled$file == "" | is.na(bundled$file)))
})

# --- Checks against the real catalog (skipped until the other agent's
# catalog.csv / sysdata.rda exist) ---------------------------------------

real_catalog_path <- function() {
  system.file("extdata", "catalog.csv", package = "markerdata")
}

test_that("real catalog: names are unique and SNP55K bundled row matches data()", {
  path <- real_catalog_path()
  skip_if(!nzchar(path) || !file.exists(path), "real catalog.csv not present yet")

  cat <- marker_catalog()
  expect_false(anyDuplicated(cat$name) > 0)

  row <- cat[cat$name == "SNP55K_maize282_maf04", ]
  expect_equal(nrow(row), 1L)
  expect_true(row$bundled)

  geno <- marker_data("SNP55K_maize282_maf04")
  expect_equal(dim(geno), c(10650L, 285L))

  e <- new.env()
  utils::data(list = "SNP55K_maize282_maf04", package = "markerdata", envir = e)
  expect_identical(geno, get("SNP55K_maize282_maf04", envir = e))
})

test_that("real references database: every catalog name has an entry", {
  path <- real_catalog_path()
  skip_if(!nzchar(path) || !file.exists(path), "real catalog.csv not present yet")

  has_refs <- tryCatch(
    {
      exists(".marker_refs", envir = asNamespace("markerdata"), inherits = FALSE)
    },
    error = function(e) FALSE
  )
  skip_if_not(has_refs, "sysdata .marker_refs not present yet")

  cat <- marker_catalog()
  refs <- get(".marker_refs", envir = asNamespace("markerdata"))
  expect_true(all(cat$name %in% names(refs)))

  r <- marker_references()
  expect_s3_class(r, "bibentry")

  one <- marker_references("SNP55K_maize282_maf04")
  txt <- format(one, style = "text")
  expect_true(length(txt) > 0)
  expect_true(any(nzchar(txt)))

  expect_error(print(one), NA)
})
