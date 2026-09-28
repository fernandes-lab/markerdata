# Internal packed-file format (format 1) and marker-data.frame contract
# checks. Not exported; documented with @noRd so no Rd pages are generated.

#' Check that a data frame satisfies the marker-data contract
#'
#' @param x A `data.frame` to validate.
#' @param name `character(1)` used to identify `x` in error messages.
#' @return `invisible(TRUE)` if `x` satisfies the contract; otherwise this
#'   function stops with an informative message.
#' @noRd
.check_marker <- function(x, name = "geno") {
  if (!is.data.frame(x)) {
    stop(sprintf("`%s` must be a data.frame.", name), call. = FALSE)
  }
  if (ncol(x) < 6L) {
    stop(sprintf(
      "`%s` must have at least 6 columns (5 map columns + 1 individual).",
      name
    ), call. = FALSE)
  }
  req <- c("snp", "allele", "chr", "pos", "cm")
  got <- colnames(x)[1:5]
  if (!identical(got, req)) {
    stop(sprintf(
      "`%s`: the first five columns must be named %s, got %s.",
      name,
      paste(sprintf("'%s'", req), collapse = ", "),
      paste(sprintf("'%s'", got), collapse = ", ")
    ), call. = FALSE)
  }
  if (!is.character(x$snp)) {
    stop(sprintf("`%s$snp` must be character.", name), call. = FALSE)
  }
  if (!is.character(x$allele)) {
    stop(sprintf("`%s$allele` must be character.", name), call. = FALSE)
  }
  if (!(is.character(x$chr) || is.integer(x$chr))) {
    stop(sprintf("`%s$chr` must be character or integer.", name), call. = FALSE)
  }
  if (anyNA(x$chr)) {
    stop(sprintf("`%s$chr` must not contain NA.", name), call. = FALSE)
  }
  if (!is.numeric(x$pos)) {
    stop(sprintf("`%s$pos` must be integer/numeric.", name), call. = FALSE)
  }
  if (anyNA(x$pos)) {
    stop(sprintf("`%s$pos` must not contain NA.", name), call. = FALSE)
  }
  if (!is.numeric(x$cm)) {
    stop(sprintf("`%s$cm` must be numeric.", name), call. = FALSE)
  }
  if (anyNA(x$cm)) {
    stop(sprintf("`%s$cm` must not contain NA.", name), call. = FALSE)
  }
  chr_chr <- as.character(x$chr)
  ok_order <- unsplit_ok <- TRUE
  for (lev in unique(chr_chr)) {
    v <- x$cm[chr_chr == lev]
    if (length(v) > 1L && any(diff(v) < 0)) {
      ok_order <- FALSE
      break
    }
  }
  if (!ok_order) {
    stop(sprintf(
      "`%s$cm` must be non-decreasing within each chromosome.", name
    ), call. = FALSE)
  }

  ind_cols <- seq.int(6L, ncol(x))
  ind_names <- colnames(x)[ind_cols]
  if (any(!nzchar(ind_names))) {
    stop(sprintf("`%s`: individual column names must be non-empty.", name),
      call. = FALSE
    )
  }
  if (anyDuplicated(ind_names)) {
    stop(sprintf("`%s`: individual column names must be unique.", name),
      call. = FALSE
    )
  }
  for (j in ind_cols) {
    col <- x[[j]]
    if (!is.numeric(col)) {
      stop(sprintf(
        "`%s`: individual column '%s' must be numeric (integer or double).",
        name, colnames(x)[j]
      ), call. = FALSE)
    }
    if (anyNA(col)) {
      stop(sprintf(
        "`%s`: individual column '%s' must not contain NA.",
        name, colnames(x)[j]
      ), call. = FALSE)
    }
    if (any(col != -1 & col != 0 & col != 1)) {
      stop(sprintf(
        "`%s`: individual column '%s' must be coded -1/0/1.",
        name, colnames(x)[j]
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Pack a genotype vector of {-1,0,1} values into a raw vector
#'
#' Four calls per byte, 2 bits each (byte-coded value = genotype + 1, so
#' `0, 1, 2`), least-significant pair first. Padded with zero calls
#' (genotype 0) to a multiple of 4.
#'
#' @param v Integer vector of genotype values in `{-1, 0, 1}`.
#' @return A raw vector of `ceiling(length(v) / 4)` bytes.
#' @noRd
.pack_calls <- function(v) {
  n <- length(v)
  pad <- (4L - n %% 4L) %% 4L
  if (pad > 0L) v <- c(v, rep.int(0L, pad))
  m <- matrix(v + 1L, nrow = 4L)
  as.raw(m[1, ] + 4L * m[2, ] + 16L * m[3, ] + 64L * m[4, ])
}

#' Unpack a raw vector into a vector of {0,1,2} values
#'
#' @param r Raw vector produced by `.pack_calls()`.
#' @param n Number of calls to return.
#' @return Integer vector of length `n` with values in `{0, 1, 2}`.
#' @noRd
.unpack_calls <- function(r, n) {
  b <- as.integer(r)
  v <- c(rbind(
    b %% 4L,
    (b %/% 4L) %% 4L,
    (b %/% 16L) %% 4L,
    b %/% 64L
  ))
  v[seq_len(n)] - 1L
}

#' Pack a marker data.frame into format-1 list
#'
#' @param x A marker `data.frame` satisfying the contract checked by
#'   `.check_marker()`.
#' @param name `character(1)` dataset name, stored in the packed object.
#' @return A `list` with components `format`, `name`, `meta`, `ind`, `n_snp`,
#'   `n_ind`, `layout`, `geno`, `geno_mode`, `row_names` -- see package
#'   internal documentation ("packed file format, format 1").
#' @noRd
.pack <- function(x, name) {
  n_snp <- nrow(x)
  n_ind <- ncol(x) - 5L
  meta <- x[, 1:5]
  rownames(meta) <- NULL
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  ind <- colnames(x)[-(1:5)]

  gm <- as.matrix(x[, -(1:5), drop = FALSE])
  geno_mode <- storage.mode(gm)
  g_int <- as.integer(gm) # -1/0/1, column-major (SNP x individual)

  # layout "ind": column-major over SNP x individual (all SNPs of
  #   individual 1, then individual 2, ...) -- this is exactly gm's own
  #   column-major storage order.
  raw_ind <- .pack_calls(g_int)
  # layout "snp": row-major over SNP x individual (all individuals at SNP 1,
  #   then SNP 2, ...). gm is stored column-major, so transpose before
  #   flattening.
  g_snp <- as.integer(t(matrix(g_int, nrow = n_snp, ncol = n_ind)))
  raw_snp <- .pack_calls(g_snp)

  c_snp <- memCompress(raw_snp, type = "xz")
  c_ind <- memCompress(raw_ind, type = "xz")

  if (length(c_snp) <= length(c_ind)) {
    layout <- "snp"
    geno <- raw_snp
  } else {
    layout <- "ind"
    geno <- raw_ind
  }

  row_names <- if (.row_names_info(x) < 0) {
    NULL
  } else {
    as.character(rownames(x))
  }

  list(
    format = 1L,
    name = as.character(name),
    meta = meta,
    ind = ind,
    n_snp = as.integer(n_snp),
    n_ind = as.integer(n_ind),
    layout = layout,
    geno = geno,
    geno_mode = geno_mode,
    row_names = row_names
  )
}

#' Unpack a format-1 packed object into a marker data.frame
#'
#' @param p A list produced by `.pack()`.
#' @return The reconstructed marker `data.frame`.
#' @noRd
.unpack <- function(p) {
  if (is.null(p$format) || p$format != 1L) {
    stop(sprintf(
      "Unsupported packed format: %s (only format 1 is supported).",
      if (is.null(p$format)) "NULL" else as.character(p$format)
    ), call. = FALSE)
  }
  n_snp <- p$n_snp
  n_ind <- p$n_ind
  v <- .unpack_calls(p$geno, n_snp * n_ind)
  M <- matrix(v, nrow = n_snp, ncol = n_ind, byrow = (p$layout == "snp"))
  colnames(M) <- p$ind
  storage.mode(M) <- p$geno_mode

  out <- data.frame(p$meta, as.data.frame(M, check.names = FALSE),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  # A freshly constructed data.frame already has automatic (compact) row
  # names, matching .row_names_info(x) < 0; only override when the packed
  # object carried explicit row names.
  if (!is.null(p$row_names)) {
    rownames(out) <- p$row_names
  }
  out
}
