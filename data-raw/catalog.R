#!/usr/bin/env Rscript
# data-raw/catalog.R
#
# Single source of truth for markerdata's dataset metadata and literature
# references. Run from the package root:
#
#   Rscript data-raw/catalog.R
#
# It writes two files:
#   inst/extdata/catalog.csv  -- one row per dataset, read at runtime by
#                                 marker_catalog() (and by marker_download() /
#                                 marker_cache_clear() to resolve file/release
#                                 names, sizes and checksums).
#   R/sysdata.rda             -- one object, `.marker_refs`, a named list of
#                                 utils::bibentry objects consumed by
#                                 marker_references().
#
# n_snp, n_ind, het and chr_labels are NOT typed by hand: they are computed
# here from the built objects in data-raw/built/<name>.rda (see
# data-raw/<script>.R for how those are produced). bytes and md5 are left
# blank -- data-raw/release.R fills them in once the packed .rds assets are
# built and ready to upload as a GitHub release. Re-running this script
# preserves any bytes/md5 already present in inst/extdata/catalog.csv, so it
# never wipes checksums that release.R has already computed.
#
# WheatHex_div (the hexaploid bread-wheat panel built by
# data-raw/dryad_wheat_zhanggao.R) is intentionally NOT included below: it is
# a second, not-yet-released wheat option (see that script's header).

suppressWarnings(suppressMessages({
  ok_utils <- requireNamespace("utils", quietly = TRUE)
}))
stopifnot(ok_utils)

built_dir <- "data-raw/built"
if (!dir.exists("data-raw") || !file.exists("DESCRIPTION")) {
  stop("Run this script from the markerdata package root (data-raw/catalog.R).")
}

## ---------------------------------------------------------------------
## 1. Dataset metadata (everything that is NOT derived from the built .rda)
## ---------------------------------------------------------------------

p <- function(family, given) person(given, family)

datasets <- list(
  list(
    name = "SNP55K_maize282_maf04",
    species = "Zea mays",
    common_name = "maize",
    panel = "282 maize inbred association panel, Illumina MaizeSNP50 (55K) array, MAF>0.4 subset",
    ploidy = "diploid",
    mating = "inbred",
    license = "Panzea public data",
    source = "Panzea, maize282set MaizeSNP50 genotypes (MAF>0.4 subset)",
    source_doi = "",
    reference = "Cook et al. (2012) Plant Physiology 158(2):824-834",
    doi = "10.1104/pp.111.185033",
    bundled = TRUE,
    refs = c(
      bibentry(bibtype = "Article",
        key = "cook2012maize282",
        title = "Genetic architecture of maize kernel composition in the nested association mapping and inbred association panels",
        author = c(p("Cook", "J. P."), p("McMullen", "M. D."), p("Holland", "J. B."),
                   p("Tian", "F."), p("Bradbury", "P."), p("Ross-Ibarra", "J."),
                   p("Buckler", "E. S."), p("Flint-Garcia", "S. A.")),
        journal = "Plant Physiology", year = 2012, volume = 158, number = 2,
        pages = "824-834", doi = "10.1104/pp.111.185033")
    )
  ),
  list(
    name = "SoySNP50K_div",
    species = "Glycine max",
    common_name = "soybean",
    panel = "USDA Soybean Germplasm Collection diversity subset, SoySNP50K iSelect BeadChip",
    ploidy = "diploid",
    mating = "inbred",
    license = "Public USDA/SoyBase data (cite)",
    source = "SoyBase, glyma.Wm82.gnm2.div.Song_Hyten_2015.hmp.gz (HapMap)",
    source_doi = "",
    reference = "Song et al. (2013) PLoS ONE 8(1):e54985; Song et al. (2015) G3 5(10):1999-2006",
    doi = "10.1371/journal.pone.0054985; 10.1534/g3.115.019000",
    bundled = FALSE,
    refs = c(
      bibentry(bibtype = "Article",
        key = "song2013soysnp50k",
        title = "Development and evaluation of SoySNP50K, a high-density genotyping array for soybean",
        author = c(p("Song", "Q."), p("Hyten", "D. L."), p("Jia", "G."),
                   p("Quigley", "C. V."), p("Fickus", "E. W."), p("Nelson", "R. L."),
                   p("Cregan", "P. B.")),
        journal = "PLoS ONE", year = 2013, volume = 8, number = 1, pages = "e54985",
        doi = "10.1371/journal.pone.0054985"),
      bibentry(bibtype = "Article",
        key = "song2015soyfingerprint",
        title = "Fingerprinting soybean germplasm and its utility in genomic research",
        author = c(p("Song", "Q."), p("Hyten", "D. L."), p("Jia", "G."),
                   p("Quigley", "C. V."), p("Fickus", "E. W."), p("Nelson", "R. L."),
                   p("Cregan", "P. B.")),
        journal = "G3: Genes, Genomes, Genetics", year = 2015, volume = 5, number = 10,
        pages = "1999-2006", doi = "10.1534/g3.115.019000"),
      bibentry(bibtype = "Misc",
        key = "soybase2015div",
        title = "SoySNP50K diversity panel HapMap (Wm82.gnm2.div.Song_Hyten_2015)",
        author = person(family = "SoyBase"), year = 2015,
        url = "https://data.soybase.org/Glycine/max/diversity/Wm82.gnm2.div.Song_Hyten_2015/")
    )
  ),
  list(
    name = "EUChip60K_div",
    species = "Eucalyptus grandis x E. urophylla",
    common_name = "eucalyptus",
    panel = "synthetic breeding population, EUChip60K Illumina Infinium array",
    ploidy = "diploid",
    mating = "outbred",
    license = "CC0-1.0",
    source = "Dryad, SNP_Data.txt",
    source_doi = "10.5061/dryad.ms580",
    reference = "Resende et al. (2017) Heredity 119(4):245-255; Silva-Junior et al. (2015) New Phytologist 206(4):1527-1540",
    doi = "10.1038/hdy.2017.37; 10.1111/nph.13322",
    bundled = FALSE,
    refs = c(
      bibentry(bibtype = "Article",
        key = "resende2017euchip",
        title = "Assessing the expected response to genomic selection of individuals and families in Eucalyptus breeding with an additive-dominant model",
        author = c(p("Resende", "R. T."), p("Resende", "M. D. V."), p("Silva", "F. F."),
                   p("Azevedo", "C. F."), p("Takahashi", "E. K."), p("Silva-Junior", "O. B."),
                   p("Grattapaglia", "D.")),
        journal = "Heredity", year = 2017, volume = 119, number = 4, pages = "245-255",
        doi = "10.1038/hdy.2017.37"),
      bibentry(bibtype = "Article",
        key = "silvajunior2015euchip",
        title = "A flexible multispecies genome-wide 60K SNP chip developed from pooled resequencing of 240 Eucalyptus tree genomes across 12 species",
        author = c(p("Silva-Junior", "O. B."), p("Faria", "D. A."), p("Grattapaglia", "D.")),
        journal = "New Phytologist", year = 2015, volume = 206, number = 4,
        pages = "1527-1540", doi = "10.1111/nph.13322"),
      bibentry(bibtype = "Misc",
        key = "dryad2017euchip",
        title = "Data from: Assessing the expected response to genomic selection of individuals and families in Eucalyptus breeding with an additive-dominant model",
        author = c(p("Resende", "R. T."), p("Resende", "M. D. V."), p("Silva", "F. F."),
                   p("Azevedo", "C. F."), p("Takahashi", "E. K."), p("Silva-Junior", "O. B."),
                   p("Grattapaglia", "D.")),
        year = 2017, doi = "10.5061/dryad.ms580", note = "Dryad Digital Repository")
    )
  ),
  list(
    name = "Sorghum_div",
    species = "Sorghum bicolor",
    common_name = "sorghum",
    panel = "WEST Diversity Panel, genotyping-by-sequencing (GBS)",
    ploidy = "diploid",
    mating = "inbred",
    license = "CC-BY-4.0",
    source = "Zenodo, WEST Diversity Panel GBS (rawGBS_imputed_850x5mi.vcf.gz)",
    source_doi = "10.5281/zenodo.5019227",
    reference = "Ferguson et al. (2021) Plant Physiology 187(3):1481-1500",
    doi = "10.1093/plphys/kiab346",
    bundled = FALSE,
    refs = c(
      bibentry(bibtype = "Article",
        key = "ferguson2021sorghumwue",
        title = "Machine learning-enabled phenotyping for GWAS and TWAS of WUE traits in 869 field-grown sorghum accessions",
        author = c(p("Ferguson", "J. N."), p("Fernandes", "S. B."), p("Monier", "B."),
                   p("Miller", "N. D."), p("Allen", "D."), p("Dmitrieva", "A."),
                   p("Schmuker", "P."), p("Lozano", "R."), p("Valluru", "R."),
                   p("Buckler", "E. S."), p("Gore", "M. A."), p("Brown", "P. J."),
                   p("Spalding", "E. P."), p("Leakey", "A. D. B.")),
        journal = "Plant Physiology", year = 2021, volume = 187, number = 3,
        pages = "1481-1500", doi = "10.1093/plphys/kiab346"),
      bibentry(bibtype = "Misc",
        key = "ferguson2020westgbs",
        title = "WEST Diversity Panel GBS information",
        author = c(p("Ferguson", "J. N."), p("Fernandes", "S. B."), p("Monier", "B.")),
        year = 2020, doi = "10.5281/zenodo.5019227", note = "Zenodo (CC-BY-4.0)")
    )
  ),
  list(
    name = "Rice_div",
    species = "Oryza sativa",
    common_name = "rice",
    panel = "IRRI elite tropical rice (indica / indica-admixed) genomic-selection panel, GBS",
    ploidy = "diploid",
    mating = "inbred",
    license = "CC0-1.0",
    source = "Dryad, MET_crfilt_.75_allchrom.hmp.txt",
    source_doi = "10.5061/dryad.7369p",
    reference = "Spindel et al. (2015) PLoS Genetics 11(2):e1004982",
    doi = "10.1371/journal.pgen.1004982",
    bundled = FALSE,
    refs = c(
      bibentry(bibtype = "Article",
        key = "spindel2015rice",
        title = "Genomic Selection and Association Mapping in Rice (Oryza sativa): Effect of Trait Genetic Architecture, Training Population Composition, Marker Number and Statistical Model on Accuracy of Rice Genomic Selection in Elite, Tropical Rice Breeding Lines",
        author = c(p("Spindel", "J."), p("Begum", "H."), p("Akdemir", "D."),
                   p("Virk", "P."), p("Collard", "B."), p("Redona", "E."),
                   p("Atlin", "G."), p("Jannink", "J-L."), p("McCouch", "S. R.")),
        journal = "PLoS Genetics", year = 2015, volume = 11, number = 2,
        pages = "e1004982", doi = "10.1371/journal.pgen.1004982"),
      bibentry(bibtype = "Misc",
        key = "dryad2015rice",
        title = "Data from: Genomic selection and association mapping in rice (Oryza sativa)",
        author = c(p("Spindel", "J."), p("Begum", "H."), p("Akdemir", "D."),
                   p("Virk", "P."), p("Collard", "B."), p("Redona", "E."),
                   p("Atlin", "G."), p("Jannink", "J-L."), p("McCouch", "S. R.")),
        year = 2015, doi = "10.5061/dryad.7369p", note = "Dryad Digital Repository")
    )
  ),
  list(
    name = "Cattle_div",
    species = "Bos taurus",
    common_name = "cattle",
    panel = "worldwide cattle diversity panel (134 breeds), Illumina BovineSNP50 array",
    ploidy = "diploid",
    mating = "outbred",
    license = "CC0-1.0",
    source = "Dryad, WorldWideCattleDryad.zip (PLINK .ped/.map)",
    source_doi = "10.5061/dryad.th092",
    reference = "Decker et al. (2014) PLoS Genetics 10(3):e1004254",
    doi = "10.1371/journal.pgen.1004254",
    bundled = FALSE,
    refs = c(
      bibentry(bibtype = "Article",
        key = "decker2014cattle",
        title = "Worldwide Patterns of Ancestry, Divergence, and Admixture in Domesticated Cattle",
        author = c(p("Decker", "J. E."), p("McKay", "S. D."), p("Rolf", "M. M."),
                   p("Kim", "J."), p("Molina Alcala", "A."), p("Sonstegard", "T. S."),
                   p("Hanotte", "O."), p("Taylor", "J. F.")),
        journal = "PLoS Genetics", year = 2014, volume = 10, number = 3,
        pages = "e1004254", doi = "10.1371/journal.pgen.1004254"),
      bibentry(bibtype = "Misc",
        key = "dryad2014cattle",
        title = "Data from: Worldwide patterns of ancestry, divergence, and admixture in domesticated cattle",
        author = c(p("Decker", "J. E."), p("McKay", "S. D."), p("Rolf", "M. M."),
                   p("Kim", "J."), p("Molina Alcala", "A."), p("Sonstegard", "T. S."),
                   p("Hanotte", "O."), p("Taylor", "J. F.")),
        year = 2014, doi = "10.5061/dryad.th092", note = "Dryad Digital Repository")
    )
  ),
  list(
    name = "Poplar_div",
    species = "Populus trichocarpa",
    common_name = "poplar",
    panel = "landscape-genomics panel, PoplarConsortium-ORNL 34K SNP array",
    ploidy = "diploid",
    mating = "outbred",
    license = "CC0-1.0",
    source = "Dryad, FOR_DRYAD_TandB_fwd_june2013_434.txt",
    source_doi = "10.5061/dryad.7s848",
    reference = "Geraldes et al. (2014) Evolution 68(11):3260-3280; Geraldes et al. (2013) Molecular Ecology Resources 13(2):306-323",
    doi = "10.1111/evo.12497; 10.1111/1755-0998.12056",
    bundled = FALSE,
    refs = c(
      bibentry(bibtype = "Article",
        key = "geraldes2014poplar",
        title = "Landscape genomics of Populus trichocarpa: the role of hybridization, limited gene flow, and natural selection in shaping patterns of population structure",
        author = c(p("Geraldes", "A."), p("Farzaneh", "N."), p("Grassa", "C. J."),
                   p("McKown", "A. D."), p("Guy", "R. D."), p("Mansfield", "S. D."),
                   p("Douglas", "C. J."), p("Cronk", "Q. C. B.")),
        journal = "Evolution", year = 2014, volume = 68, number = 11,
        pages = "3260-3280", doi = "10.1111/evo.12497"),
      bibentry(bibtype = "Article",
        key = "geraldes2013poplararray",
        title = "A 34K SNP genotyping array for Populus trichocarpa: design, application to the study of natural populations and transferability to other Populus species",
        author = c(p("Geraldes", "A."), p("DiFazio", "S. P."), p("Slavov", "G. T.")),
        journal = "Molecular Ecology Resources", year = 2013, volume = 13, number = 2,
        pages = "306-323", doi = "10.1111/1755-0998.12056"),
      bibentry(bibtype = "Misc",
        key = "dryad2013poplar",
        title = "Data from: Landscape genomics of Populus trichocarpa",
        author = c(p("Geraldes", "A."), p("Farzaneh", "N."), p("Grassa", "C. J."),
                   p("McKown", "A. D."), p("Guy", "R. D."), p("Mansfield", "S. D."),
                   p("Douglas", "C. J."), p("Cronk", "Q. C. B.")),
        year = 2013, doi = "10.5061/dryad.7s848", note = "Dryad Digital Repository")
    )
  ),
  list(
    name = "Wheat_div",
    species = "Triticum turgidum",
    common_name = "wheat",
    panel = "Tetraploid wheat Global Collection (TGC), 90K iSelect SNP array",
    ploidy = "allotetraploid (AABB, disomic; coded per subgenome)",
    mating = "inbred",
    license = "T3/Wheat public data (cite; terms to confirm)",
    source = "T3/Wheat (Triticeae Toolbox), ..._TGC-GTC_collection_SNP90K_mapped_filtered.vcf",
    source_doi = "",
    reference = "Maccaferri et al. (2019) Nature Genetics 51:885-895",
    doi = "10.1038/s41588-019-0381-3",
    bundled = FALSE,
    refs = c(
      bibentry(bibtype = "Article",
        key = "maccaferri2019durum",
        title = "Durum wheat genome highlights past domestication signatures and future improvement targets",
        author = c(p("Maccaferri", "M."), p("Harris", "N. S."), p("Twardziok", "S. O."),
                   p("Cattivelli", "L."), p("Tuberosa", "R.")),
        journal = "Nature Genetics", year = 2019, volume = 51, pages = "885-895",
        doi = "10.1038/s41588-019-0381-3"),
      bibentry(bibtype = "Misc",
        key = "t3wheat_tgc",
        title = "Tetraploid wheat Global Collection (TGC), SNP90K genotypes",
        author = person(family = "Triticeae Toolbox (T3/Wheat)"), year = 2020,
        note = "T3/Wheat public data; redistribution terms to confirm")
    )
  )
)

names(datasets) <- vapply(datasets, `[[`, character(1), "name")

## ---------------------------------------------------------------------
## 2. Chromosome-label range, computed from the built object's `chr` column
## ---------------------------------------------------------------------

.lcp <- function(x) {
  x <- x[nzchar(x)]
  if (length(x) <= 1L) return("")
  minlen <- min(nchar(x))
  if (minlen == 0L) return("")
  i <- 0L
  repeat {
    if (i >= minlen) break
    ch <- substr(x[1L], i + 1L, i + 1L)
    if (!all(substr(x, i + 1L, i + 1L) == ch)) break
    i <- i + 1L
  }
  substr(x[1L], 1L, i)
}

chr_range <- function(chr) {
  u <- unique(as.character(chr))
  is_scaf <- grepl("scaffold", u, ignore.case = TRUE)
  main <- u[!is_scaf]
  n_scaf <- sum(is_scaf)
  if (length(main) == 0L) main <- u

  m <- regexpr("[0-9]+[A-Za-z]*$", main)
  has_tail <- m > 0L
  tail <- rep("", length(main))
  tail[has_tail] <- regmatches(main, m)
  prefix <- main
  prefix[has_tail] <- substr(main[has_tail], 1L, m[has_tail] - 1L)

  dm <- regmatches(tail, regexec("^([0-9]+)([A-Za-z]*)$", tail))
  digits <- vapply(dm, function(x) if (length(x) == 3L) suppressWarnings(as.integer(x[2])) else NA_integer_, integer(1))
  suffix <- vapply(dm, function(x) if (length(x) == 3L) x[3] else "", character(1))

  if (all(has_tail) && length(unique(prefix)) == 1L && !anyNA(digits)) {
    # common literal prefix across the (main) labels and the scaffold labels,
    # if any, is stripped for a cleaner display (e.g. an assembly-name prefix)
    common <- if (n_scaf > 0L) .lcp(c(prefix[1L], .lcp(u[is_scaf]))) else ""
    disp <- if (nchar(common) > 0L) substring(main, nchar(common) + 1L) else main
    ord <- order(digits, suffix)
    rng <- paste0(disp[ord[1L]], "-", disp[ord[length(ord)]])
  } else {
    s <- sort(main)
    rng <- paste0(s[1L], "-", s[length(s)])
  }
  if (n_scaf > 0L) {
    rng <- paste0(rng, " (+", n_scaf, " scaffold", if (n_scaf != 1L) "s" else "", ")")
  }
  rng
}

## ---------------------------------------------------------------------
## 3. Compute n_snp / n_ind / het / chr_labels from data-raw/built/<name>.rda
## ---------------------------------------------------------------------

computed <- lapply(datasets, function(d) {
  rda <- file.path(built_dir, paste0(d$name, ".rda"))
  if (!file.exists(rda)) {
    warning(sprintf("catalog.R: %s not found; leaving n_snp/n_ind/het/chr_labels blank for %s.",
                    rda, d$name), call. = FALSE)
    return(list(n_snp = NA_integer_, n_ind = NA_integer_, het = NA_real_, chr_labels = ""))
  }
  e <- new.env()
  nm <- load(rda, envir = e)
  obj <- get(nm[1L], envir = e)
  n_snp <- nrow(obj)
  n_ind <- ncol(obj) - 5L
  G <- as.matrix(obj[, 6:ncol(obj), drop = FALSE])
  het <- round(sum(G == 0L, na.rm = TRUE) / length(G), 3)
  list(n_snp = n_snp, n_ind = n_ind, het = het, chr_labels = chr_range(obj$chr))
})
names(computed) <- names(datasets)

## ---------------------------------------------------------------------
## 4. Preserve bytes/md5 from an existing catalog.csv, if present
## ---------------------------------------------------------------------

catalog_path <- "inst/extdata/catalog.csv"
prev_bytes <- prev_md5 <- setNames(rep("", length(datasets)), names(datasets))
if (file.exists(catalog_path)) {
  prev <- utils::read.csv(catalog_path, stringsAsFactors = FALSE, na.strings = "",
                          colClasses = "character")
  for (i in seq_len(nrow(prev))) {
    nm <- prev$name[i]
    if (nm %in% names(datasets)) {
      if (!is.na(prev$bytes[i])) prev_bytes[[nm]] <- prev$bytes[i]
      if (!is.na(prev$md5[i]))   prev_md5[[nm]]   <- prev$md5[i]
    }
  }
}

## ---------------------------------------------------------------------
## 5. Assemble inst/extdata/catalog.csv
## ---------------------------------------------------------------------

rows <- lapply(datasets, function(d) {
  cm <- computed[[d$name]]
  data.frame(
    name = d$name,
    species = d$species,
    common_name = d$common_name,
    panel = d$panel,
    ploidy = d$ploidy,
    mating = d$mating,
    n_snp = cm$n_snp,
    n_ind = cm$n_ind,
    het = cm$het,
    chr_labels = cm$chr_labels,
    license = d$license,
    source = d$source,
    source_doi = d$source_doi,
    reference = d$reference,
    doi = d$doi,
    bundled = d$bundled,
    file = if (d$bundled) "" else paste0(d$name, ".rds"),
    release = if (d$bundled) "" else "data-v1",
    bytes = prev_bytes[[d$name]],
    md5 = prev_md5[[d$name]],
    stringsAsFactors = FALSE, check.names = FALSE
  )
})
catalog <- do.call(rbind, rows)
rownames(catalog) <- NULL

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(catalog, file = catalog_path, row.names = FALSE, na = "")
message("Wrote ", catalog_path, " (", nrow(catalog), " rows).")

## ---------------------------------------------------------------------
## 6. Assemble R/sysdata.rda (.marker_refs)
## ---------------------------------------------------------------------

.marker_refs <- lapply(datasets, `[[`, "refs")
names(.marker_refs) <- names(datasets)

dir.create("R", recursive = TRUE, showWarnings = FALSE)
save(.marker_refs, file = "R/sysdata.rda", compress = "xz", version = 2)
message("Wrote R/sysdata.rda (.marker_refs, ", length(.marker_refs), " entries).")

## ---------------------------------------------------------------------
## 7. Self-check: every bibentry prints and converts to BibTeX
## ---------------------------------------------------------------------

for (nm in names(.marker_refs)) {
  x <- .marker_refs[[nm]]
  base::print(x, style = "text")
  bt <- utils::toBibtex(x)
  stopifnot(is.character(bt), length(bt) > 0L)
}

message("catalog.R: OK.")
