# data-raw/

This folder builds and packages the marker panels distributed by markerdata.
Nothing here is part of the installed package except `catalog.R`'s output
(`inst/extdata/catalog.csv`); the rest is gitignored and exists only to
(re)produce that output and the GitHub release assets that `marker_data()`
downloads at runtime.

## Layout

```
data-raw/
  catalog.R                 single source of truth for dataset metadata and
                             literature references -> inst/extdata/catalog.csv
                             and R/sysdata.rda (.marker_refs)
  release.R                 packs data-raw/built/*.rda into the release
                             format and stages it under data-raw/release/,
                             then fills in catalog.csv's file/release/bytes/md5
  <species>_<source>.R      one build script per dataset (see below)
  README.md                 this file

  <species>/                raw input files for one build script, e.g.
  soybase_cache/             data-raw/sorghum/, data-raw/wheat/,
                             data-raw/soybase_cache/, ... -- gitignored,
                             fetched by hand (see "Raw inputs" below)

  built/                     gitignored. One <name>.rda per dataset, written
                             by its build script; the working, uncompressed
                             (well, xz-compressed .rda, but not packed)
                             data.frame

  release/                   gitignored. One <name>.rds per dataset, written
                             by release.R; the packed, xz-compressed asset
                             that actually gets uploaded to a GitHub release
```

Only `catalog.R`'s two outputs (`inst/extdata/catalog.csv` and
`R/sysdata.rda`) and the build/release scripts themselves are tracked in
git. `data-raw/built/`, `data-raw/release/`, and every `data-raw/<species>/`
raw-input folder are listed in `.gitignore` and never committed -- rebuild
them locally when you need them.

## Datasets and their build scripts

| dataset                 | build script                    |
|--------------------------|---------------------------------|
| `SNP55K_maize282_maf04`  | (bundled; ships from `simplePHENOTYPES`-style input, see `data/`) |
| `SoySNP50K_div`          | `soybase_soysnp50k.R`           |
| `EUChip60K_div`          | `dryad_euchip60k.R`             |
| `Sorghum_div`            | `sorghum_wue.R`                 |
| `Rice_div`               | `dryad_rice_spindel.R`          |
| `Cattle_div`             | `dryad_cattle_decker.R`         |
| `Poplar_div`             | `dryad_poplar_geraldes.R`       |
| `Wheat_div`              | `t3_wheat_tgc.R`                |
| `WheatHex_div` (not yet released) | `dryad_wheat_zhanggao.R` |

`WheatHex_div` is a second, hexaploid bread-wheat option kept alongside the
tetraploid `Wheat_div` panel. It builds successfully but has no entry in
`inst/extdata/catalog.csv` yet and is not downloadable with `marker_data()`
-- add a `catalog.R` entry and run `release.R` for it once it is ready.

## Raw inputs

Each build script's raw input is a large file from a public repository
(Dryad, Zenodo, SoyBase, T3/Wheat) that is **not** fetched automatically.
**Dryad in particular blocks automated/scripted downloads**, so every raw
file listed in a build script's header must be downloaded by hand (through
the browser, or a signed-in `curl`/API session where the repository
requires one) into the corresponding `data-raw/<species>/` (or
`data-raw/soybase_cache/`) folder before running that script. See each
script's header comment for the exact source URL/DOI, file name, and
license.

## (Re)building a dataset

1. Fetch the raw input(s) for the dataset into its `data-raw/<species>/`
   folder (see the script's header comment for the source).
2. Run the build script from the package root, e.g.:

   ```sh
   Rscript data-raw/dryad_rice_spindel.R
   ```

   Each script takes optional positional arguments (input path, subset
   size, output path); run it with no arguments to use the defaults, which
   write `data-raw/built/<name>.rda`. Every script sanity-checks its output
   with `stopifnot(... %in% c(-1L, 0L, 1L))` and
   `simplePHENOTYPES::as_population()` before saving.
3. Re-run `Rscript data-raw/catalog.R` so `inst/extdata/catalog.csv` picks
   up the freshly built object's `n_snp`, `n_ind`, `het` and `chr_labels`
   (computed from `data-raw/built/<name>.rda`, not typed by hand).

## The packed release format (format 1)

`marker_data()` downloads and caches a small, packed `.rds` per dataset,
not the plain `data.frame` in `data-raw/built/`. Packing (`R/pack.R`,
`.pack()` / `.unpack()`) turns a marker `data.frame` (columns `snp`,
`allele`, `chr`, `pos`, `cm`, then one `-1`/`0`/`1` column per individual)
into a compact `list` with these fields:

| field       | contents |
|-------------|----------|
| `format`    | packed-format version (currently `1L`) |
| `name`      | dataset name (character) |
| `meta`      | the 5 map columns (`snp`, `allele`, `chr`, `pos`, `cm`) as a data.frame |
| `ind`       | individual (column) names, in order |
| `n_snp`     | number of SNPs (rows) |
| `n_ind`     | number of individuals (columns) |
| `layout`    | `"snp"` (row-major: all individuals at SNP 1, then SNP 2, ...) or `"ind"` (column-major: all SNPs of individual 1, then individual 2, ...) -- whichever compresses smaller is kept |
| `geno`      | the genotype calls packed 2 bits per call (byte value = genotype + 1, so `0`/`1`/`2` for `-1`/`0`/`1`), 4 calls per byte with the least-significant pair first, padded with zero calls to a multiple of 4, then `xz`-compressed |
| `geno_mode` | original storage mode of the genotype matrix (`"integer"` or `"double"`), so `.unpack()` restores the same type |
| `row_names` | explicit row names if the source data.frame had any, else `NULL` |

`.pack(x, name)` builds this list (choosing whichever of the `"snp"` /
`"ind"` layouts compresses smaller) and `.unpack(p)` reconstructs the exact
original data.frame; `identical(.unpack(.pack(x, name)), x)` holds for every
released dataset (asserted by `release.R`).

## Release workflow

1. `Rscript data-raw/catalog.R` -- refresh `inst/extdata/catalog.csv` and
   `R/sysdata.rda` from the current metadata and the objects in
   `data-raw/built/`.
2. Build (or rebuild) whichever datasets changed, as above.
3. `Rscript data-raw/release.R` -- for every non-bundled catalog row with a
   built `.rda`, validate it (`.check_marker()`), pack it (`.pack()`),
   round-trip-check it (`identical(.unpack(p), x)`), write
   `data-raw/release/<name>.rds`, and fill in that row's `file`, `release`,
   `bytes` and `md5` in `inst/extdata/catalog.csv`. Rows with no built
   `.rda` are reported and left unchanged.
4. Publish the packed assets as a GitHub release (first time):

   ```sh
   gh release create data-v1 --title "markerdata data v1" \
     --notes "Packed marker panels for markerdata 0.1.0" \
     data-raw/release/*.rds
   ```

   or update an existing release with newly packed assets:

   ```sh
   gh release upload data-v1 data-raw/release/*.rds --clobber
   ```

5. Commit the updated `inst/extdata/catalog.csv` (and `R/sysdata.rda` if
   references changed). Do not commit `data-raw/built/` or
   `data-raw/release/` -- they are gitignored and reproducible from the
   scripts above.

## Checklist before a CRAN release

- **T3/Wheat (`Wheat_div`)**: redistribution terms for T3/Wheat
  (Triticeae Toolbox) public data have not been confirmed. Get an explicit
  answer from T3/Wheat before shipping this dataset's release asset
  publicly, or drop it from the catalog if redistribution is not allowed.
- **SoyBase (`SoySNP50K_div`)**: likewise confirm SoyBase's redistribution
  terms for the SoySNP50K diversity HapMap before a CRAN release; the
  license field is deliberately left as "cite" rather than a specific
  license identifier until this is confirmed.
- All other datasets (`EUChip60K_div`, `Sorghum_div`, `Rice_div`,
  `Cattle_div`, `Poplar_div`) are CC0 or CC-BY-4.0 and already clear to
  redistribute with attribution.
