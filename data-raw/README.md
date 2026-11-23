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
  s1database_panels.R       builds the sample-only Arabidopsis_div and Human_div
  thin.R                    the shared thinning rule (sourced by samples.R and
                             s1database_panels.R)
  samples.R                 thins every built panel to the in-package sample
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
| `Arabidopsis_div` (sample only) | `s1database_panels.R` |
| `Human_div` (sample only)        | `s1database_panels.R` |
| `SugarBeet_div` (data-v2) | `zenodo_sugarbeet_chu.R`        |
| `Tomato_div` (data-v2)    | `mendeley_tomato_51k.R`         |
| `WheatHex_div` (data-v2)  | `dryad_wheat_zhanggao.R`        |
| `Salmon_div` (data-v2)    | `dryad_salmon_kijas.R`          |
| `Sheep_div` (data-v2)     | `dryad_sheep_ciani.R`           |
| `Potato_div` (data-v2)    | `dryad_potato_endelman.R`       |
| `Sugarcane_div` (data-v2) | `usda_sugarcane_ratoon.R`       |

Raw-only scripts with no data set yet (see `TODO.md`): `figshare_chicken_li.R`.

The catalog (`inst/extdata/catalog.csv`) has three kinds of rows: the bundled
maize panel, hosted panels (sample in the package, whole panel in a GitHub
release: tag `data-v1` for the first seven, `data-v2` for the seven added
later) and two **sample-only** data sets. `Arabidopsis_div` and `Human_div` are built by
`s1database_panels.R` from the processed pilot panels of the s1database
project (`data-panels/<panel>/pilot-v1/`: `variants.tsv`,
`dosages.simulation.tsv.gz`), thinned with the same rule as the other samples
(`thin.R`, shared with `samples.R`). markerdata hosts no whole panel for
them: their catalog rows have no `file`/`release`/`bytes`/`md5`, `source_url`
holds the original download, and `n_snp`/`n_ind`/`het` describe the sample
(computed by `catalog.R` from `inst/extdata/samples/<name>.rds`, since there
is no built `.rda`). `release.R` skips them. `source_url` is filled for every
row (empty for the maize and T3/Wheat rows, where no URL is on record).

`WheatHex_div` is the hexaploid bread-wheat panel, kept alongside the
tetraploid `Wheat_div` panel.

Two panels are not what they look like, and their help topics say so:
`Potato_div` is an autotetraploid stored **diploidized** (dosage lost), and
`Sugarcane_div` holds the diploid-style calls of its source (no dosage ever
existed) at *Sorghum bicolor* positions.

Each hosted row names its own release tag (`release` field in `catalog.R`).
`release.R` packs only rows that still have no checksum and never re-tags or
re-packs the others, so the published `data-v1` assets stay valid.

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
3. `Rscript data-raw/release.R` -- for every hosted catalog row with a
   built `.rda` and no checksum yet, validate it (`.check_marker()`), pack it (`.pack()`),
   round-trip-check it (`identical(.unpack(p), x)`), write
   `data-raw/release/<name>.rds`, and fill in that row's `file`, `release`,
   `bytes` and `md5` in `inst/extdata/catalog.csv`. Rows that already have a
   checksum are skipped (set `FORCE <- TRUE` at the top of the script to
   re-pack everything); rows with no built `.rda` are reported and left
   unchanged.
4. `Rscript data-raw/samples.R` (and `Rscript data-raw/s1database_panels.R`
   for the two sample-only data sets) -- rebuild the samples shipped inside the
   package (`inst/extdata/samples/<name>.rds`, 2,000 SNPs for
   every individual, tracked in git) from the same built `.rda` files. These are
   what `marker_data("<name>")` returns without `full = TRUE`.
5. Publish the packed assets, one release per tag, listing that tag's files
   (never `data-raw/release/*.rds`, which mixes tags). `data-v1` (first seven
   panels, already published):

   ```sh
   gh release create data-v1 --title "markerdata data v1" \
     --notes "Packed marker panels for markerdata 0.1.0" \
     data-raw/release/{SoySNP50K,EUChip60K,Sorghum,Rice,Cattle,Poplar,Wheat}_div.rds
   ```

   `data-v2` (seven panels, NOT published yet; do this before the next
   package release):

   ```sh
   gh release create data-v2 --title "markerdata data v2" \
     --notes "Packed marker panels for markerdata 0.1.0 (sugar beet, tomato, hexaploid wheat, Atlantic salmon, sheep, potato, sugarcane)" \
     data-raw/release/SugarBeet_div.rds data-raw/release/Tomato_div.rds \
     data-raw/release/WheatHex_div.rds data-raw/release/Salmon_div.rds \
     data-raw/release/Sheep_div.rds data-raw/release/Potato_div.rds \
     data-raw/release/Sugarcane_div.rds
   ```

   Until `data-v2` exists, `marker_data(name, full = TRUE)` for these seven
   warns and returns `NULL` (the sample still works). To add or replace
   assets in an existing tag, use `gh release upload <tag> <files> --clobber`.

6. Commit the updated `inst/extdata/catalog.csv` (and `R/sysdata.rda` if
   references changed). Do not commit `data-raw/built/` or
   `data-raw/release/` -- they are gitignored and reproducible from the
   scripts above.

## Checklist before a CRAN release

- **T3/Wheat (`Wheat_div`)**: public data under the Toronto Statement (cite
  the data source and T3; respect the producers' right to publish first
  analyses); USDA catalogs label it public domain and the same files are on
  GrainGenes. No formal license text exists. Optionally ask the T3 curators
  to confirm redistribution in writing.
- **SoyBase (`SoySNP50K_div`)**: U.S. public domain (USDA catalog label);
  cite the papers and write "SoyBase" with its URL.
- The `data-v2` panels are CC0-1.0 or CC-BY-4.0 (`Sugarcane_div` is a US
  Government work, public domain) and clear to redistribute with
  attribution.
- All other datasets (`EUChip60K_div`, `Sorghum_div`, `Rice_div`,
  `Cattle_div`, `Poplar_div`) are CC0 or CC-BY-4.0 and already clear to
  redistribute with attribution.
