
# markerdata: Marker Data Sets for Breeding Simulation

[![R-CMD-check](https://github.com/fernandes-lab/markerdata/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fernandes-lab/markerdata/actions/workflows/R-CMD-check.yaml)

markerdata provides curated single-nucleotide polymorphism (SNP) marker
panels drawn from public genomic resources across several species
(maize, soybean, eucalyptus, sorghum, rice, cattle, poplar, wheat). Each
panel is a plain data frame with five marker-metadata columns (`snp`,
`allele`, `chr`, `pos`, `cm`) followed by one `-1`/`0`/`1` genotype
column per individual, so it can be used directly as founder genotypes
for breeding-program simulation with `simplePHENOTYPES::as_population()`
or with breedingDesigner.

One small maize panel is bundled with the package and loads instantly;
the other seven panels are downloaded on first use with `marker_data()`
and cached locally, so later calls are instant too. Altogether the seven
downloadable panels total about 7.5 MB, and none of them bloat the
installed package.

## Installation

``` r
# install.packages("markerdata")  # once on CRAN
remotes::install_github("fernandes-lab/markerdata")
```

## Usage

List every data set in the catalog and see which ones are already
available on this machine:

``` r
library(markerdata)
marker_catalog()[, c("name", "common_name", "n_snp", "n_ind", "het", "mating", "available")]
#>                    name common_name n_snp n_ind   het  mating available
#> 1 SNP55K_maize282_maf04       maize 10650   280 0.004  inbred      TRUE
#> 2         SoySNP50K_div     soybean 12000   300 0.000  inbred     FALSE
#> 3         EUChip60K_div  eucalyptus 12000   768 0.463 outbred     FALSE
#> 4           Sorghum_div     sorghum 12000   842 0.051  inbred     FALSE
#> 5              Rice_div        rice 12000   370 0.033  inbred     FALSE
#> 6            Cattle_div      cattle 12000  1543 0.287 outbred     FALSE
#> 7            Poplar_div      poplar 12000   434 0.318 outbred     FALSE
#> 8             Wheat_div       wheat 12000  1854 0.076  inbred     FALSE
```

The bundled maize panel loads without any download:

``` r
head(marker_data("SNP55K_maize282_maf04")[, 1:8])
#>           snp allele chr     pos        cm 4226 4722 33-16
#> 1 ss196422159    A/G   1  379844 0.0000000   -1   -1    -1
#> 2 ss196422171    G/A   1  613257 0.2482394   -1   -1    -1
#> 3 ss196422173    G/A   1  659354 0.2972609    1    1    -1
#> 4 ss196422186    G/C   1  992572 0.6515820    1   -1    -1
#> 5 ss196500940    G/A   1 2044264 1.7694445    1    1     1
#> 6 ss196422234    A/G   1 2044555 1.7697537   -1   -1    -1
```

A remote panel is downloaded on first use and cached afterward. Its
`-1`/`0`/`1` genotype coding is exactly what
`simplePHENOTYPES::as_population()` expects for founder genotypes:

``` r
geno <- marker_data("Cattle_div")
pop <- simplePHENOTYPES::as_population(geno, individuals = 1:100)
```

## Data sets

| Data set | Species | Panel | SNPs x individuals | Heterozygous calls | Mating | License | Download | Reference |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| SNP55K_maize282_maf04 | Zea mays | 282 maize inbred association panel, Illumina MaizeSNP50 (55K) array, MAF\>0.4 subset | 10,650 x 280 | 0.4% | inbred | Panzea public data | bundled | Cook et al. (2012) Plant Physiology 158(2):824-834 [10.1104/pp.111.185033](https://doi.org/10.1104/pp.111.185033) |
| SoySNP50K_div | Glycine max | USDA Soybean Germplasm Collection diversity subset, SoySNP50K iSelect BeadChip | 12,000 x 300 | 0.0% | inbred | Public USDA/SoyBase data (cite) | 0.3 MB | Song et al. (2013) PLoS ONE 8(1):e54985 [10.1371/journal.pone.0054985](https://doi.org/10.1371/journal.pone.0054985); Song et al. (2015) G3 5(10):1999-2006 [10.1534/g3.115.019000](https://doi.org/10.1534/g3.115.019000) |
| EUChip60K_div | Eucalyptus grandis x E. urophylla | synthetic breeding population, EUChip60K Illumina Infinium array | 12,000 x 768 | 46.3% | outbred | CC0-1.0 | 0.5 MB | Resende et al. (2017) Heredity 119(4):245-255 [10.1038/hdy.2017.37](https://doi.org/10.1038/hdy.2017.37); Silva-Junior et al. (2015) New Phytologist 206(4):1527-1540 [10.1111/nph.13322](https://doi.org/10.1111/nph.13322) |
| Sorghum_div | Sorghum bicolor | WEST Diversity Panel, genotyping-by-sequencing (GBS) | 12,000 x 842 | 5.1% | inbred | CC-BY-4.0 | 0.9 MB | Ferguson et al. (2021) Plant Physiology 187(3):1481-1500 [10.1093/plphys/kiab346](https://doi.org/10.1093/plphys/kiab346) |
| Rice_div | Oryza sativa | IRRI elite tropical rice (indica / indica-admixed) genomic-selection panel, GBS | 12,000 x 370 | 3.3% | inbred | CC0-1.0 | 0.3 MB | Spindel et al. (2015) PLoS Genetics 11(2):e1004982 [10.1371/journal.pgen.1004982](https://doi.org/10.1371/journal.pgen.1004982) |
| Cattle_div | Bos taurus | worldwide cattle diversity panel (134 breeds), Illumina BovineSNP50 array | 12,000 x 1,543 | 28.7% | outbred | CC0-1.0 | 3.0 MB | Decker et al. (2014) PLoS Genetics 10(3):e1004254 [10.1371/journal.pgen.1004254](https://doi.org/10.1371/journal.pgen.1004254) |
| Poplar_div | Populus trichocarpa | landscape-genomics panel, PoplarConsortium-ORNL 34K SNP array | 12,000 x 434 | 31.8% | outbred | CC0-1.0 | 1.0 MB | Geraldes et al. (2014) Evolution 68(11):3260-3280 [10.1111/evo.12497](https://doi.org/10.1111/evo.12497); Geraldes et al. (2013) Molecular Ecology Resources 13(2):306-323 [10.1111/1755-0998.12056](https://doi.org/10.1111/1755-0998.12056) |
| Wheat_div | Triticum turgidum | Tetraploid wheat Global Collection (TGC), 90K iSelect SNP array | 12,000 x 1,854 | 7.6% | inbred | T3/Wheat public data (cite; terms to confirm) | 1.6 MB | Maccaferri et al. (2019) Nature Genetics 51:885-895 [10.1038/s41588-019-0381-3](https://doi.org/10.1038/s41588-019-0381-3) |

## Cache

Downloaded data sets are cached under `marker_cache_dir()`, which
resolves (in order) `getOption("markerdata.cache")`, the environment
variable `MARKERDATA_CACHE`, or else
`tools::R_user_dir("markerdata", "cache")`. Use `marker_download()` to
fetch one or more data sets ahead of time (for example before working
offline), and `marker_cache_clear()` to remove cached files again.

## Citing

Cite the package itself with `citation("markerdata")`. Every data set
keeps its own original license and must be cited on its own terms; get
its citation with `marker_references("<name>")`, e.g.
`marker_references("Rice_div")`.

## License

The markerdata package code is released under the MIT license. Each
bundled or downloadable data set keeps its own original license, listed
in the License column of the table above.
