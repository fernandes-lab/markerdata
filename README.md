
# markerdata: Marker Data Sets for Breeding Simulation

[![R-CMD-check](https://github.com/fernandes-lab/markerdata/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fernandes-lab/markerdata/actions/workflows/R-CMD-check.yaml)

markerdata provides curated single-nucleotide polymorphism (SNP) marker
panels drawn from public genomic resources across crops, trees,
livestock, fish and two sample-only reference species (the table below
lists every data set and species). Each panel is a plain data frame with
five marker-metadata columns (`snp`, `allele`, `chr`, `pos`, `cm`)
followed by one `-1`/`0`/`1` genotype column per individual, so it can
be used directly as founder genotypes for breeding-program simulation
with `simplePHENOTYPES::as_population()` or with breedingDesigner.

Every data set works offline. One small maize panel is bundled whole;
each of the other data sets ships a sample (2,000 SNPs, every
individual; 1,000 SNPs for the four largest) that loads instantly. For
all but two of them, `marker_data(name, full = TRUE)` downloads the
whole panel on first use and caches it locally, and `marker_download()`
fetches every hosted whole panel at once (about 19 MB).
`Arabidopsis_div` (1001 Genomes) and `Human_div` (1000 Genomes Project)
are **sample only**: markerdata does not host their whole data sets, and
`marker_catalog()$source_url` gives the original download.

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
#>                     name     common_name n_snp n_ind   het           mating
#> 1  SNP55K_maize282_maf04           maize 10650   280 0.004           inbred
#> 2          SoySNP50K_div         soybean 12000 20087 0.000           inbred
#> 3          EUChip60K_div      eucalyptus 12000   768 0.463          outbred
#> 4            Sorghum_div         sorghum 12000   842 0.051           inbred
#> 5               Rice_div            rice 12000   370 0.033           inbred
#> 6             Cattle_div          cattle 12000  1543 0.287          outbred
#> 7             Poplar_div          poplar 12000   434 0.318          outbred
#> 8              Wheat_div           wheat 12000  1854 0.076           inbred
#> 9        Arabidopsis_div     arabidopsis  2000  1097 0.000 inbred (selfing)
#> 10             Human_div           human  2000  2504 0.195          outbred
#> 11         SugarBeet_div      sugar beet 12000  1928 0.139          outbred
#> 12            Tomato_div          tomato 12000   353 0.077  inbred (mostly)
#> 13          WheatHex_div     bread wheat 12000   159 0.070           inbred
#> 14            Salmon_div Atlantic salmon 12000   777 0.298          outbred
#> 15             Sheep_div           sheep 12000   671 0.347          outbred
#> 16            Potato_div          potato  8659   586 0.690 outbred (clonal)
#> 17         Sugarcane_div       sugarcane  6615   432 0.121 outbred (clonal)
#>    available
#> 1       TRUE
#> 2      FALSE
#> 3      FALSE
#> 4       TRUE
#> 5      FALSE
#> 6      FALSE
#> 7      FALSE
#> 8      FALSE
#> 9      FALSE
#> 10     FALSE
#> 11     FALSE
#> 12     FALSE
#> 13     FALSE
#> 14     FALSE
#> 15     FALSE
#> 16     FALSE
#> 17     FALSE
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

Every other panel returns its bundled sample by default, also without
any download:

``` r
dim(marker_data("Cattle_div"))
#> 'Cattle_div': returning the bundled sample (1000 SNPs x 1543 individuals of 12000 x 1543). Use full = TRUE for the whole panel.
#> [1] 1000 1548
```

Ask for the whole panel with `full = TRUE`; it is downloaded on first
use and cached afterward. To download every whole panel in one go, call
`marker_download()`. The `-1`/`0`/`1` genotype coding is exactly what
`simplePHENOTYPES::as_population()` expects for founder genotypes:

``` r
geno <- marker_data("Cattle_div", full = TRUE)
pop <- simplePHENOTYPES::as_population(geno, individuals = 1:100)

marker_download() # every hosted whole panel
```

## Data sets

| Data set | Species | Panel | SNPs x individuals | Heterozygous calls | Mating | License | Full download | Reference |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| SNP55K_maize282_maf04 | Zea mays | 282 maize inbred association panel, Illumina MaizeSNP50 (55K) array, MAF\>0.4 subset | 10,650 x 280 | 0.4% | inbred | Panzea public data | bundled | Cook et al. (2012) Plant Physiology 158(2):824-834 [10.1104/pp.111.185033](https://doi.org/10.1104/pp.111.185033) |
| SoySNP50K_div | Glycine max | USDA Soybean Germplasm Collection diversity subset, SoySNP50K iSelect BeadChip | 12,000 x 20,087 | 0.0% | inbred | Public USDA/SoyBase data (cite) | 6.2 MB | Song et al. (2013) PLoS ONE 8(1):e54985 [10.1371/journal.pone.0054985](https://doi.org/10.1371/journal.pone.0054985); Song et al. (2015) G3 5(10):1999-2006 [10.1534/g3.115.019000](https://doi.org/10.1534/g3.115.019000) |
| EUChip60K_div | Eucalyptus grandis x E. urophylla | synthetic breeding population, EUChip60K Illumina Infinium array | 12,000 x 768 | 46.3% | outbred | CC0-1.0 | 0.5 MB | Resende et al. (2017) Heredity 119(4):245-255 [10.1038/hdy.2017.37](https://doi.org/10.1038/hdy.2017.37); Silva-Junior et al. (2015) New Phytologist 206(4):1527-1540 [10.1111/nph.13322](https://doi.org/10.1111/nph.13322) |
| Sorghum_div | Sorghum bicolor | WEST Diversity Panel, genotyping-by-sequencing (GBS) | 12,000 x 842 | 5.1% | inbred | CC-BY-4.0 | 0.9 MB | Ferguson et al. (2021) Plant Physiology 187(3):1481-1500 [10.1093/plphys/kiab346](https://doi.org/10.1093/plphys/kiab346) |
| Rice_div | Oryza sativa | IRRI elite tropical rice (indica / indica-admixed) genomic-selection panel, GBS | 12,000 x 370 | 3.3% | inbred | CC0-1.0 | 0.3 MB | Spindel et al. (2015) PLoS Genetics 11(2):e1004982 [10.1371/journal.pgen.1004982](https://doi.org/10.1371/journal.pgen.1004982) |
| Cattle_div | Bos taurus | worldwide cattle diversity panel (134 breeds), Illumina BovineSNP50 array | 12,000 x 1,543 | 28.7% | outbred | CC0-1.0 | 3.0 MB | Decker et al. (2014) PLoS Genetics 10(3):e1004254 [10.1371/journal.pgen.1004254](https://doi.org/10.1371/journal.pgen.1004254) |
| Poplar_div | Populus trichocarpa | landscape-genomics panel, PoplarConsortium-ORNL 34K SNP array | 12,000 x 434 | 31.8% | outbred | CC0-1.0 | 1.0 MB | Geraldes et al. (2014) Evolution 68(11):3260-3280 [10.1111/evo.12497](https://doi.org/10.1111/evo.12497); Geraldes et al. (2013) Molecular Ecology Resources 13(2):306-323 [10.1111/1755-0998.12056](https://doi.org/10.1111/1755-0998.12056) |
| Wheat_div | Triticum turgidum | Tetraploid wheat Global Collection (TGC), 90K iSelect SNP array | 12,000 x 1,854 | 7.6% | inbred | T3/Wheat public data (cite; terms to confirm) | 1.6 MB | Maccaferri et al. (2019) Nature Genetics 51:885-895 [10.1038/s41588-019-0381-3](https://doi.org/10.1038/s41588-019-0381-3) |
| Arabidopsis_div | Arabidopsis thaliana | 1001 Genomes, 1,097 accessions, 6,000-SNP pilot panel, TAIR10; SAMPLE ONLY (2,000 SNPs x 1,097 accessions): the whole data set is at source_url | 2,000 x 1,097 | 0.0% | inbred (selfing) | 1001 Genomes Data Usage Policy (cite; acknowledge the Weigel lab, MPI; redistribution terms to confirm) | sample only; [original data](https://1001genomes.org/data/GMI-MPI/releases/v3.1/) | 1001 Genomes Consortium (2016) Cell 166(2):481-491 [10.1016/j.cell.2016.05.063](https://doi.org/10.1016/j.cell.2016.05.063) |
| Human_div | Homo sapiens | 1000 Genomes Project phase 3, 2,504 people, 6,000-SNP pilot panel, GRCh37; SAMPLE ONLY (2,000 SNPs x 2,504 people): the whole data set is at source_url | 2,000 x 2,504 | 19.5% | outbred | IGSR open data, no embargo (cite; see the IGSR data reuse statement) | sample only; [original data](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) | 1000 Genomes Project Consortium (2015) Nature 526:68-74 [10.1038/nature15393](https://doi.org/10.1038/nature15393) |
| SugarBeet_div | Beta vulgaris | 1,928 public germplasm accessions (cultivated and wild sea beet), genotyping-by-sequencing (GBS), imputed by the source | 12,000 x 1,928 | 13.9% | outbred | CC-BY-4.0 | 2.9 MB | Tehseen et al. (2023) Crop Science 63(4):2255-2273 [10.1002/csc2.20978](https://doi.org/10.1002/csc2.20978) |
| Tomato_div | Solanum lycopersicum | two germplasm collections (162 + 191 accessions), Axiom 51K tomato array; allele column is major/minor (the source has no reference allele) | 12,000 x 353 | 7.7% | inbred (mostly) | CC-BY-4.0 | 0.3 MB | Yeon et al. (2024) BMC Plant Biology 24:222 [10.1186/s12870-024-04934-8](https://doi.org/10.1186/s12870-024-04934-8) |
| WheatHex_div | Triticum aestivum | 159 hexaploid wheat lines, Axiom 55K-type array | 12,000 x 159 | 7.0% | inbred | CC0-1.0 | 0.3 MB | Zhang & Gao (2024) Dryad (data-only deposit) [10.5061/dryad.4xgxd25ht](https://doi.org/10.5061/dryad.4xgxd25ht) |
| Salmon_div | Salmo salar | farmed Tasmanian breeding population, Affymetrix Axiom array | 12,000 x 777 | 29.8% | outbred | CC0-1.0 | 0.8 MB | Kijas et al. (2017) Animal Genetics 48(2):237-241 [10.1111/age.12513](https://doi.org/10.1111/age.12513) |
| Sheep_div | Ovis aries | 37 populations (Merino and Merino-derived breeds, other domestic breeds, a few wild Ovis), Illumina OvineSNP50 | 12,000 x 671 | 34.7% | outbred | CC0-1.0 | 1.5 MB | Ciani et al. (2015) Genetics Selection Evolution 47:64 [10.1186/s12711-015-0139-z](https://doi.org/10.1186/s12711-015-0139-z) |
| Potato_div | Solanum tuberosum | 586 Wisconsin breeding clones, SNP array, DMv4.03 positions; DIPLOIDIZED (dosage 0 -\> -1, dosage 4 -\> 1, dosages 1-3 -\> 0); allele column is N/N | 8,659 x 586 | 69.0% | outbred (clonal) | CC0-1.0 | 0.4 MB | Caraza-Harter & Endelman (2022) Theoretical and Applied Genetics 135(9):2943-2951 [10.1007/s00122-022-04159-z](https://doi.org/10.1007/s00122-022-04159-z) |
| Sugarcane_div | Saccharum spp. hybrids | 432 clones (USDA-ARS Canal Point breeding clones and commercial cultivars), capture sequencing; positions are Sorghum bicolor v3.1 coordinates | 6,615 x 432 | 12.1% | outbred (clonal) | Public domain (US Government work) | 0.4 MB | Islam et al. (2023) Frontiers in Plant Science 14:1205999 [10.3389/fpls.2023.1205999](https://doi.org/10.3389/fpls.2023.1205999) |

Sizes are those of the whole panels, except for `Arabidopsis_div` and
`Human_div`, where the counts describe the sample because the whole data
set is not hosted (follow the link in the last-but-one column). The
human data are from publicly released, consented 1000 Genomes samples
and are included for method development and teaching. The bundled sample
of each panel is 2,000 SNPs (1,000 for `SoySNP50K_div`, `Cattle_div`,
`Wheat_div` and `SugarBeet_div`, to keep the package small) and every
individual of the panel: SNPs evenly spaced over those polymorphic among
the individuals.

Two data sets need care. `Potato_div` is an autotetraploid stored
**diploidized** (dosage 0 becomes `-1`, dosage 4 becomes `1`, dosages 1
to 3 become `0`), so its allele dosage is lost; `Sugarcane_div` holds
the diploid-style calls of its source (no dosage exists) at *Sorghum
bicolor* positions. Their help topics explain what is lost. Whole panels
are hosted in the GitHub releases `data-v1` and `data-v2` (the `release`
column of `marker_catalog()` says which).

## Cache

Downloaded whole panels are cached under `marker_cache_dir()`, which
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
