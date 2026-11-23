# Documentation for markerdata's marker data sets.
#
# SNP55K_maize282_maf04 is the only object actually bundled with the package
# (data/SNP55K_maize282_maf04.rda), so it gets a standard data-documentation
# block ending in the string "SNP55K_maize282_maf04". Every other data set
# below is a help TOPIC only (each block ends in NULL and carries an explicit
# @name): they are not package objects. marker_data("<name>") returns the
# sample (2,000 SNPs, every individual) shipped in inst/extdata/samples/
# (built by data-raw/samples.R); marker_data("<name>", full = TRUE) downloads the whole
# panel from the package's GitHub release and caches it locally (see
# marker_cache_dir() / marker_cache_clear()). Arabidopsis_div and Human_div
# are sample-only: no whole panel is hosted, and the catalog's source_url
# points to the original download (built by data-raw/s1database_panels.R).
#
# data-raw/catalog.R is the single source of truth for each dataset's
# metadata (species, panel, size, license, source, literature reference) and
# is also exposed at runtime via marker_catalog() and marker_references().

#' Maize 282 association panel genotypes (MaizeSNP50 subset)
#'
#' A maize marker matrix for the 282 inbred association panel, genotyped on
#' the Illumina MaizeSNP50 (55K) array and filtered to SNPs with
#' minor-allele frequency (MAF) greater than 0.4 to keep the bundled file
#' small: 10,650 SNPs by 280 inbred lines. Genotypes are coded `-1`/`0`/`1`
#' (reference homozygote / heterozygote / alternate homozygote). This is
#' the only dataset actually bundled whole with markerdata; every other
#' data set in the package (see [marker_catalog()]) ships a sample, and all
#' but [Arabidopsis_div] and [Human_div], which are samples only, can be
#' downloaded whole on request with [marker_data()].
#'
#' @format A data frame with 10,650 rows (SNPs) and 285 columns. The first
#'   five columns are marker metadata (`snp`, `allele`, `chr`, `pos`, `cm`);
#'   the remaining 280 columns are inbred lines, each a numeric (double)
#'   genotype in `{-1, 0, 1}`.
#' @section Columns:
#' \describe{
#'   \item{snp}{marker name.}
#'   \item{allele}{allele pair, e.g. "A/G".}
#'   \item{chr}{chromosome, 1-10.}
#'   \item{pos}{physical position in base pairs (B73 AGPv2).}
#'   \item{cm}{genetic position in centiMorgans (see the section below).}
#' }
#' The remaining columns, one per inbred line and named by accession, each
#' hold a numeric (double) genotype in `{-1, 0, 1}` (reference homozygote /
#' heterozygote / alternate homozygote).
#' @section The genetic map is synthetic:
#' `cm` holds a modelled genetic map, generated from the physical positions
#' with `simplePHENOTYPES::synthetic_map()`; it is not taken from any
#' published maize linkage map (IBM, NAM, or otherwise) and should not be
#' quoted as a measured recombination distance. It exists so that
#' meiosis-based simulation functions (crossing, selfing, doubled haploids)
#' have usable recombination distances on this bundled example data.
#' @source Panzea, maize282set genotypes (Illumina MaizeSNP50 array, MAF >
#'   0.4 subset).
#' @references
#' Cook JP, McMullen MD, Holland JB, Tian F, Bradbury P, Ross-Ibarra J,
#' Buckler ES, Flint-Garcia SA (2012). Genetic architecture of maize kernel
#' composition in the nested association mapping and inbred association
#' panels. *Plant Physiology* 158(2):824-834. \doi{10.1104/pp.111.185033}
#' @keywords datasets
#' @family marker datasets
#' @examples
#' head(SNP55K_maize282_maf04[, 1:8])
#' identical(marker_data("SNP55K_maize282_maf04"), SNP55K_maize282_maf04)
"SNP55K_maize282_maf04"

#' Soybean SoySNP50K genotypes (SoyBase diversity subset)
#'
#' @description
#' A soybean marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from the SoySNP50K genotyping of the
#' USDA Soybean Germplasm Collection: all 20,087 accessions and 12,000 SNPs
#' (MAF >= 0.05 and at most 5% missing calls; remaining missing calls
#' imputed to the major homozygote) on the Wm82.a2 (`Wm82.gnm2`) assembly. The
#' package ships a sample of 1,000 SNPs (fewer than the other samples, to keep
#' the package small) and every individual, returned by
#' `marker_data("SoySNP50K_div")`; the whole panel described here is
#' downloaded on first use with `marker_data("SoySNP50K_div", full = TRUE)`
#' (a packed, xz-compressed download of about 6 MB) and cached locally
#' afterward (see [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1` (reference homozygote / heterozygote /
#' alternate homozygote). Like the inbred [Sorghum_div] accessions (and
#' unlike the outbred [EUChip60K_div] trees), these accessions are largely
#' homozygous, so the dataset suits additive schemes (line development,
#' genomic selection) rather than heterosis/combining-ability schemes.
#'
#' The SoyBase HapMap spans not just the 20 soybean chromosomes but also a
#' handful of unplaced scaffold contigs: 68 of the 12,000 SNPs in this
#' subset fall on scaffolds, so the `chr` column actually has 35 distinct
#' labels (chromosomes Gm01 through Gm20, plus 15 `scaffold_*` contigs),
#' not 20.
#'
#' The genetic-map column `cm` is a synthetic map built from physical
#' positions with `simplePHENOTYPES::synthetic_map()` (the SoyBase files
#' carry no cM); replace it with a real map if you have one.
#'
#' @format A data frame of 12,000 rows (SNPs) by 20,092 columns: the first
#'   five are `snp`, `allele`, `chr`, `pos`, `cm`; the remaining 20,087 are
#'   accessions (USDA PI numbers), each an integer genotype in
#'   `{-1, 0, 1}`.
#' @source SoyBase,
#'   <https://data.soybase.org/Glycine/max/diversity/Wm82.gnm2.div.Song_Hyten_2015/>
#'   (HapMap file `glyma.Wm82.gnm2.div.Song_Hyten_2015.hmp.gz`). Public
#'   USDA/SoyBase data; please cite the references below. Built by
#'   `data-raw/soybase_soysnp50k.R`.
#' @references
#' Song Q, Hyten DL, Jia G, Quigley CV, Fickus EW, Nelson RL, Cregan PB
#' (2013). Development and evaluation of SoySNP50K, a high-density
#' genotyping array for soybean. *PLoS ONE* 8(1):e54985.
#' \doi{10.1371/journal.pone.0054985}
#'
#' Song Q, Hyten DL, Jia G, Quigley CV, Fickus EW, Nelson RL, Cregan PB
#' (2015). Fingerprinting soybean germplasm and its utility in genomic
#' research. *G3: Genes, Genomes, Genetics* 5(10):1999-2006.
#' \doi{10.1534/g3.115.019000}
#' @name SoySNP50K_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("SoySNP50K_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("SoySNP50K_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Eucalyptus EUChip60K genotypes (grandis x urophylla diversity subset)
#'
#' @description
#' A eucalyptus marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from an *Eucalyptus grandis* x
#' *E. urophylla* synthetic breeding population genotyped on the EUChip60K
#' Illumina Infinium array: 768 trees and 12,000 SNPs (MAF >= 0.05 in the
#' subset). The package ships a
#' sample of 2,000 SNPs and every individual, returned by
#' `marker_data("EUChip60K_div")`; the whole panel described here is downloaded
#' on first use with `marker_data("EUChip60K_div", full = TRUE)` (a packed,
#' xz-compressed download of about 2 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1` (reference homozygote / heterozygote /
#' alternate homozygote). Unlike the inbred [SoySNP50K_div] accessions,
#' these trees are **outbred**: about 46% of genotype calls are
#' heterozygous, so a dominance layer can be simulated on the population.
#' This makes the dataset suitable for schemes that exploit heterosis /
#' combining ability, such as reciprocal recurrent selection.
#'
#' The genetic-map column `cm` is a synthetic map built from physical
#' positions with `simplePHENOTYPES::synthetic_map()` (the array carries no
#' cM); replace it with a real map if you have one.
#'
#' @format A data frame of 12,000 rows (SNPs) by 773 columns: the first
#'   five are `snp`, `allele`, `chr`, `pos`, `cm`; the remaining 768 are
#'   trees, each an integer genotype in `{-1, 0, 1}`. SNP IDs and
#'   chromosome labels follow the array's `EuBR<chr>s<pos>` convention.
#' @source Dryad Digital Repository, \doi{10.5061/dryad.ms580} (license
#'   CC0-1.0, public domain), file `SNP_Data.txt`. Built by
#'   `data-raw/dryad_euchip60k.R`.
#' @references
#' Resende RT, Resende MDV, Silva FF, Azevedo CF, Takahashi EK,
#' Silva-Junior OB, Grattapaglia D (2017). Assessing the expected response
#' to genomic selection of individuals and families in Eucalyptus breeding
#' with an additive-dominant model. *Heredity* 119(4):245-255.
#' \doi{10.1038/hdy.2017.37}
#'
#' Silva-Junior OB, Faria DA, Grattapaglia D (2015). A flexible
#' multispecies genome-wide 60K SNP chip developed from pooled
#' resequencing of 240 Eucalyptus tree genomes across 12 species.
#' *New Phytologist* 206(4):1527-1540. \doi{10.1111/nph.13322}
#' @name EUChip60K_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("EUChip60K_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("EUChip60K_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Sorghum WEST Diversity Panel genotypes (GBS subset)
#'
#' @description
#' A sorghum marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from the sorghum **WEST Diversity
#' Panel** genotyped by genotyping-by-sequencing (GBS) for the
#' water-use-efficiency GWAS/TWAS study of Ferguson et al. (2021): 842
#' accessions and 12,000 SNPs (MAF >= 0.05 in the subset). The package ships a
#' sample of 2,000 SNPs and every individual, returned by
#' `marker_data("Sorghum_div")`; the whole panel described here is downloaded
#' on first use with `marker_data("Sorghum_div", full = TRUE)` (a packed,
#' xz-compressed download of about 3 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1` (reference homozygote / heterozygote /
#' alternate homozygote). Like the inbred [SoySNP50K_div] accessions (and
#' unlike the outbred [EUChip60K_div] trees), these accessions are largely
#' homozygous: about 5% of genotype calls are heterozygous, so the dataset
#' suits additive schemes (line development, genomic selection) rather than
#' heterosis/combining-ability schemes.
#'
#' The genetic-map column `cm` is a synthetic map built from physical
#' positions with `simplePHENOTYPES::synthetic_map()`; `snp` IDs are
#' synthesized as `S<chr>_<pos>` (the source matrix carries none). Replace
#' the map with a real one if you have it.
#'
#' @format A data frame of 12,000 rows (SNPs) by 847 columns: the first
#'   five are `snp`, `allele`, `chr`, `pos`, `cm`; the remaining 842 are
#'   accessions, each an integer genotype in `{-1, 0, 1}`.
#' @source WEST Diversity Panel GBS, Zenodo \doi{10.5281/zenodo.5019227}
#'   (license CC-BY-4.0). Built by `data-raw/sorghum_wue.R` from a recoded
#'   `-1`/`0`/`1` genotype matrix derived from that deposit.
#' @references
#' Ferguson JN, Fernandes SB, Monier B, Miller ND, Allen D, Dmitrieva A,
#' Schmuker P, Lozano R, Valluru R, Buckler ES, Gore MA, Brown PJ,
#' Spalding EP, Leakey ADB (2021). Machine learning-enabled phenotyping for
#' GWAS and TWAS of WUE traits in 869 field-grown sorghum accessions.
#' *Plant Physiology* 187(3):1481-1500. \doi{10.1093/plphys/kiab346}
#'
#' Ferguson JN, Fernandes SB, Monier B, et al. (2020). WEST Diversity Panel
#' GBS information. *Zenodo*, \doi{10.5281/zenodo.5019227} (CC-BY-4.0).
#' @name Sorghum_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("Sorghum_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("Sorghum_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Rice IRRI tropical genomic-selection panel genotypes (GBS subset)
#'
#' @description
#' A rice marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from the IRRI elite tropical rice
#' (indica / indica-admixed) genomic-selection breeding panel genotyped by
#' GBS (Spindel et al. 2015): 370 lines and 12,000 SNPs (MAF >= 0.05 in the
#' subset). The package ships a
#' sample of 2,000 SNPs and every individual, returned by
#' `marker_data("Rice_div")`; the whole panel described here is downloaded
#' on first use with `marker_data("Rice_div", full = TRUE)` (a packed,
#' xz-compressed download of about 1 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1`. These are inbred lines (about 3%
#' heterozygous), so the dataset suits additive schemes (line development,
#' genomic selection).
#'
#' The genetic-map column `cm` is a synthetic map built from physical
#' positions with `simplePHENOTYPES::synthetic_map()`.
#'
#' @format A data frame of 12,000 rows (SNPs) by 375 columns (5 map columns
#'   + 370 lines), each genotype in `{-1, 0, 1}`.
#' @source Dryad \doi{10.5061/dryad.7369p} (CC0), file
#'   `MET_crfilt_.75_allchrom.hmp.txt`; built by
#'   `data-raw/dryad_rice_spindel.R`.
#' @references
#' Spindel J, Begum H, Akdemir D, Virk P, Collard B, Redona E, Atlin G,
#' Jannink J-L, McCouch SR (2015). Genomic Selection and Association
#' Mapping in Rice (*Oryza sativa*). *PLoS Genetics* 11(2):e1004982.
#' \doi{10.1371/journal.pgen.1004982}
#' @name Rice_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("Rice_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("Rice_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Cattle worldwide diversity panel genotypes (BovineSNP50 subset)
#'
#' @description
#' A cattle marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from a worldwide cattle diversity panel
#' (134 breeds) genotyped on the Illumina BovineSNP50 array (Decker et al.
#' 2014): 1,543 animals and 12,000 SNPs (MAF >= 0.05 in the subset, UMD3.1
#' positions). The package ships a
#' sample of 1,000 SNPs and every individual, returned by
#' `marker_data("Cattle_div")`; the whole panel described here is downloaded
#' on first use with `marker_data("Cattle_div", full = TRUE)` (a packed,
#' xz-compressed download of about 5 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1`. Cattle are **outbred** (about 29%
#' heterozygous), so a dominance layer can be simulated (cf.
#' [EUChip60K_div]) --- useful for animal-breeding schemes exploiting
#' heterosis / combining ability.
#'
#' The genetic-map column `cm` is a synthetic map built from physical
#' positions with `simplePHENOTYPES::synthetic_map()`.
#'
#' @format A data frame of 12,000 rows (SNPs) by 1,548 columns (5 map
#'   columns + 1,543 animals), each genotype in `{-1, 0, 1}`.
#' @source Dryad \doi{10.5061/dryad.th092} (CC0), PLINK
#'   `ww_cattle_dryad_20Sept2013.ped/.map`; built by
#'   `data-raw/dryad_cattle_decker.R`.
#' @references
#' Decker JE, McKay SD, Rolf MM, Kim J, Molina Alcala A, Sonstegard TS,
#' Hanotte O, et al. (2014). Worldwide Patterns of Ancestry, Divergence,
#' and Admixture in Domesticated Cattle. *PLoS Genetics* 10(3):e1004254.
#' \doi{10.1371/journal.pgen.1004254}
#' @name Cattle_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("Cattle_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("Cattle_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Poplar (Populus trichocarpa) landscape-genomics panel genotypes (34K subset)
#'
#' @description
#' A poplar marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from a *Populus trichocarpa*
#' landscape-genomics panel genotyped on the PoplarConsortium-ORNL 34K SNP
#' array (Geraldes et al. 2014): 434 accessions and 12,000 SNPs (MAF >=
#' 0.05 in the subset). The package ships a
#' sample of 2,000 SNPs and every individual, returned by
#' `marker_data("Poplar_div")`; the whole panel described here is downloaded
#' on first use with `marker_data("Poplar_div", full = TRUE)` (a packed,
#' xz-compressed download of about 1 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1`. These forest trees are **outbred**
#' (about 32% heterozygous), so a dominance layer can be simulated (cf.
#' [EUChip60K_div]) --- useful for tree-breeding schemes.
#'
#' The genetic-map column `cm` is a synthetic map built from physical
#' positions with `simplePHENOTYPES::synthetic_map()`.
#'
#' @format A data frame of 12,000 rows (SNPs) by 439 columns (5 map columns
#'   + 434 accessions), each genotype in `{-1, 0, 1}`. SNP IDs are the
#'   array's `scaffold_<chr>_<pos>` names.
#' @source Dryad \doi{10.5061/dryad.7s848} (CC0), file
#'   `FOR_DRYAD_TandB_fwd_june2013_434.txt`; built by
#'   `data-raw/dryad_poplar_geraldes.R`.
#' @references
#' Geraldes A, Farzaneh N, Grassa CJ, McKown AD, Guy RD, Mansfield SD,
#' Douglas CJ, Cronk QCB (2014). Landscape genomics of *Populus
#' trichocarpa*. *Evolution* 68(11):3260-3280. \doi{10.1111/evo.12497}
#'
#' Geraldes A, DiFazio SP, Slavov GT, et al. (2013). A 34K SNP genotyping
#' array for *Populus trichocarpa*. *Molecular Ecology Resources*
#' 13(2):306-323. \doi{10.1111/1755-0998.12056}
#' @name Poplar_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("Poplar_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("Poplar_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Wheat Tetraploid Global Collection genotypes (90K SNP subset)
#'
#' @description
#' A wheat marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from the **Tetraploid wheat Global
#' Collection (TGC)** --- wild and domesticated emmer, durum wheat
#' landraces and tetraploid subspecies --- genotyped on the wheat 90K
#' iSelect SNP array (Maccaferri et al. 2019): 1,854 accessions and 12,000
#' SNPs (MAF >= 0.05 in the subset, spanning the 14 nuclear chromosomes of
#' the A and B subgenomes, 1A--7B). The package ships a
#' sample of 1,000 SNPs and every individual, returned by
#' `marker_data("Wheat_div")`; the whole panel described here is downloaded
#' on first use with `marker_data("Wheat_div", full = TRUE)` (a packed,
#' xz-compressed download of about 6 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1`. Tetraploid wheat lines are largely
#' homozygous (about 8% heterozygous), so the dataset suits additive
#' schemes (cf. [SoySNP50K_div]).
#'
#' **Polyploidy note.** Durum/emmer wheat is an *allotetraploid*
#' (2n = 4x = 28, AABB) with **disomic** inheritance: homologues pair
#' within a subgenome, and homoeologues do not (under the *Ph1* system). A
#' subgenome-specific SNP on chromosome 1A therefore segregates as an
#' ordinary diploid locus linked only to other 1A markers, which is exactly
#' what coding each subgenome-specific 90K SNP as a diploid locus (`-1`/
#' `0`/`1`, from the deposited VCF's `0/0`, `0/1`, `1/1` calls) with one
#' genetic map per chromosome 1A..7B represents: the diploid-per-subgenome
#' model captures the inheritance correctly. The only approximation is at
#' the assay level: the 90K array's probes are subgenome-specific, so a
#' marker that happened to co-amplify both homoeologues would need a
#' composite 0..4 dosage call, which is a data-curation matter rather than
#' a meiosis-model limitation.
#'
#' The genetic-map column `cm` is a synthetic map built from physical
#' positions with `simplePHENOTYPES::synthetic_map()`, one map per
#' chromosome.
#'
#' @format A data frame of 12,000 rows (SNPs) by 1,859 columns (5 map
#'   columns + 1,854 accessions), each genotype in `{-1, 0, 1}`. SNP IDs
#'   are the 90K array's `IWB` probe names (made unique across
#'   subgenomes).
#' @source T3/Wheat (Triticeae Toolbox), file
#'   `..._TGC-GTC_collection_SNP90K_mapped_filtered.vcf`; built by
#'   `data-raw/t3_wheat_tgc.R`. Reference per USDA GrainGenes ("Global
#'   Durum Genomic Resources").
#' @references
#' Maccaferri M, Harris NS, Twardziok SO, Pasam RK, Gundlach H, Spannagl M,
#' et al., Cattivelli L, Tuberosa R (2019). Durum wheat genome highlights
#' past domestication signatures and future improvement targets. *Nature
#' Genetics* 51:885-895. \doi{10.1038/s41588-019-0381-3}
#' @name Wheat_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("Wheat_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("Wheat_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Sugar beet public germplasm genotypes (GBS subset)
#'
#' @description
#' A sugar beet (*Beta vulgaris*) marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from the genotyping-by-sequencing (GBS)
#' data of 1,928 publicly available US *Beta vulgaris* germplasm accessions
#' (cultivated beet and wild sea beet; Tehseen et al. 2023): 12,000 SNPs
#' (MAF >= 0.05 in the subset). The package ships a sample of 1,000 SNPs and every individual, returned by `marker_data("SugarBeet_div")`; the whole
#' panel described here is downloaded on first use with
#' `marker_data("SugarBeet_div", full = TRUE)` (a packed, xz-compressed
#' download of about 3 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1`. Sugar beet is outcrossing, so the
#' accessions are not inbred: about 14% of genotype calls are heterozygous.
#' The source file was already imputed by its authors; calls still missing
#' were set to the rounded per-SNP mean.
#'
#' The `chr` column has 17 labels: chromosomes `1` to `9` plus 8 unplaced
#' scaffold labels (`SFD_*`), which are kept as in the source. The genome
#' assembly is not stated by the source. The genetic-map column `cm` is a
#' synthetic map built from physical positions with
#' `simplePHENOTYPES::synthetic_map()`.
#'
#' @format A data frame of 12,000 rows (SNPs) by 1,933 columns (5 map
#'   columns + 1,928 accessions), each genotype in `{-1, 0, 1}`.
#' @source Zenodo, <https://zenodo.org/records/11551725> (license
#'   CC-BY-4.0), file `GBS_imputed_MA05_100bp_082421.hmp.txt`; built by
#'   `data-raw/zenodo_sugarbeet_chu.R`.
#' @references
#' Tehseen MM, Poore RC, Fugate KK, Bolton MD, Ramachandran V, Wyatt NA, Li
#' X, Chu C (2023). Potential of publicly available *Beta vulgaris*
#' germplasm for sustainable sugarbeet improvement indicated by combining
#' analysis of genetic diversity and historic resistance evaluation. *Crop
#' Science* 63(4):2255-2273. \doi{10.1002/csc2.20978}
#' @name SugarBeet_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("SugarBeet_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("SugarBeet_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Tomato 51K array genotypes (two germplasm collections)
#'
#' @description
#' A tomato (*Solanum lycopersicum*) marker matrix in the same numeric
#' format as [SNP55K_maize282_maf04], curated from the Axiom 51K tomato
#' array genotyping of two germplasm collections (162 and 191 accessions,
#' 353 in all; Yeon et al. 2024): 12,000 SNPs (MAF >= 0.05 in the subset).
#' The package ships a sample of 2,000 SNPs and every individual, returned by
#' `marker_data("Tomato_div")`; the whole panel described here is downloaded
#' on first use with `marker_data("Tomato_div", full = TRUE)` (a packed,
#' xz-compressed download of about 0.3 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1`. Cultivated tomato is a selfing
#' species, so the accessions are mostly homozygous (about 8% heterozygous
#' calls) and the data set suits additive schemes (cf. [SoySNP50K_div]).
#'
#' **The `allele` column is `major/minor`.** The source has no reference
#' allele, so for each SNP the allele seen more often across the 353
#' accessions is coded as the first allele (`-1` when homozygous) and the
#' other as the second (`1`). The genome assembly is not stated by the
#' source. The genetic-map column `cm` is a synthetic map built from
#' physical positions with `simplePHENOTYPES::synthetic_map()`.
#'
#' @format A data frame of 12,000 rows (SNPs) by 358 columns (5 map columns
#'   + 353 accessions), each genotype in `{-1, 0, 1}`; chromosomes `01` to
#'   `12`.
#' @source Mendeley Data, <https://data.mendeley.com/datasets/bxcpc274fh/1>
#'   (\doi{10.17632/bxcpc274fh.1}; license CC-BY-4.0), files
#'   `TGC1 51K SNP chip genotypic data.xlsx` and `TGC2 51K SNP chip
#'   genotypic data.xlsx`; built by `data-raw/mendeley_tomato_51k.R`.
#' @references
#' Yeon J, Nguyen TTP, Kim M, Sim S-C (2024). Prediction accuracy of
#' genomic estimated breeding values for fruit traits in cultivated tomato
#' (*Solanum lycopersicum* L.). *BMC Plant Biology* 24:222.
#' \doi{10.1186/s12870-024-04934-8}
#' @name Tomato_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("Tomato_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("Tomato_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Hexaploid bread wheat genotypes (159 lines, 55K-type array)
#'
#' @description
#' A bread wheat (*Triticum aestivum*) marker matrix in the same numeric
#' format as [SNP55K_maize282_maf04], curated from an Axiom array
#' genotyping of 159 hexaploid wheat lines (Zhang and Gao 2024): 12,000
#' SNPs (MAF >= 0.05 in the subset) on the 21 nuclear chromosomes 1A to 7D.
#' The package ships a sample of 2,000 SNPs and every individual,
#' returned by `marker_data("WheatHex_div")`; because this panel has only
#' 159 lines, its sample keeps all of them. The whole panel is downloaded on
#' first use with `marker_data("WheatHex_div", full = TRUE)` (a packed,
#' xz-compressed download of about 0.3 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1`. The lines are largely homozygous
#' (about 7% heterozygous calls), so the data set suits additive schemes
#' (cf. [Wheat_div], the tetraploid wheat panel).
#'
#' **Polyploidy note.** Bread wheat is an *allohexaploid* (2n = 6x = 42,
#' AABBDD) with **disomic** inheritance, as in [Wheat_div]: each
#' subgenome-specific SNP is coded as an ordinary diploid locus, with one
#' genetic map per chromosome 1A..7D.
#'
#' This is a data-only deposit: it has no journal article. The genetic-map
#' column `cm` is a synthetic map built from physical positions with
#' `simplePHENOTYPES::synthetic_map()`.
#'
#' @format A data frame of 12,000 rows (SNPs) by 164 columns (5 map columns
#'   + 159 lines), each genotype in `{-1, 0, 1}`.
#' @source Zenodo, <https://zenodo.org/records/10791194> (mirror of Dryad
#'   \doi{10.5061/dryad.4xgxd25ht}; license CC0-1.0), file
#'   `original_genotype_data.xlsb`; built by
#'   `data-raw/dryad_wheat_zhanggao.R`.
#' @references
#' Zhang X, Gao Z (2024). Original genotype data of 159 wheat samples.
#' *Dryad Digital Repository*. \doi{10.5061/dryad.4xgxd25ht}
#' @name WheatHex_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("WheatHex_div")) # bundled sample: all 159 lines
#' \donttest{
#' geno <- marker_data("WheatHex_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Atlantic salmon farmed breeding population genotypes (Axiom array)
#'
#' @description
#' An Atlantic salmon (*Salmo salar*) marker matrix in the same numeric
#' format as [SNP55K_maize282_maf04], curated from a farmed Tasmanian
#' breeding population of 777 fish genotyped on an Affymetrix Axiom array
#' (Kijas et al. 2017): 12,000 SNPs (MAF >= 0.05 in the subset) from the
#' 106,492 SNPs in the deposited files, on chromosomes 01 to 29. The package
#' ships a sample of 2,000 SNPs and every individual, returned by
#' `marker_data("Salmon_div")`; the whole panel described here is downloaded
#' on first use with `marker_data("Salmon_div", full = TRUE)` (a packed,
#' xz-compressed download of about 0.8 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1`. The population is **outbred**: about
#' 30% of genotype calls are heterozygous, so a dominance layer can be
#' simulated (cf. [Cattle_div]). Missing calls were set to the rounded
#' per-SNP mean. The genome assembly is not stated by the source. The
#' genetic-map column `cm` is a synthetic map built from physical positions
#' with `simplePHENOTYPES::synthetic_map()`.
#'
#' @format A data frame of 12,000 rows (SNPs) by 782 columns (5 map columns
#'   + 777 fish), each genotype in `{-1, 0, 1}`.
#' @source Zenodo, <https://zenodo.org/records/5013066> (mirror of Dryad
#'   \doi{10.5061/dryad.pm354}; license CC0-1.0), file `Salmon_SNP.zip`
#'   (PLINK); built by `data-raw/dryad_salmon_kijas.R`.
#' @references
#' Kijas J, Elliot N, Kube P, Evans B, Botwright N, King H, Primmer CR,
#' Verbyla K (2017). Diversity and linkage disequilibrium in farmed
#' Tasmanian Atlantic salmon. *Animal Genetics* 48(2):237-241.
#' \doi{10.1111/age.12513}
#' @name Salmon_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("Salmon_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("Salmon_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Sheep diversity panel genotypes (OvineSNP50, 37 populations)
#'
#' @description
#' A sheep (*Ovis aries*) marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from 671 animals of 37 populations
#' genotyped on the Illumina OvineSNP50 array (Ciani et al. 2015): 12,000
#' SNPs (MAF >= 0.05 in the subset, from 45,847 SNPs in the source) on
#' autosomes 01 to 26. The populations are Merino and Merino-derived
#' breeds, other domestic breeds, and a few wild *Ovis* (Argali, Urial,
#' Mouflon), which are kept. The package ships a sample of 2,000 SNPs and every individual, returned by `marker_data("Sheep_div")`; the whole panel
#' described here is downloaded on first use with
#' `marker_data("Sheep_div", full = TRUE)` (a packed, xz-compressed download
#' of about 1.5 MB) and cached locally afterward (see [marker_cache_dir()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1`. The animals are **outbred**: about 35%
#' of genotype calls are heterozygous, so a dominance layer can be
#' simulated (cf. [Cattle_div]). Because the panel mixes domestic breeds
#' and wild relatives, it has strong population structure. SNPs with no
#' chromosome (label 0) and the label 27 (assumed to be the X chromosome,
#' not confirmed) were dropped. Missing calls were set to the rounded
#' per-SNP mean. The genome assembly is not stated by the source. The
#' genetic-map column `cm` is a synthetic map built from physical positions
#' with `simplePHENOTYPES::synthetic_map()`.
#'
#' @format A data frame of 12,000 rows (SNPs) by 676 columns (5 map columns
#'   + 671 animals), each genotype in `{-1, 0, 1}`.
#' @source Dryad, <https://datadryad.org/dataset/doi:10.5061/dryad.2p0qf>
#'   (\doi{10.5061/dryad.2p0qf}; license CC0-1.0), file `CIANI.rar`
#'   (PLINK); built by `data-raw/dryad_sheep_ciani.R`.
#' @references
#' Ciani E, Lasagna E, D'Andrea M, Alloggio I, Marroni F, Ceccobelli S,
#' Delgado Bermejo JV, Sarti FM, Kijas J, Lenstra JA, Pilla F (2015).
#' Merino and Merino-derived sheep breeds: a genome-wide intercontinental
#' study. *Genetics Selection Evolution* 47:64.
#' \doi{10.1186/s12711-015-0139-z}
#' @name Sheep_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("Sheep_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("Sheep_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Potato breeding clone genotypes (diploidized tetraploid dosages)
#'
#' @description
#' A potato (*Solanum tuberosum*) marker matrix in the same numeric format
#' as [SNP55K_maize282_maf04], curated from the SNP array genotyping of 586
#' Wisconsin breeding clones (Caraza-Harter and Endelman 2022): 8,659 SNPs
#' (all those with MAF >= 0.05, computed from the original dosages), on
#' chromosomes 01 to 12 with DMv4.03 positions. The package ships a sample of 2,000 SNPs and every individual, returned by `marker_data("Potato_div")`;
#' the whole panel described here is downloaded on first use with
#' `marker_data("Potato_div", full = TRUE)` (a packed, xz-compressed
#' download of about 0.4 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @section Polyploid coding: what is lost:
#' Cultivated potato is an **autotetraploid**, and the source gives each
#' genotype as an allele dosage from 0 to 4. The `-1`/`0`/`1` format of this
#' package cannot hold five dosage classes, so the dosages were
#' **diploidized**: dosage 0 became `-1`, dosage 4 became `1`, and the
#' three intermediate dosages 1, 2 and 3 were all collapsed to `0`. About
#' 69% of the calls are therefore `0`, and `0` here means "any mixed
#' dosage", not a diploid heterozygote.
#'
#' What is lost: the dosage itself. Simplex (1), duplex (2) and triplex (3)
#' genotypes can no longer be told apart, so the data cannot be used to
#' simulate tetrasomic inheritance (double reduction, multivalent pairing)
#' or dosage-based additive and dominance effects. The allele frequency
#' used for the MAF filter was computed from the original dosages, but the
#' frequencies you compute from the stored `-1`/`0`/`1` calls will differ.
#' Use this data set to test diploid-style pipelines on a polyploid crop's
#' marker structure (map, linkage disequilibrium, population structure),
#' not to study autopolyploid genetics.
#'
#' @details
#' The `allele` column is `N/N` because the source gives dosages, not
#' nucleotides. The source has no missing calls. The genetic-map column
#' `cm` is a synthetic map built from physical positions with
#' `simplePHENOTYPES::synthetic_map()`.
#'
#' @format A data frame of 8,659 rows (SNPs) by 591 columns (5 map columns
#'   + 586 clones), each genotype in `{-1, 0, 1}` (see the section above).
#' @source Zenodo, <https://zenodo.org/records/6818259> (also Dryad
#'   \doi{10.5061/dryad.3n5tb2rk0}; license CC0-1.0), file
#'   `marker_data.csv`; built by `data-raw/dryad_potato_endelman.R`.
#' @references
#' Caraza-Harter MV, Endelman JB (2022). The genetic architectures of vine
#' and skin maturity in tetraploid potato. *Theoretical and Applied
#' Genetics* 135(9):2943-2951. \doi{10.1007/s00122-022-04159-z}
#' @name Potato_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("Potato_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("Potato_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Sugarcane breeding clone genotypes (diploid calls, sorghum positions)
#'
#' @description
#' A sugarcane (*Saccharum* spp. hybrids) marker matrix in the same numeric
#' format as [SNP55K_maize282_maf04], curated from the capture-sequencing
#' genotypes of 432 clones (414 USDA-ARS Canal Point breeding clones and 18
#' commercial cultivars; Islam et al. 2023): 6,615 SNPs (all those with
#' MAF >= 0.05) on chromosomes 01 to 10. The package ships a sample of 2,000 SNPs and every individual, returned by `marker_data("Sugarcane_div")`;
#' the whole panel described here is downloaded on first use with
#' `marker_data("Sugarcane_div", full = TRUE)` (a packed, xz-compressed
#' download of about 0.4 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
#'
#' @section Polyploid coding: what is lost:
#' Sugarcane is a **high polyploid** (the source study assumed ten
#' copies of each chromosome). The source file, however, holds
#' **diploid-style genotype calls** (`0/0`, `0/1`, `1/1`), because the
#' sequencing depth did not allow the allele dosage to be estimated. These
#' were recoded as `-1`, `0` and `1`. So no dosage exists in the data
#' to begin with: a `0` means only "both alleles seen", and the number of
#' copies of each allele is unknown (as is the dosage-dependent additive
#' effect). This data set cannot be used to simulate polyploid inheritance;
#' use it for diploid-style pipelines that need the marker structure of a
#' sugarcane breeding population.
#'
#' **Positions are sorghum coordinates.** No sugarcane reference genome was
#' used: reads were aligned to *Sorghum bicolor* v3.1, so `chr` (01 to 10)
#' and `pos` are sorghum positions, not sugarcane ones.
#'
#' @details
#' About 12% of the calls are heterozygous (`0`). The reference allele of
#' the VCF is the major allele. The source has no missing calls. The
#' genetic-map column `cm` is a synthetic map built from the (sorghum)
#' physical positions with `simplePHENOTYPES::synthetic_map()`.
#'
#' @format A data frame of 6,615 rows (SNPs) by 437 columns (5 map columns
#'   + 432 clones), each genotype in `{-1, 0, 1}`.
#' @source USDA Ag Data Commons, \doi{10.15482/USDA.ADC/1528986}
#'   (public domain, US Government work), file
#'   `sugarcane.10501.SNPs_.432.Inds_.vcf`; built by
#'   `data-raw/usda_sugarcane_ratoon.R`.
#' @references
#' Islam MS, Corak K, McCord P, Hulse-Kemp AM, Lipka AE (2023). A first
#' look at the ability to use genomic prediction for improving the
#' ratooning ability of sugarcane. *Frontiers in Plant Science* 14:1205999.
#' \doi{10.3389/fpls.2023.1205999}
#' @name Sugarcane_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' dim(marker_data("Sugarcane_div")) # bundled sample
#' \donttest{
#' geno <- marker_data("Sugarcane_div", full = TRUE) # whole panel (download)
#' if (!is.null(geno)) dim(geno)
#' }
NULL

#' Arabidopsis thaliana 1001 Genomes genotypes (sample only)
#'
#' @description
#' An *Arabidopsis thaliana* marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from the 1001 Genomes project (1,135
#' accessions, TAIR10 assembly; 1001 Genomes Consortium 2016). It is a
#' **sample only**: 2,000 SNPs by all 1,097 accessions, thinned from a
#' 6,000-SNP pilot panel by the same rule as the other samples (see
#' [marker_data()]), and returned by `marker_data("Arabidopsis_div")`.
#' markerdata does not host the whole data set, so
#' `marker_data("Arabidopsis_div", full = TRUE)` warns, returns `NULL` and
#' downloads nothing. **The original, whole data set is at
#' <https://1001genomes.org/data/GMI-MPI/releases/v3.1/>** (also the
#' `source_url` column of [marker_catalog()]).
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1` (reference homozygote / heterozygote /
#' alternate homozygote). *A. thaliana* is a selfing species and its
#' accessions are inbred, so the sample has no heterozygous calls and suits
#' additive schemes (cf. [SoySNP50K_div]). Missing calls in the source were
#' already imputed (the imputed "simulation" file of the pilot panel was
#' used); the sample has no `NA`.
#'
#' The genetic-map column `cm` is a synthetic map built from physical
#' positions with `simplePHENOTYPES::synthetic_map()` (the 1001 Genomes
#' release carries no cM); replace it with a real map if you have one.
#'
#' @format A data frame of 2,000 rows (SNPs, 5 chromosomes) by 1,102 columns:
#'   `snp` (`chr:pos:ref:alt`), `allele` (`REF/ALT`), `chr`, `pos`, `cm`,
#'   then 1,097 accessions (1001 Genomes accession ids, as character column
#'   names), each an integer genotype in `{-1, 0, 1}`.
#' @source 1001 Genomes, GMI-MPI release v3.1,
#'   <https://1001genomes.org/data/GMI-MPI/releases/v3.1/>
#'   (`1001genomes_snp-short-indel_only_ACGTN.vcf.gz`), through the
#'   already-processed pilot panel; built by `data-raw/s1database_panels.R`.
#'   Use is subject to the 1001 Genomes data usage policy (cite the paper
#'   below; acknowledge the Weigel laboratory, Max Planck Institute).
#' @references
#' 1001 Genomes Consortium (2016). 1,135 genomes reveal the global pattern
#' of polymorphism in *Arabidopsis thaliana*. *Cell* 166(2):481-491.
#' \doi{10.1016/j.cell.2016.05.063}
#' @name Arabidopsis_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' geno <- marker_data("Arabidopsis_div") # the sample (no download)
#' dim(geno)
#' marker_catalog()[marker_catalog()$name == "Arabidopsis_div", "source_url"]
NULL

#' Human 1000 Genomes Project phase 3 genotypes (sample only)
#'
#' @description
#' A human (*Homo sapiens*) marker matrix in the same numeric format as
#' [SNP55K_maize282_maf04], curated from the 1000 Genomes Project phase 3
#' release (2,504 people from 26 populations, GRCh37; 1000 Genomes Project
#' Consortium 2015). It is a **sample only**: 2,000 SNPs by all 2,504
#' people, thinned from a 6,000-SNP pilot panel by the same rule as
#' the other samples (see [marker_data()]), and returned by
#' `marker_data("Human_div")`. markerdata does not host the whole data set,
#' so `marker_data("Human_div", full = TRUE)` warns, returns `NULL` and
#' downloads nothing. **The original, whole data set is at
#' <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>** (also the
#' `source_url` column of [marker_catalog()]).
#'
#' The data are from publicly released, consented 1000 Genomes samples and
#' are included for method development and teaching.
#'
#' @details
#' Genotypes are coded `-1`/`0`/`1` (reference homozygote / heterozygote /
#' alternate homozygote). Humans are outbred: about 20% of calls in the
#' sample are heterozygous. The pilot panel has no missing calls.
#'
#' The genetic-map column `cm` is a synthetic map built from physical
#' positions with `simplePHENOTYPES::synthetic_map()`; it is not a measured
#' human recombination map (use a published one for real work).
#'
#' @format A data frame of 2,000 rows (SNPs, 22 autosomes) by 2,509 columns:
#'   `snp` (`chr:pos:ref:alt`), `allele` (`REF/ALT`), `chr`, `pos`, `cm`,
#'   then 2,504 people (1000 Genomes sample ids such as `HG00096`), each an
#'   integer genotype in `{-1, 0, 1}`.
#' @source 1000 Genomes Project phase 3, release 20130502,
#'   <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>
#'   (`ALL.chr*.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`),
#'   through the already-processed pilot panel; built by
#'   `data-raw/s1database_panels.R`. The IGSR makes these data available
#'   without embargo; cite the paper below.
#' @references
#' The 1000 Genomes Project Consortium (2015). A global reference for human
#' genetic variation. *Nature* 526:68-74. \doi{10.1038/nature15393}
#' @name Human_div
#' @seealso [marker_data()], [marker_catalog()], [marker_references()]
#' @family marker datasets
#' @keywords datasets
#' @examples
#' geno <- marker_data("Human_div") # the sample (no download)
#' dim(geno)
#' marker_catalog()[marker_catalog()$name == "Human_div", "source_url"]
NULL
