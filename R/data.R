# Documentation for markerdata's marker data sets.
#
# SNP55K_maize282_maf04 is the only object actually bundled with the package
# (data/SNP55K_maize282_maf04.rda), so it gets a standard data-documentation
# block ending in the string "SNP55K_maize282_maf04". The other seven panels
# below are help TOPICS only (each block ends in NULL and carries an explicit
# @name): they are not package objects, and their data is downloaded on first
# use with marker_data("<name>") from the package's GitHub release, then
# cached locally (see marker_cache_dir() / marker_cache_clear()).
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
#' the only dataset actually bundled with markerdata; the seven larger
#' panels documented elsewhere in this package ([SoySNP50K_div],
#' [EUChip60K_div], [Sorghum_div], [Rice_div], [Cattle_div], [Poplar_div],
#' [Wheat_div]) are downloaded on first use with [marker_data()].
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
#' USDA Soybean Germplasm Collection: 300 accessions and 12,000 SNPs
#' (MAF >= 0.05 in the subset) on the Wm82.a2 (`Wm82.gnm2`) assembly. This
#' dataset is **not** shipped inside the markerdata package; it is
#' downloaded on first use with `marker_data("SoySNP50K_div")` (a packed,
#' xz-compressed download of about 1 MB) and cached locally afterward (see
#' [marker_cache_dir()]).
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
#' @format A data frame of 12,000 rows (SNPs) by 305 columns: the first
#'   five are `snp`, `allele`, `chr`, `pos`, `cm`; the remaining 300 are
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
#' \donttest{
#' geno <- marker_data("SoySNP50K_div")
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
#' subset). This dataset is **not** shipped inside the markerdata package;
#' it is downloaded on first use with `marker_data("EUChip60K_div")` (a
#' packed, xz-compressed download of about 2 MB) and cached locally
#' afterward.
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
#' \donttest{
#' geno <- marker_data("EUChip60K_div")
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
#' accessions and 12,000 SNPs (MAF >= 0.05 in the subset). This dataset is
#' **not** shipped inside the markerdata package; it is downloaded on first
#' use with `marker_data("Sorghum_div")` (a packed, xz-compressed download
#' of about 3 MB) and cached locally afterward.
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
#' \donttest{
#' geno <- marker_data("Sorghum_div")
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
#' subset). This dataset is **not** shipped inside the markerdata package;
#' it is downloaded on first use with `marker_data("Rice_div")` (a packed,
#' xz-compressed download of about 1 MB) and cached locally afterward.
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
#' \donttest{
#' geno <- marker_data("Rice_div")
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
#' positions). This dataset is **not** shipped inside the markerdata
#' package; it is downloaded on first use with `marker_data("Cattle_div")`
#' (a packed, xz-compressed download of about 5 MB) and cached locally
#' afterward.
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
#' \donttest{
#' geno <- marker_data("Cattle_div")
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
#' 0.05 in the subset). This dataset is **not** shipped inside the
#' markerdata package; it is downloaded on first use with
#' `marker_data("Poplar_div")` (a packed, xz-compressed download of about
#' 1 MB) and cached locally afterward.
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
#' \donttest{
#' geno <- marker_data("Poplar_div")
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
#' the A and B subgenomes, 1A--7B). This dataset is **not** shipped inside
#' the markerdata package; it is downloaded on first use with
#' `marker_data("Wheat_div")` (a packed, xz-compressed download of about
#' 6 MB) and cached locally afterward.
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
#' \donttest{
#' geno <- marker_data("Wheat_div")
#' if (!is.null(geno)) dim(geno)
#' }
NULL
