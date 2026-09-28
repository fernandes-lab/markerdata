# Moved verbatim from breedingDesigner R/data.R on 2026-09-28, with the .rda files
# (data/) and their build scripts (data-raw/). To update for this package:
# - each @section Distribution still describes breedingDesigner's .gitignore policy;
# - \link{SoySNP50K_div} / \link{EUChip60K_div} point at objects that stay in
#   breedingDesigner (use \link[breedingDesigner]{...} or move those two here too).
# - data-raw/t3_wheat_tgc.R (tetraploid TGC, the documented Wheat_div) and
#   data-raw/dryad_wheat_zhanggao.R (older hexaploid bread-wheat panel) both write
#   data/Wheat_div.rda; running the second overwrites the documented object.

#' Sorghum WEST Diversity Panel genotypes (GBS subset)
#'
#' A sorghum marker matrix in the same numeric format as
#' \code{simplePHENOTYPES::SNP55K_maize282_maf04}, so it can seed a \code{founders}
#' node for simulation. It is a curated subset of the sorghum \strong{WEST Diversity
#' Panel} genotyped by genotyping-by-sequencing (GBS) for the water-use-efficiency
#' GWAS/TWAS study of Ferguson et al. (2021): 842 accessions and 12,000 SNPs
#' (MAF >= 0.05 in the subset, thinned evenly across the 10 sorghum chromosomes).
#' Genotypes are coded \code{-1}/\code{0}/\code{1} (reference homozygote /
#' heterozygote / alternate homozygote).
#'
#' Like the inbred \code{\link{SoySNP50K_div}} accessions (and unlike the outbred
#' \code{\link{EUChip60K_div}} trees), these accessions are largely homozygous:
#' about 5\% of genotype calls are heterozygous, so the dataset suits additive
#' schemes (line development, genomic selection) rather than heterosis/combining-
#' ability schemes.
#'
#' The genetic-map column \code{cm} is a synthetic map built from physical
#' positions with \code{simplePHENOTYPES::synthetic_map()}; \code{snp} IDs are
#' synthesized as \code{S<chr>_<pos>} (the source matrix carries none). Replace the
#' map with a real one if you have it.
#'
#' @section Distribution: The source genotypes are openly licensed
#'   \strong{CC-BY-4.0} (redistributable with attribution). The built \code{.rda}
#'   is nonetheless kept \strong{out of git for now} (see \code{.gitignore}) so
#'   \pkg{breedingDesigner} does not grow large with genotype files: the example
#'   marker datasets are slated to move to a separate companion data package. Build
#'   it locally with \code{data-raw/sorghum_wue.R}.
#'
#' @format A data frame of 12,000 rows (SNPs) by 847 columns: the first five are
#'   \code{snp}, \code{allele}, \code{chr}, \code{pos}, \code{cm}; the remaining 842
#'   are accessions, each an integer genotype in \code{{-1, 0, 1}}.
#' @source WEST Diversity Panel GBS, Zenodo \doi{10.5281/zenodo.5019227}
#'   (license CC-BY-4.0). The bundled subset is built by
#'   \code{data-raw/sorghum_wue.R} from a recoded \code{-1}/\code{0}/\code{1}
#'   genotype matrix derived from that deposit.
#' @references
#' Ferguson JN, Fernandes SB, Monier B, Miller ND, Allen D, Dmitrieva A, Schmuker P,
#' Lozano R, Valluru R, Buckler ES, Gore MA, Brown PJ, Spalding EP, Leakey ADB
#' (2021). Machine learning-enabled phenotyping for GWAS and TWAS of WUE traits in
#' 869 field-grown sorghum accessions. \emph{Plant Physiology} 187(3):1481-1500.
#' \doi{10.1093/plphys/kiab346}
#'
#' Ferguson JN, Fernandes SB, Monier B, et al. (2020). WEST Diversity Panel GBS
#' information. \emph{Zenodo}, \doi{10.5281/zenodo.5019227} (CC-BY-4.0).
#' @keywords datasets
#' @examples
#' \dontrun{
#' pop <- simplePHENOTYPES::as_population(Sorghum_div, individuals = 1:100)
#' }
"Sorghum_div"

#' Rice IRRI tropical genomic-selection panel genotypes (GBS subset)
#'
#' A rice marker matrix in the same numeric format as
#' \code{simplePHENOTYPES::SNP55K_maize282_maf04}. Curated subset of the IRRI elite
#' tropical rice (indica / indica-admixed) genomic-selection breeding panel
#' genotyped by GBS (Spindel et al. 2015): 370 lines and 12,000 SNPs (MAF >= 0.05,
#' thinned evenly across the 12 rice chromosomes). Genotypes are coded \code{-1}/
#' \code{0}/\code{1}. These are inbred lines (~3\% heterozygous), so the dataset
#' suits additive schemes (line development, genomic selection).
#'
#' @section Distribution: Source genotypes are CC0 (Dryad; redistributable). The
#'   built \code{.rda} is kept out of git for now (see \code{.gitignore}); build it
#'   locally with \code{data-raw/dryad_rice_spindel.R}. Slated for the companion
#'   data package so \pkg{breedingDesigner} stays small.
#' @format A data frame of 12,000 rows (SNPs) by 375 columns (5 map columns +
#'   370 lines), each genotype in \code{{-1, 0, 1}}.
#' @source Dryad \doi{10.5061/dryad.7369p} (CC0), file
#'   \code{MET_crfilt_.75_allchrom.hmp.txt}; built by
#'   \code{data-raw/dryad_rice_spindel.R}.
#' @references
#' Spindel J, Begum H, Akdemir D, Virk P, Collard B, Redona E, Atlin G, Jannink J-L,
#' McCouch SR (2015). Genomic Selection and Association Mapping in Rice
#' (\emph{Oryza sativa}). \emph{PLoS Genetics} 11(2):e1004982.
#' \doi{10.1371/journal.pgen.1004982}
#' @keywords datasets
#' @examples
#' \dontrun{
#' pop <- simplePHENOTYPES::as_population(Rice_div, individuals = 1:100)
#' }
"Rice_div"

#' Cattle worldwide diversity panel genotypes (BovineSNP50 subset)
#'
#' A cattle marker matrix in the same numeric format as
#' \code{simplePHENOTYPES::SNP55K_maize282_maf04}. Curated subset of a worldwide
#' cattle diversity panel (134 breeds) genotyped on the Illumina BovineSNP50 array
#' (Decker et al. 2014): 1,543 animals and 12,000 SNPs (MAF >= 0.05, thinned evenly
#' across the 29 autosomes; UMD3.1 positions). Genotypes are coded \code{-1}/
#' \code{0}/\code{1}. Cattle are \strong{outbred} (~29\% heterozygous), so a
#' dominance layer can be simulated (cf. \code{\link{EUChip60K_div}}) --- useful for
#' animal-breeding schemes exploiting heterosis / combining ability.
#'
#' @section Distribution: Source genotypes are CC0 (Dryad; redistributable). The
#'   built \code{.rda} is kept out of git for now (see \code{.gitignore}); build it
#'   locally with \code{data-raw/dryad_cattle_decker.R}. Slated for the companion
#'   data package.
#' @format A data frame of 12,000 rows (SNPs) by 1,548 columns (5 map columns +
#'   1,543 animals), each genotype in \code{{-1, 0, 1}}.
#' @source Dryad \doi{10.5061/dryad.th092} (CC0), PLINK
#'   \code{ww_cattle_dryad_20Sept2013.ped/.map}; built by
#'   \code{data-raw/dryad_cattle_decker.R}.
#' @references
#' Decker JE, McKay SD, Rolf MM, Kim J, Molina Alcala A, Sonstegard TS, Hanotte O,
#' et al. (2014). Worldwide Patterns of Ancestry, Divergence, and Admixture in
#' Domesticated Cattle. \emph{PLoS Genetics} 10(3):e1004254.
#' \doi{10.1371/journal.pgen.1004254}
#' @keywords datasets
#' @examples
#' \dontrun{
#' pop <- simplePHENOTYPES::as_population(Cattle_div, individuals = 1:100)
#' }
"Cattle_div"

#' Poplar (Populus trichocarpa) landscape-genomics panel genotypes (34K subset)
#'
#' A poplar marker matrix in the same numeric format as
#' \code{simplePHENOTYPES::SNP55K_maize282_maf04}. Curated subset of a
#' \emph{Populus trichocarpa} landscape-genomics panel genotyped on the
#' PoplarConsortium-ORNL 34K SNP array (Geraldes et al. 2014): 434 accessions and
#' 12,000 SNPs (MAF >= 0.05, thinned evenly across the 19 nuclear chromosomes).
#' Genotypes are coded \code{-1}/\code{0}/\code{1}. These forest trees are
#' \strong{outbred} (~32\% heterozygous), so a dominance layer can be simulated
#' (cf. \code{\link{EUChip60K_div}}) --- useful for tree-breeding schemes.
#'
#' @section Distribution: Source genotypes are CC0 (Dryad; redistributable). The
#'   built \code{.rda} is kept out of git for now (see \code{.gitignore}); build it
#'   locally with \code{data-raw/dryad_poplar_geraldes.R}. Slated for the companion
#'   data package.
#' @format A data frame of 12,000 rows (SNPs) by 439 columns (5 map columns +
#'   434 accessions), each genotype in \code{{-1, 0, 1}}. SNP IDs are the array's
#'   \code{scaffold_<chr>_<pos>} names.
#' @source Dryad \doi{10.5061/dryad.7s848} (CC0), file
#'   \code{FOR_DRYAD_TandB_fwd_june2013_434.txt}; built by
#'   \code{data-raw/dryad_poplar_geraldes.R}.
#' @references
#' Geraldes A, Farzaneh N, Grassa CJ, McKown AD, Guy RD, Mansfield SD, Douglas CJ,
#' Cronk QCB (2014). Landscape genomics of \emph{Populus trichocarpa}.
#' \emph{Evolution} 68(11):3260-3280. \doi{10.1111/evo.12497}
#'
#' Geraldes A, DiFazio SP, Slavov GT, et al. (2013). A 34K SNP genotyping array for
#' \emph{Populus trichocarpa}. \emph{Molecular Ecology Resources} 13(2):306-323.
#' \doi{10.1111/1755-0998.12056}
#' @keywords datasets
#' @examples
#' \dontrun{
#' pop <- simplePHENOTYPES::as_population(Poplar_div, individuals = 1:100)
#' }
"Poplar_div"

#' Wheat Tetraploid Global Collection genotypes (90K SNP subset)
#'
#' A wheat marker matrix in the same numeric format as
#' \code{simplePHENOTYPES::SNP55K_maize282_maf04}. Curated subset of the
#' \strong{Tetraploid wheat Global Collection (TGC)} --- wild and domesticated
#' emmer, durum wheat landraces and tetraploid subspecies --- genotyped on the
#' wheat 90K iSelect SNP array (Maccaferri et al. 2019): 1,854 accessions and
#' 12,000 SNPs (MAF >= 0.05, thinned evenly across the 14 nuclear chromosomes of
#' the A and B subgenomes, 1A--7B). Genotypes are coded \code{-1}/\code{0}/\code{1}.
#' Tetraploid wheat lines are largely homozygous (~8\% heterozygous), so the
#' dataset suits additive schemes (cf. \code{\link{SoySNP50K_div}}).
#'
#' \strong{Polyploidy caveat.} Durum/emmer wheat is an \emph{allotetraploid}
#' (2n = 4x = 28, AABB). The 90K SNPs are subgenome-specific, so each is treated
#' here as an ordinary \strong{diploid} locus (the deposited VCF already calls
#' genotypes \code{0/0}, \code{0/1}, \code{1/1}, recoded to \code{-1}/\code{0}/
#' \code{1}) --- the standard treatment for wheat SNP arrays, but an
#' \emph{approximation} that ignores homoeologous relationships and tetrasomic
#' inheritance. A faithful representation needs a tetrasomic / polyploid dosage
#' model; that belongs to the \strong{polyploid framework planned in
#' \pkg{simplePHENOTYPES}} (the engine), to which this dataset should migrate.
#' Until then, use \code{Wheat_div} for additive, diploid-style simulation only.
#'
#' @section Distribution: Distributed via T3/Wheat (Triticeae Toolbox), a public
#'   USDA small-grains resource; cite Maccaferri et al. (2019). The built
#'   \code{.rda} is kept out of git for now (see \code{.gitignore}); build it
#'   locally with \code{data-raw/t3_wheat_tgc.R}. Slated for the companion data
#'   package (confirm T3 redistribution terms at migration).
#' @format A data frame of 12,000 rows (SNPs) by 1,859 columns (5 map columns +
#'   1,854 accessions), each genotype in \code{{-1, 0, 1}}. SNP IDs are the 90K
#'   array's \code{IWB} probe names (made unique across subgenomes).
#' @source T3/Wheat (Triticeae Toolbox), file
#'   \code{..._TGC-GTC_collection_SNP90K_mapped_filtered.vcf}; built by
#'   \code{data-raw/t3_wheat_tgc.R}. Reference per USDA GrainGenes
#'   ("Global Durum Genomic Resources").
#' @references
#' Maccaferri M, Harris NS, Twardziok SO, Pasam RK, Gundlach H, Spannagl M, ...
#' Cattivelli L, Tuberosa R (2019). Durum wheat genome highlights past
#' domestication signatures and future improvement targets. \emph{Nature Genetics}
#' 51:885-895. \doi{10.1038/s41588-019-0381-3}
#' @keywords datasets
#' @examples
#' \dontrun{
#' pop <- simplePHENOTYPES::as_population(Wheat_div, individuals = 1:100)
#' }
"Wheat_div"
