# markerdata 0.1.0

* Initial release. Bundles `SNP55K_maize282_maf04` whole and ships a sample
  (2,000 SNPs, 1,000 for the four largest; every individual) of every other data set; see
  `marker_catalog()` for the full list (crops, trees, cattle, sheep, Atlantic
  salmon, and the sample-only human and Arabidopsis sets). The whole panels
  are hosted in two GitHub releases, `data-v1` and `data-v2`.
  `Potato_div` (autotetraploid) is stored diploidized and `Sugarcane_div`
  holds the diploid-style calls of its source; their help topics explain
  what is lost.
  `marker_data(name)` returns the sample offline; `marker_data(name,
  full = TRUE)` downloads and caches the whole panel, and
  `marker_download()` fetches all of them. `Arabidopsis_div` (1001 Genomes)
  and `Human_div` (1000 Genomes Project phase 3) are sample only: markerdata
  hosts no whole panel for them, `full = TRUE` warns and gives the original
  download, and the new `source_url` catalog column holds that link for every
  data set. Also provides a catalog
  (`marker_catalog()`) and citation objects (`marker_references()`).
