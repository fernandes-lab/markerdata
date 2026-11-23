# TODO

## Before the next package release

- [ ] Publish GitHub release data-v2 (seven assets) before the next package release.
  The exact `gh release create data-v2 ...` command is in `data-raw/README.md`.
  Until then `marker_data(name, full = TRUE)` for the seven data-v2 panels
  warns and returns `NULL` (their samples work).

## Autopolyploids and allele dosage

The `-1`/`0`/`1` format cannot hold dosages `0..ploidy` for autopolyploids
(potato, alfalfa, sugarcane, blueberry). Until a dosage format exists, any
autopolyploid panel must be stored diploidized (all heterozygous dosage
classes collapsed to `0`) and say so in its documentation.

Done so far:

- `Potato_div` (autotetraploid, 586 clones) is stored diploidized: dosage 0 ->
  -1, dosage 4 -> 1, dosages 1-3 -> 0, so 69% of its calls are `0`. Its help
  topic has a section on what is lost.
- `Sugarcane_div` (high polyploid, 432 clones) is diploid-called by its
  source (`0/0`, `0/1`, `1/1`): no dosage exists in the source. Positions are
  Sorghum bicolor v3.1 coordinates. Its help topic says so.

- [ ] Open: design a dosage-capable format. It affects `.check_marker()`, the
  2-bit packing in `R/pack.R`, `simplePHENOTYPES` and `breedingDesigner`.
- [ ] Open: `Arabidopsis_div` and `Human_div` are sample-only; decide whether
  to host a larger thinned panel in the GitHub release.

## Species still without a data set

- [ ] blueberry: Dryad doi:10.5061/dryad.kd4jq6h (CC0). Needs a logged-in
  manual download into `data-raw/blueberry/`. The two 56-byte error files
  there are not data.
- [ ] chicken: Figshare 16627099. The `.bed` holds 1,512 samples but the
  `.fam` lists 839: ask the author or find another set. The script
  `data-raw/figshare_chicken_li.R` is ready.
- [ ] pig: candidate Dryad doi:10.5061/dryad.t1r3d (CC0); positions
  unconfirmed.
- [ ] cassava: candidate Figshare doi:10.6084/m9.figshare.14773005.v1
  (CC BY 4.0).
- [ ] conifer: candidate Sitka spruce, Dryad doi:10.5061/dryad.ghx3ffbv6;
  linkage-group positions only.
- [ ] canola: candidate Zenodo 4302088; contents unknown.
- [ ] alfalfa: candidate doi:10.57745/L0FLJD; license unknown.
- [ ] barley: no open set with 200+ lines found.
- [ ] hexaploid wheat larger than 159 lines: none found with positions
  (`WheatHex_div` has 159).
