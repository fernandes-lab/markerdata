## Submission

Initial submission of markerdata 0.1.0.

## Test environments

* local macOS, R 4.6.1
* GitHub Actions (ubuntu-latest release/devel/oldrel, windows-latest, macOS-latest)

## R CMD check results

0 errors | 0 warnings | 1 note (new submission)

The package tarball is about 4.6 MB: it ships a small sample of each data set
(2,000 SNPs, or 1,000 for the four largest, for every individual).

## Notes for CRAN

* Large data sets are not shipped whole in the package: a small sample of
  each is included, and the whole panel is downloaded only when the user
  asks for it (`full = TRUE`) from the package's GitHub releases (`data-v1`,
  `data-v2`) and cached under
  `tools::R_user_dir("markerdata", "cache")`. Download failures produce a
  warning and a NULL return, never an error, so examples and tests run offline.
