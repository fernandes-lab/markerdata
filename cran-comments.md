## Submission

Initial submission of markerdata 0.1.0.

## Test environments

* local macOS, R 4.6.1
* GitHub Actions (ubuntu-latest release/devel/oldrel, windows-latest, macOS-latest)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for CRAN

* Large data sets are not shipped in the package; they are downloaded on
  first use from the package's GitHub release and cached under
  `tools::R_user_dir("markerdata", "cache")`. Download failures produce a
  warning and a NULL return, never an error, so examples and tests run offline.
