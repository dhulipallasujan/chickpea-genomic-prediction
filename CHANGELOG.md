# Changelog

## 2.2.0

- Reran the two-season analysis with the expanded ridge-penalty grid and regenerated the results table and figures from the saved predictions.
- Added all-site yield-per-plant analysis for 11 site-year sheets across six locations, including grouped CV1, known-line CV2, and leave-one-location-out validation.
- Saved leave-one-location and leave-one-season out-of-fold predictions; added independent checks that reproduce validation summaries from saved predictions.
- Generalized marker-by-environment features to support more than two environments and repaired contradictory ridge-boundary test expectations.
- Added a compact all-site marker panel, base-R figures, and GitHub exclusions for raw data, local libraries, caches, and R session files.

## 2.1.0

- Fixed the ridge penalty grid (1 to 1e5, log-spaced); v2.0 selected the 1000 ceiling in every fold. Added grid-edge warnings and QC reporting.
- Reported all CV2 gains against the environment mean; the known-line baseline is shown separately because it is worse than the mean.
- Withdrew the marker-density conclusion pending rerun.
- Added `rrBLUP::mixed.solve` cross-check, season-to-season yield correlation, generated results tables and ggplot2 figures, and ridge-behaviour tests.
- Replaced the truncated MIT licence with the standard text; removed leftover scripts from an earlier unrelated project; stopped tracking the 47 MB prediction file (now compressed and git-ignored); CI now syntax-checks every script.

## 2.0.0

- Added two-season grouped CV1, known-line CV2, marker-by-season ridge models, and leave-one-season-out checks.
- Added genotype-derived group folds, fold-logic tests, marker-density sensitivity, and multi-environment outputs.

## 1.0.0

- Single-season CicerSeq genomic-prediction benchmark with marker QC, nested ridge tuning, repeated CV1 and a training-mean comparator.
