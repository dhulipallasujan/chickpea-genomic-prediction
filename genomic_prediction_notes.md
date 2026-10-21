# Genomic prediction result notes

Target: GY in ICRISAT_2014-15.
Matched accessions: 2807; retained genome-wide SNPs: 194.
Five-fold cross-validation repeated 5 times; split-to-split SD is not a population confidence interval.
Phenotypes are accession means, not plot-level records. This is a within-environment CV1 benchmark, not independent external validation.
Read one marker in every 500 HapMap rows, applied MAF/missingness QC, then retained every 5th QC-passed marker for the compact model.
Mean RMSE reduction versus training-mean baseline: 8.4%.
Mean Pearson r: 0.402; Spearman rho: 0.402; calibration slope (observed ~ predicted): 0.909.
RR-BLUP used standardized markers; ridge penalty was selected from 1, 10, 100, 1000 by inner three-fold CV within each training fold.
CicerSeq asks users to follow the spirit of the Toronto Statement; cite the source and do not redistribute the source files.
