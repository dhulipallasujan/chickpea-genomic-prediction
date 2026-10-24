# Genomic prediction result notes

Target: GY in ICRISAT_2014-15.
Matched accessions: 2807; retained genome-wide SNPs: 194.
Five-fold cross-validation repeated 5 times; split-to-split SD is not a population confidence interval.
Phenotypes are accession means, not plot-level records. This is a within-environment CV1 benchmark, not independent external validation.
Read one marker in every 500 HapMap rows, applied MAF/missingness QC, then retained every 5th QC-passed marker for the compact model.
Mean RMSE reduction versus training-mean baseline: 8.9%.
Mean Pearson r: 0.412; Spearman rho: 0.409; calibration slope (observed ~ predicted): 1.02.
RR-BLUP used standardized markers; ridge penalty was selected from 1, 3.16227766016838, 10, 31.6227766016838, 100, 316.227766016838, 1000, 3162.27766016838, 10000, 31622.7766016838, 1e+05 by inner three-fold CV within each training fold.
CicerSeq asks users to follow the spirit of the Toronto Statement; cite the source and do not redistribute the source files.
