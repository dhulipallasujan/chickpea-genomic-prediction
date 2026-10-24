# All-site yield-per-plant validation notes

Trait harmonization: YPP in all sheets except ICRISAT, where the supplied column is YPPlnt.
The workbook supplied 11 site-year environments across 6 locations.
3 duplicate normalized accession-by-site-year rows were averaged before validation.
Grouped CV1 holds all records from marker-derived genotype clusters in one fold; CV2 holds out line-by-site-year cells while retaining that line in another environment.
The GxE model uses marker main effects plus marker-by-site-year interactions and fixed effects for observed site-year environments.
Leave-one-location-out holds out every year from one location. Its model has no target-site effect or target-site interaction, because neither can be estimated without target-site phenotypes.
Location-specific scale and management differences are therefore a substantial challenge; interpret both RMSE and ranking/calibration metrics.
Three repeated split sets summarize partition variability, not population confidence intervals. Related groups are approximate PCA/k-means clusters, not recorded pedigree families.
Yield columns are named differently by site in the source workbook; site-specific trait mapping and counts are reported in all_site_environment_counts.csv.
The ridge penalty grid was 1 to 100,000; 0% of nested CV tuned folds selected a grid endpoint. Interpret tuning accordingly.
CicerSeq source files are downloaded locally and are not redistributed in this repository.
The all-site marker model uses every 20th QC-passed marker (50 markers); the full sampled QC set has 953 markers.
