# Multi-environment validation notes

Trait: GY; environments: ICRISAT_2014-15 and ICRISAT_2015-16.
2800 accessions have grain-yield phenotypes in both seasons.
Grouped CV1 holds out all records from marker-derived PC clusters; CV2 holds out line-by-season cells while retaining the same line in the other season.
The GxE model is a ridge model with marker main effects plus marker-by-season interactions; environment is fixed.
Leave-one-season-out tests the same ICRISAT location across two years, using paired lines only. It is not leave-one-location-out validation.
Split SDs summarize repeated random partitions, not population confidence intervals. Marker PCA clusters are defined from genotype data without phenotypes.
Additional-year/site GxE, independent-location prediction, and pedigree/family external validation are not supported by these data alone.
CicerSeq source files are downloaded locally and are not redistributed in the repository.
