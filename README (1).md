# Public inputs (downloaded locally)

The genomic-prediction workflow uses files from the [ICRISAT CicerSeq portal](https://cegresources.icrisat.org/cicerseq/?page_id=3444):

- Phenotypes: `3000_Phenotyping_data_QC.xlsx` ([direct file](https://cegresources.icrisat.org/cicerseq/phenotypic_files/3000_Phenotyping_data_QC.xlsx))
- Cultivated HapMap markers: `Cultivated_Ca1.hmp.txt.gz` through `Cultivated_Ca8.hmp.txt.gz` ([Ca1 example](https://cegresources.icrisat.org/cicerseq/hapmap_files/Cultivated/Cultivated_Ca1.hmp.txt.gz))

Run `Rscript run_genomic_prediction.R` from the repository root. Missing input files are downloaded automatically. These large raw files are excluded from Git and the project ZIP; do not redistribute them. The portal asks users to follow the spirit of the Toronto Statement and cite the data source and associated study. See the root README for the study citation, analysis details, and limitations.

The workbook and marker files are not stored in this repository. The project analyzes grain yield from ICRISAT_2014-15 and ICRISAT_2015-16, which are two seasons at one location..
