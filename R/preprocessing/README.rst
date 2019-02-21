Preprocessing
=============

How to run the preprocessing scripts.

1. **01STOPGAP2_functions.R**. Contains all the functions needed to run the next script
2. **02STOPGAP2_run.R**. Run STOPGAP2 pipeline to generate processed-GWAS data.
3. **03mesh2efoSTOPGAP.R**. Converts MESH terms to EFO identifiers in the processed-GWAS data for harmonizing purposes.
4. **04getOpenTargetsData.R**. Downloading all data from Open Targets using its API.
5. **05processHarmonizome.R**. Download the LINCS L1000 dataset and maps LINCS identifiers to ChEMBL ids for harmonizing purposes.
