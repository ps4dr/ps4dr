Preprocessing
=============
How to run the pre-processing scripts:

1. **01STOPGAP2_functions.R**. Contains all the functions needed to run the next script
2. **02STOPGAP2_run.R**. Run STOPGAP2 pipeline to generate processed-GWAS data.
3. **03mesh2efoSTOPGAP.R**. Converts MESH terms to EFO identifiers in the processed-GWAS data for harmonizing purposes.
4. **04getOpenTargetsData.R**. Downloading all data from Open Targets using its API.
5. **05processHarmonizome.R**. Download the LINCS L1000 dataset and maps LINCS identifiers to ChEMBL ids for harmonizing purposes.

How to Download the Data
------------------------
GWAS data sets for STOPGAP
~~~~~~~~~~~~~~~~~~~~~~~~~~
We deployed STOPGAP pipeline (https://github.com/StatGenPRD/STOPGAP) to merge and harmonize GWAS data sets from publicly
available data sources.

1. GRASP: Genome-Wide Repository of Associations Between SNPs and Phenotypes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Associations between SNPs and phenotypes in different disease conditions, curated by the NIH

- Web Page: https://grasp.nhlbi.nih.gov/Overview.aspx
- Download: https://s3.amazonaws.com/NHLBI_Public/GRASP/GraspFullDataset2.zip
- License: ???

.. [1] Leslie R, *et al.* (2014) GRASP: analysis of genotype-phenotype results from 1,390
       genome-wide association studies and corresponding open access database. Bioinformatics 30(12), i185-94.

2. GWAS Catalog: The NHGRI-EBI Catalog of published genome-wide association studies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Web Page: https://www.ebi.ac.uk/gwas/home
- Download: https://www.ebi.ac.uk/gwas/api/search/downloads/full
- License: The NHGRI-EBI GWAS Catalog and all its contents are available under the general Terms of Use for
  EMBL-EBI Services.

.. [2] Buniello A, *et al.* The NHGRI-EBI GWAS Catalog of published genome-wide association studies,
       targeted arrays and summary statistics 2019. Nucleic Acids Research, 2019, Vol. 47 (Database issue):
       D1005-D1012.

3. GWAS DB: an update database for human genetic variants identified by genome-wide association studies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GWASdb is a one stop shop which combines collections of traits/diseases associated SNP (TASs) from current GWAS and
their comprehensive functional annotations, as well as disease classifications.

- Web Page: http://jjwanglab.org/gwasdb
- Download: ftp://147.8.193.36/GWASdb/gwasdb_20150819_snp_trait.gz
- License: Data of GWASdb is ONLY provided for research purposes and MUST NOT be used as other purpose, such as
  commercial usage and resources abuse, before contacting us.

.. [3] Li, Mulin Jun, *et al.* (2015) GWASdb v2: an update database for human genetic variants identified by
       genome-wide association studies." Nucleic acids research 44.D1: D869-D876.

4. PHEWAS Catalog: The NHGRI-EBI Catalog of published genome-wide association studies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Web Page: https://phewascatalog.org/phewas
- Download: https://phewascatalog.org/files/phewas-catalog.csv.zip
- License:

.. [4] Denny JC, *et al.* (2013) Systematic comparison of phenome-wide association study of
       electronic medical record data and genome-wide association study data. Nat Biotechnol. Dec;31(12):1102-10.
