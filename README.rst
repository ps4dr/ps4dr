MSDRP (Multi-scale Drug Repositioning Pipeline)
===============================================


Installation
------------

The required R libraries and dependencies to run the workflow are detailed in the
requirements.txt file.

How to Download the Data
========================

GWAS data sets for STOPGAP:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We deployed STOPGAP pipeline (https://github.com/StatGenPRD/STOPGAP) to merge and harmonize GWAS data sets from publicly
available data sources.


1. GRASP: Genome-Wide Repository of Associations Between SNPs and Phenotypes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Associations between SNPs and phenotypes in different disease conditions, curated by the NIH

- Web Page: https://grasp.nhlbi.nih.gov/Overview.aspx
- Download: https://s3.amazonaws.com/NHLBI_Public/GRASP/GraspFullDataset2.zip
- License: ???
- Citation: Leslie R, O’Donnell CJ, Johnson AD (2014) GRASP: analysis of genotype-phenotype results from 1,390
  genome-wide association studies and corresponding open access database. Bioinformatics 30(12), i185-94. GRASP
  Build 2.0.0.0

2. GWAS Catalog: The NHGRI-EBI Catalog of published genome-wide association studies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Web Page: https://www.ebi.ac.uk/gwas/home
- Download: https://www.ebi.ac.uk/gwas/api/search/downloads/full
- License: The NHGRI-EBI GWAS Catalog and all its contents are available under the general Terms of Use for
  EMBL-EBI Services.
- Citation: Buniello A, MacArthur JAL, Cerezo M, Harris LW, Hayhurst J, Malangone C, McMahon A, Morales J, Mountjoy E,
  Sollis E, Suveges D, Vrousgou O, Whetzel PL, Amode R, Guillen JA, Riat HS, Trevanion SJ, Hall P, Junkins H, Flicek P,
  Burdett T, Hindorff LA, Cunningham F and Parkinson H.
  The NHGRI-EBI GWAS Catalog of published genome-wide association studies, targeted arrays and summary statistics 2019.
  Nucleic Acids Research, 2019, Vol. 47 (Database issue): D1005-D1012.

3. GWAS DB: an update database for human genetic variants identified by genome-wide association studies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GWASdb is a one stop shop which combines collections of traits/diseases associated SNP (TASs) from current GWAS and
their comprehensive functional annotations, as well as disease classifications.

- Web Page: http://jjwanglab.org/gwasdb
- Download: ftp://147.8.193.36/GWASdb/gwasdb_20150819_snp_trait.gz
- License: Data of GWASdb is ONLY provided for research purposes and MUST NOT be used as other purpose, such as
  commercial usage and resources abuse, before contacting us.
- Citation: Li, Mulin Jun, et al. "GWASdb v2: an update database for human genetic variants identified by genome-wide
  association studies." Nucleic acids research 44.D1 (2015): D869-D876.

4. PHEWAS Catalog: The NHGRI-EBI Catalog of published genome-wide association studies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Web Page: https://phewascatalog.org/phewas
- Download: https://phewascatalog.org/files/phewas-catalog.csv.zip
- License:
- Citation: Denny JC, Bastarache L, Ritchie MD et al. Systematic comparison of phenome-wide association study of
  electronic medical record data and genome-wide association study data. Nat Biotechnol. 2013 Dec;31(12):1102-10.


Data Preprocessing
------------------

Data sets:

1.
2.
3.

The scripts to process the data sets used in the case scenario of the manuscript are located in the
/data folder.

