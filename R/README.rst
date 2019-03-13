How to Use
==========

Preprocessing
-------------

This directory contains all the preprocessing scripts to download and harmonize the dataset showcased in the manuscript.

Pipeline
--------

All the R scripts in this directory run the MSDRP pipeline. The scripts expect the following input:

1. Gene expression data in a given disease context.
2. Gene expression data comparing normal condition with perturbed agents such as chemicals and drugs.
3. GWAS data in a given disease context.

Scripts
-------
1. **SignificantOverlap.R:** calculates the significant overlap between/among different data sets using Fisher's exact test. 
2. **DiseaseSPIA.R:** calculates pathway enrichment analysis using Signaling Pathway Impact (SPIA) for each of the Disease Genes. 
3. **DrugSPIA.R:** conducts pathway enrichment analysis using Signaling Pathway Impact (SPIA) for Drug perturbed genes in each Disease. 
4. **CheckDistribution.R:** generates correlation and dissimilarity scores for all drug-disease pathway pairs.
5. **DrugCombination.R:** create pathway activities with the combination of two drugs. 
