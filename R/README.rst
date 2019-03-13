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
1. `SignificantOverlap.R: <https://github.com/asifemon/msdrp/blob/master/R/1-SignificantOverlap.R>`_ calculates the significant overlap between/among different data sets using Fisher's exact test. 
2. `DiseaseSPIA.R: <https://github.com/asifemon/msdrp/blob/master/R/2-DiseaseSPIA.R>`_ calculates pathway enrichment analysis using Signaling Pathway Impact (SPIA) for each of the Disease Genes. 
3. `DrugSPIA.R: <https://github.com/asifemon/msdrp/blob/master/R/3-DrugSPIA.R>`_ conducts pathway enrichment analysis using Signaling Pathway Impact (SPIA) for Drug perturbed genes in each Disease. 
4. `CheckDistribution.R: <https://github.com/asifemon/msdrp/blob/master/R/4-CheckDistribution.R>`_ generates correlation and dissimilarity scores for all drug-disease pathway pairs.
5. `DrugCombination.R: <https://github.com/asifemon/msdrp/blob/master/R/5-DrugCombination.R>`_ create pathway activities with the combination of two drugs. 
