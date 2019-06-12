Analysis
========
All the R scripts in this directory run the MSDRP pipeline. The scripts expect the following input:

1. Gene expression data in a given disease context.
2. Gene expression data comparing normal condition with perturbed agents such as chemicals and drugs.
3. GWAS data in a given disease context.

Scripts
-------
1. Calculate the significant overlap between/among different data sets using Fisher's exact test.

.. code-block:: sh

   $ cd R/analysis/
   $ Rscript SignificantOverlap.R ../../data/

2. Run pathway enrichment analysis using Signaling Pathway Impact (SPIA) for each of the disease genes.

.. code-block:: sh

   $ cd R/analysis/
   $ Rscript DiseaseSPIA.R ../../data/

3. Conduct pathway enrichment analysis using SPIA for drug perturbed genes in each disease.

.. code-block:: sh

   $ cd R/analysis/
   $ Rscript DrugSPIA.R ../../data/

4. Generate correlation and dissimilarity scores for all drug-disease pathway pairs.

.. code-block:: sh

   $ cd R/analysis/
   $ Rscript CheckDistribution.R ../../data/

5. Create pathway activities with the combination of two drugs.

.. code-block:: sh

   $ cd R/analysis/
   $ Rscript DrugCombination.R ../../data/
