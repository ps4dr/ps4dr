PS4DR (Pathway Signatures for Drug Repositioning)
=================================================
This package comprises a modular workflow designed to identify drug repositioning candidates using multi-omics data
sets.

Installation
------------
If ``devtools`` is not installed, do:

.. code-block:: sh

   $ R -e 'install.packages("devtools", repos="http://cran.us.r-project.org")'

This library can be installed directly from GitHub using the following instructions adapted
from `R Package primer <https://kbroman.org/pkg_primer/pages/github.html>`_:

.. code-block:: sh

   $ R -e 'library(devtools); install_github("ps4dr/ps4dr")'

Alternatively, ``ps4dr`` can be cloned then installed from the source with:

.. code-block:: sh

   $ R -e 'install.packages(c("BiocManager", "data.table", "doParallel", "doSNOW", "dplyr", "foreach", "gplots", "gridExtra", "httr", "jsonlite", "pROC", "RecordLinkage", "reshape", "riverplot", "tidyr", "tidyverse"))'
   $ R -e 'BiocManager::install(c("BiocParallel", "EnsDb.Hsapiens.v86", "SPIA", "graphite"))'

Then

.. code-block:: sh

   $ git clone https://github.com/ps4dr/ps4dr.git
   $ R -e 'library(devtools); install("ps4dr")'

Reproduction
------------
To run the entire pipeline:

.. code-block:: sh

   $ sh reproduction.sh

Alternatively, see the instructions to:

1. Run all pre-processing scripts using the instructions at
   https://github.com/ps4dr/ps4dr/tree/master/R/preprocessing
2. Run all analysis scripts using the instructions at
   https://github.com/ps4dr/ps4dr/tree/master/R/analysis

Results
-------
A schematic figure of the workflow is presented below. The R scripts necessary to run the MSDRP pipeline are located
`in the R directory <https://github.com/ps4dr/ps4dr/tree/master/R>`_.

.. image:: https://github.com/ps4dr/ps4dr/blob/master/data/img/workflow.jpg
    :width: 500px

**Figure** 1. Design of the MSDRP workflow. Differentially expressed genes/proteins (i.e., DEG/DEP) from disease and
drug perturbed profiles are passed as input together with GWAS data. Once the data is correctly formatted, users can
define a custom pipeline, or series of steps in the workflow that will then be applied to the datasets. The steps
performed in this pipeline constitute the optional portion of the workflow and involve filtering the *-omics* features
coming from the dataset in order to reduce dimensionality by exclusively analyzing genes that have been associated with
GWAS studies. Next, a previously selected pathway enrichment method is applied to DEG/DEP datasets deriving from both
the disease and drug perturbed profiles to evaluate the direction of dysregulation for each affected pathway in each of
these contexts. Finally, the workflow prioritizes drugs by finding the drugs that are predicted to invert the pathway
signatures observed in the pathophysiology context.

How to Modify the Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~
Notes how to change parts of the workflow:

1. Selecting different gene sets (i.e., "gene set intersection" part in the figure)
2. Modifying the Pathway enrichment Analysis method (i.e., GSEA instead of SPIA)
