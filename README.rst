MSDRP (A Modular and Multi-Scale Drug Repositioning Workflow)
=============================================================
This package comprises a modular workflow designed to identify drug repositioning candidates using multi omics data
sets.

Workflow
--------
A schematic figure of the workflow is presented below. The R scripts necessary to run the MSDRP pipeline are located
`in the R directory <https://github.com/asifemon/msdrp/tree/master/R>`_.


.. image:: https://github.com/asifemon/msdrp/blob/master/data/img/workflow.jpg
    :width: 500px

Installation
------------
The required R libraries and dependencies to run the workflow are detailed in the
`requirements.txt file <https://github.com/asifemon/msdrp/blob/master/requirements.txt>`_.

Case Scenario
-------------
The case scenario uses different datasets that are described `in the preprocessing folder <https://github.com/asifemon/msdrp/tree/master/R/preprocessing>`_.

How to Modify the Workflow
--------------------------

Notes how to change parts of the workflow:

1. Selecting different gene sets (i.e., "gene set intersection" part in the figure)
2. Modifying the Pathway enrichment Analysis method (i.e., GSEA instead of SPIA)
