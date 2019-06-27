#!/usr/bin/env bash

####################
# 1. Preprocessing #
####################

# See: https://github.com/ps4dr/ps4dr/blob/master/R/preprocessing/README.rst#preprocess-gwas-with-stopgap2
Rscript R/preprocessing/STOPGAP2_run.R data/

# See: https://github.com/ps4dr/ps4dr/blob/master/R/preprocessing/README.rst#prepare-mesh-to-efo-mappings
Rscript R/preprocessing/MeSH2EFO.R data/

# See: https://github.com/ps4dr/ps4dr/blob/master/R/preprocessing/README.rst#download-open-targets-data
Rscript R/preprocessing/RetrieveDEGs.R data/

# See: https://github.com/ps4dr/ps4dr/blob/master/R/preprocessing/README.rst#download-and-preprocess-lincs
Rscript R/preprocessing/RetrieveDrugResponseData.R data/

###############
# 2. Analysis #
###############

Rscript R/analysis/SignificantOverlap.R data/

Rscript R/analysis/DiseaseSPIA.R data/

Rscript R/analysis/DrugSPIA.R data/

Rscript R/analysis/CheckDistribution.R data/

Rscript R/analysis/DrugCombination.R data/
