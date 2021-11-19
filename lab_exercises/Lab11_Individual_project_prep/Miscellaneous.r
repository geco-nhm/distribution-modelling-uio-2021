# BIOS 5211/9211
# Olav Skarpaas, NHM, University of Oslo, Nov 2021

# Miscellaneous
###############

# Repeating and expanding on some topics from previous lectures and labs,
# as well as questions from students, mostly providing some pointers.

# Topics roughly in order of modelling process, with some additional/advanced stuff at the end.

# R language
# R basics: Lab1
# Data types most often encountered: vectors, matrices, data.frames, lists

# Spatial data and maps: lectures Day 2 + Lab 2
# Key data type: raster
# Bioclim variables are based on historical data 1970-2000, future data on different Climate models
# and timeframes; see Documentation for bioclim variables: https://www.worldclim.org/data/bioclim.html

# Species occurrence data (GBIF) and generation of training data sets: Lab 2 + Lab 10
# Predictor variables should have same names in training data as in predictor rasters


# Data assessment:
# - correlation (collinearity): cor.test (pearson, kendall, spearman), wilcox.test
# - collinearity in maxent: correlation among transformed variables (see also single vs. marginal response plots)
# - FOPs and choice of transformations, Lectures day 4, Lab 4


# Misc model considerations:
# - methods for presence-absence (e.g. glm) vs. presence-only/presence-background (e.g. maxent)
# - interactions/communities: see lectures and Lab 7 on hmsc (Ryan Burner)
# - detection: see lectures and Lab 8 on occupancy modelling (Ryan Burner)
# - alternative methods for different response variables:
#    count data: e.g. glm with family=poisson or neg.bin
#    multinomial data (such as vegetation classes): e.g. random forest (e.g. sdm package)
# - model selection (also called variable selection):
#    a) manual (based on e.g. ecological knowledge, correlation, FOP-plot; Labs 3-4)
#    b) automated, e.g. stepwise forward selection based on AIC (Lab 3) or F-tests (Lab 4), or by
#       variance inflation:
library(usdm)
?vif
# - ensemble modelling (e.g. sdm), regularization (e.g. maxent), weighting, averaging... (be careful)


# Model validation and evaluation: (see also lectures and labs days 3-5)
# - Response curve plots: MIAmaxent::plotResp (single), plotResp2 (marginal)
# - Deviance, AIC
# - AUC: MIAmaxent::testAUC
# - subset vs. cross-validation: e.g. sdm package, Lab 5, exercise 2


# Other (advanced) stuff
# - Alternatives to MIAmaxent for fitting maxent models:
?dismo::maxent
?maxnet
# - Making R use several cores (parallelism) does not seem to be easy, but here are some hints:
#   https://www.bio21.unimelb.edu.au/hpc/index.php/guides/reference/r-parallelism



