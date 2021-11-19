# BIOS 5211/9211
# Olav Skarpaas, NHM, University of Oslo, Nov 2021


# Individual student project template
#####################################

# Read the instructions and write your code with annotations in this template.
# (Transferring the script to R markdown is ok, if you prefer.)
# In addition to this file, you will use code and data from previous labs, as well as
# some new global data, as described below.

# Please check that you have the latest version of the course material available on GitHub.

# You will work with ONE of these species:
# 1. Abutilon theophrasti
# 2. Euphorbia esula
# 3. Galanthus nivalis
# 4. Helianthus annuus
# 5. Impatiens glandulifera


# Load libraries
#---------------
# NB! the libraries listed here are the libraries you need for data preparation
# - you may need to load additional libraries necessary for your choice of methods.
library(rgbif)
library(raster)
library(sp)


# Source script with functions for data preparation
#--------------------------------------------------
source('lab_exercises/Lab10_RTPD/Dataprep_functions.r')


# Environmental data
#-------------------

# The World
# We will use 'bioclim', i.e. Biologically relevant climate variables: http://www.worldclim.org/bioclim.
# The fowllowing code takes a subset of these variables that are also available in the data set for
# Norway (bioclim 1, 10 and 12; see below). Feel free to get additional data if you like, or modify
# the code to increase the resolution etc., but it is also perfectly fine to use the code as it is.
# Load data for current climate (predictors) and for projected future climate (scenario RCP8.5) in
# the same resolution (set by the argument res):
predictors.global <- raster::getData('worldclim',var='bio',res=10)[[c(1,10,12)]]
scenarios.global <- raster::getData('CMIP5', var='bio', res=10, rcp=85, model='NO', year=50)[[c(1,10,12)]]
names(scenarios.global) <- names(predictors.global) # renaming scenarios to match names in predictors

# Norway
# The environmental data set for Norway is a scaled subset of
# the data compiled and described in Horvath et al. 2019. Distribution 
# modelling of vegetation types based on area frame survey
# data. Applied Vegetation Science, 22(4): 547-560 (see syllabus).
# See Lab 2 for definitions of variables.
# Load environmental predictor rasters for Norway, 2 km resolution (modify path if needed):
r.list <- list.files("lab_data/RTPD/TIFF_EXTENT/NORWAY", pattern="tif$", full.names=TRUE)
predictors.Norway <- stack(r.list)


# Training data
#--------------
# We prepare training data for distribution modelling by downloading species occurrence data (presences) from GBIF,
# adding random pseudo-absence/background points (as in Lab 2, but with absences coded as NA),
# and extracting environmental variables for presence/absence points from the
# environmental rasters above.
# You can work with the code and data as they are, or, if you like, modify the script or write
# your own code to download additional occurrence data or build the training data differently.

# Download GBIF data for the World
key <- name_backbone(name= "[WRITE YOUR SPECIES NAME HERE]", kingdom= "Plantae")$speciesKey
sp.global <- occ_search(taxonKey=key, hasCoordinate=TRUE, limit=1000)

# Generate global training data
training.data.global <- generate.training.data(sp.global,predictors.global,generate.absences=TRUE,absence.value=NA,n.abs=1000)   # To see how the function 'generate.training.data' works, open the file 'Dataprep.functions.r'

# Download GBIF data for Norway
key <- name_backbone(name= "[WRITE YOUR SPECIES NAME HERE]", kingdom= "Plantae")$speciesKey
sp.Norway <- occ_search(taxonKey=key, hasCoordinate=TRUE, country="NO", limit=1000)

# Generate training data for Norway
training.data.Norway <- generate.training.data(sp.Norway,predictors.Norway,generate.absences=TRUE,absence.value=NA,n.abs=1000)   # To see how the function 'generate.training.data' works, open the file 'Dataprep.functions.r'



# Independent work
#-----------------

# Fill inn the script with your own code (reuse code from labs!) to address
# the project tasks (see also lecture last day):

# - Use an appropriate distribution modelling technique to estimate
#   responses to climate variables and predict the potential
#   distribution of the species across the globe, now and with future climate change.

# - Mke a model for the species based on environmental predictors
#   and occurrences in Norway, and assess how observed and predicted occurrence
#   of the species varies with respect to climate, land cover and other environmental conditions
#   in Norway (to the extent that these are represented by environmental variables in the data set).

# - How good are the models for the global and national ranges, and how do they compare in terms of
#   modelled ecological responses and predictions?

# - Based on results from the modelling above, discuss
#     How is the distribution affected by climate?
#     How is the distribution affected by land use?
#     How might the distribution of the species respond to future changes in
#     - global climate?
#     - land use in Norway?

# - Based on existing knowledge of the species and your results in this study, discuss
#   limitations of the models and (briefly) how further research on the species may be aided by
#   other modelling approaches and/or other kinds of studies.

# Annotate the script with comments to explain what you are doing.

# When done, please rename this script 'Individual_project_[species name]_[your name].r'
# and upload the script in canvas along with your report. (Alternatively,
# if you prefer to work in R markdown, please provide the markdown file.)

# Thanks!

