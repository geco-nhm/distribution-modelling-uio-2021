# BIOS 5211/9211
# Olav Skarpaas, Eva Lieungh, Peter Horvath and Adam Naas,
# NHM, University of Oslo, Nov 2021

# Lab10. RTPD
#############

# The students are split into groups that will address effects of
# ONE of these data properties (see notes on data etc. for each
# theme below):
# 1. Extent of study area
# 2. Grain of study area
# 3. Density of presence points
# 4. Clumping of presence points

# Each person in the group works with a different species, and addresses
# the following questions:
# How does extent/resolution/density/clumping [ONE of these] affect
# a. Observed response frequencies
# b. Environmental variable covariation
# c. Model structure (selected environmental variables) 
# d. Modelled response curves (for selected environmental variables)
# e. Model predictions (spatial predictions, AUC)

# Compare around the table: how do results for different species compare?
# How/why do they differ?

# Make a presentation (about 10 min in total) of individual results
# (2 mins per person), as well as summary across the group (2 mins).

# Species list (one per person):
# ------------------------------
# Amur maple (sibirlønn) - Acer ginnala
# Hungarian lilac (ungarsk syrin) - Syringa josikaea
# American mink (mink) - Neovison vison
# Humpback (pukkellaks) - Oncorhynchus gorbuscha
# Common Pheasant (fasan) - Phasianus colchicus



# Load libraries
library(rgbif)
library(raster)
library(sp)
library(MIAmaxent)

# Source script with functions
source('lab_exercises/Lab10_RTPD/Dataprep_functions.r')



# Theme 1: EXTENT
#################

# You will work with three extents:
# 1. All of Norway
# 2. Southern or northern Norway (depending on distribution of species)
# 3. South-eastern or south-western Norway (depending on species, but extent may not be appropriate for all species)

# Here is the task in words:

# First, download a data set of occurrences for your species from GBIF.
# Then, for each extent, do the following:
# - Load predictor rasters (2km resolution) for the appropriate extent (modify paths and raster stack names as needed)
# - Generate a training data set based on the occurrences and the predictor rasters
# - Fit an appropriate distribution model and assess it with respect to a-e above

# Here is some draft code:

# Download GBIF data
key <- name_backbone(name= [WRITE YOUR SPECIES NAME HERE], kingdom= [WRITE THE APPROPORIATE KINGDOM HERE])$speciesKey
sp <- occ_search(taxonKey=key, hasCoordinate=TRUE, country="NO", limit=1000)

# Modify and repeat the following for all three extents:

# Load environmental predictor data as rasters
r.list <- list.files("lab_data/RTPD/TIFF_EXTENT/NORWAY", pattern="tif$", full.names=TRUE)
predictors.Norway <- stack(r.list)

# Generate training data for your species for the particular extent
training.data.Norway <- generate.training.data(sp,predictors.Norway,generate.absences=TRUE,absence.value=NA,n.abs=1000,factors=c("artype","norge123"))   # To see how the function 'generate.training.data' works, open the file 'Dataprep.functions.r'

# Inspect the data, select an appropriate distribution modelling method,
# fit a model, evaluate it and collect the information you need to answer the questions a-e above.
# Reuse code from other labs!

# Make sure to store the results/plots for each extent before you move on
# to the next, so that you can compare results across extents (i.e. answer the questions a-e above).

# Join with the others in your group to make a presentation of the results,
# as described above.



# Theme 2: GRAIN
################

# You will work with four different grains for the Oslo region:
# 1. 100 m
# 2. 1 km
# 3. 2 km
# 4. 10 km

# Here is the task in words:

# First, download a data set of occurrences for your species from GBIF.
# Then, for each grain, do the following:
# - Load predictor rasters for the appropriate grain (modify paths and raster stack names as needed)
# - Generate a training data set based on the occurrences and the predictor rasters
# - Fit an appropriate distribution model and assess it with respect to a-e above

# Here is some draft code:

# Download GBIF data
key <- name_backbone(name= [WRITE YOUR SPECIES NAME HERE], kingdom= [WRITE THE APPROPORIATE KINGDOM HERE])$speciesKey
sp <- occ_search(taxonKey=key, hasCoordinate=TRUE, country="NO", limit=1000)

# Modify and repeat the following for all three grains:

# Load environmental predictor data as rasters
r.list <- list.files("lab_data/RTPD/TIFF_GRAIN/100m", pattern="tif$", full.names=TRUE)
predictors.100m <- stack(r.list)

# Generate training data for your species for the particular grain
training.data.100m <- generate.training.data(sp,predictors.100m,generate.absences=TRUE,absence.value=NA,n.abs=1000,factors=c("artype","norge123"))   # To see how the function 'generate.training.data' works, open the file 'Dataprep.functions.r'

# Inspect the data, select an appropriate distribution modelling method,
# fit a model, evaluate it and collect the information you need to answer the questions a-e above.
# Reuse code from other labs!

# Make sure to store the results/plots for each grain before you move on
# to the next, so that you can compare results across grains (i.e. answer the questions a-e above).

# Join with the others in your group to make a presentation of the results,
# as described above.



# Theme 3: DENSITY
##################

# You will work with three different densities of occurrence data for each of the species:
# 1. original density (all occurrence data)
# 2. 50% density (50% of the occurrence data)
# 3. 10% density (10% of the occurrence data)

# Here is the task in words:

# First, load predictor rasters for Norway with 2km resolution (modify paths and raster stack names if needed)
# Second, download a data set of occurrences for your species from GBIF.
# Then, for each density of occurrences, do the following:
# - Resample the occurrence data set to the appropriate density
# - Generate a training data set based on the occurrences and the predictor rasters
# - Fit an appropriate distribution model and assess it with respect to a-e above

# Here is some draft code:

# Load environmental predictor data as rasters (Norway, 2km resolution)
r.list <- list.files("lab_data/RTPD/TIFF_EXTENT/NORWAY", pattern="tif$", full.names=TRUE)
predictors <- stack(r.list)

# Download GBIF data
key <- name_backbone(name= [WRITE YOUR SPECIES NAME HERE], kingdom= [WRITE THE APPROPORIATE KINGDOM HERE])$speciesKey
sp <- occ_search(taxonKey=key, hasCoordinate=TRUE, country="NO", limit=1000)

# Modify and repeat the following for all three densities:

# Thin the occurrence data to the appropriate density by sampling without replacement (.100pst indicates 100%)
# (modify 'thin' as needed)
n <- nrow(sp$data)
thin <- 1
sp.100pst <- sp 
sp.100pst$data <- sp$data[sample.int(round(n*thin),replace=FALSE),]

# Generate training data for your species for the particular density (.100pst indicates 100%)
training.data.100pst <- generate.training.data(sp.100pst,predictors,generate.absences=TRUE,absence.value=NA,n.abs=1000,factors=c("artype","norge123"))   # To see how the function 'generate.training.data' works, open the file 'Dataprep.functions.r'

# Inspect the data, select an appropriate distribution modelling method,
# fit a model, evaluate it and collect the information you need to answer the questions a-e above.
# Reuse code from other labs!

# Make sure to store the results/plots for each density before you move on
# to the next, so that you can compare results across densities (i.e. answer the questions a-e above).

# Join with the others in your group to make a presentation of the results,
# as described above.


# Theme 4: CLUMPING
###################

# You will work with three different degrees of clumping of occurrences for each of the species:
# 1. as original data
# 2. moderately clumped (around 5 focal points)
# 3. strongly clumped (around the same 5 focal points)


# Here is the task in words:

# First, load predictor rasters for Norway with 2km resolution (modify paths and raster stack names if needed)
# Second, download a data set of occurrences for your species from GBIF.
# Third, generate a training data set based on the occurrences and the predictor rasters
# Then, for each degree of clumping of occurrences, do the following:
# - Resample the training data set to the appropriate degree of clumping
# - Fit an appropriate distribution model and assess it with respect to a-e above

# Here is some draft code:

# Load environmental predictor data as rasters (Norway, 2km resolution)
r.list <- list.files("lab_data/RTPD/TIFF_EXTENT/NORWAY", pattern="tif$", full.names=TRUE)
predictors <- stack(r.list)

# Download GBIF data
key <- name_backbone(name= [WRITE YOUR SPECIES NAME HERE], kingdom= [WRITE THE APPROPORIATE KINGDOM HERE])$speciesKey
sp <- occ_search(taxonKey=key, hasCoordinate=TRUE, country="NO", limit=1000)

# Generate training data for your species based on occurrences and predictor rasters
training.data <- generate.training.data(sp,predictors,generate.absences=TRUE,absence.value=NA,n.abs=1000,factors=c("artype","norge123"))   # To see how the function 'generate.training.data' works, open the file 'Dataprep.functions.r'

# Resample the data to the appropriate level of clumping by resampling 50% of the data without replacement:
N <- nrow(training.data)
n <- round(N/2)                                                       # number of rows to sample (50% of data)
foci <- training.data[sample.int(N,size=5,replace=FALSE),c("x","y")]  # random sampling foci

# 1. No clumping
training.data.no.clump <- training.data[sample.int(N,size=n,replace=FALSE),]

# 2. Moderate clumping: resampling within 100 km of 5 focal points
dmax <- 100000 # max distance (m) from focal point
d <- sp::spDists(as.matrix(foci),as.matrix(training.data[,c("x","y")]))
selected.points <- apply(d<dmax,2,any)                           # points that are within dmax m of any focal point
training.data.mod.clump <- training.data[selected.points,]
N <- nrow(training.data)
if(N<n) n <- N
training.data.mod.clump <- training.data.mod.clump[sample.int(N,size=n,replace=FALSE),]

# 3. Strong clumping
# reuse code for point 2, but with lower dmax


# For each degree of clumping, inspect the data, select an appropriate distribution modelling method,
# fit a model, evaluate it and collect the information you need to answer the questions a-e above.
# Reuse code from other labs!

# Make sure to store the results/plots for each degree of clumping before you move on
# to the next, so that you can compare results across degrees of clumping (i.e. answer the questions a-e above).

# Join with the others in your group to make a presentation of the results,
# as described above.

