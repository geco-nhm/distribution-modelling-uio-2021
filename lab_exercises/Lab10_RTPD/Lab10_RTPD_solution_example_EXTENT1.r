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
# key <- name_backbone(name= [WRITE YOUR SPECIES NAME HERE], kingdom= [WRITE THE APPROPORIATE KINGDOM HERE])$speciesKey
 key <- name_backbone(name= "Acer ginnala", kingdom="Plantae")$speciesKey
# key <- name_backbone(name= "Syringa josikaea", kingdom="Plantae")$speciesKey
# key <- name_backbone(name= "Neovison vison", kingdom="Animalia")$speciesKey
# key <- name_backbone(name= "Oncorhyncus gorbuscha", kingdom="Animalia")$speciesKey
# key <- name_backbone(name= "Phasianus colchincus", kingdom="Animalia")$speciesKey
sp <- occ_search(taxonKey=key, hasCoordinate=TRUE, country="NO", limit=1000)

# Modify and repeat the following for all three extents:

# Load environmental predictor data as rasters
# r.list <- list.files("lab_data/RTPD/TIFF_EXTENT/NORWAY", pattern="tif$", full.names=TRUE)
r.list <- list.files("lab_data/RTPD/TIFF_GRAIN/2km", pattern="tif$", full.names=TRUE)
predictors.Norway <- stack(r.list)

# Generate training data for your species for the particular extent
training.data.Norway <- generate.training.data(sp,predictors.Norway,generate.absences=TRUE,absence.value=NA,n.abs=1000,factors=c("artype","norge123"))   # To see how the function 'generate.training.data' works, open the file 'Dataprep.functions.r'

# Inspect the data, select an appropriate distribution modelling method,
# fit a model, evaluate it and collect the information you need to answer the questions a-e above.
# Reuse code from other labs!

# a. Observed response frequencies
summary(training.data.Norway)
training.data.Norway <- training.data.Norway[complete.cases(training.data.Norway[,-1]),]  # Remove missing values(NAs), except for in 'presence' (where NA means absent)
summary(training.data.Norway)
par(mfrow=c(4,4))
for(i in 2:15) plotFOP(training.data.Norway,i)

# b. Environmental variable covariation
plot(training.data.Norway)
str(training.data.Norway)
par(mfrow=c(1,1))
corrplot::corrplot.mixed(cor(training.data.Norway[,-c(1,4,10)],use="complete.obs"),lower="number",upper="ellipse",outline=TRUE,tl.pos="lt",tl.cex=0.90,tl.col = "black",cl.cex = 0.90)
table(training.data.Norway$ar50.artype.2km, training.data.Norway$geo.norge123.2km)

# c. Model structure (selected environmental variables) 
# Here a rather random selection of variables - you should select based on correlations etc. for each specific species and data set!
DV <- deriveVars(training.data.Norway[,c("presence","ar50.artype.2km","bioclim.12.2km",
                                         "Growing.season.length.2km","swe.5.2km","Topographic.Wetness.Index.2km",
                                         "Total.insolation.2km")])
DVforEV <- selectDVforEV(DV$dvdata)
EV <- selectEV(DVforEV$dvdata)
ms.extent.1 <- EV$selectedmodel  # Store models for other extents by other name
summary(ms.extent.1)

# d. Modelled response curves (for selected environmental variables)
# Code for plotting some responses - all variables will not be included in models for all species, check your model!
par(mfrow=c(3,2))
plotResp(ms.extent.1,DV$transformations,"Growing.season.length.2km")
plotResp(ms.extent.1,DV$transformations,"Total.insolation.2km")
plotResp(ms.extent.1,DV$transformations,"bioclim.12.2km")
plotResp(ms.extent.1,DV$transformations,"swe.5.2km")
plotResp(ms.extent.1,DV$transformations,"ar50.artype.2km")

# e. Model predictions (spatial predictions, AUC)
par(mfrow=c(1,1))
ps.extent.1 <- projectModel(ms.extent.1,DV$transformations,predictors.Norway) # May take some time...
#plot(ps.extent.1$output,main="Model predictions") # run this if you want to plot again without projecting the model once more
i <- training.data.Norway$presence==1
points(training.data.Norway$x[i],training.data.Norway$y[i],pch=".",cex=2,col="red") # Presence points
i <- is.na(training.data.Norway$presence)
points(training.data.Norway$x[i],training.data.Norway$y[i],pch=".",cex=2,col="blue") # Absence points

# testAUC with TRAINING DATA
# Ideally we would like an independent data set, here we instead
# use the training data to test the model - You could also subset the data before fitting, and do this repeatedly (cross-validation)
test_data <- training.data.Norway
test_data$presence[is.na(test_data$presence)] <- 0  # recoding background points (NA) as "true" absences (0)
summary(test_data)
testAUC(ms.extent.1,transformations=DV$transformations,data=test_data)



# Make sure to store the results/plots for each extent before you move on
# to the next, so that you can compare results across extents (i.e. answer the questions a-e above).

# Join with the others in your group to make a presentation of the results,
# as described above.



