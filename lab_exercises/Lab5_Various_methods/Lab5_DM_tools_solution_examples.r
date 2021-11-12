# BIOS 5211/9211
# Olav Skarpaas, NHM, University of Oslo, Nov 2021

# Lab5. Additional tools for classical distribution modelling
#############################################################

# In this lab we will explore some tools for development, validation and
# evaluation of statistical distribution models in various R packages.

# Load libraries
#----------------
library(raster)


# National data for Norway
#-------------------------
# The environmental data set for Norway in this lab is a rescaled subset of
# the data compiled and described in Horvath et al. 2019. Distribution 
# modelling of vegetation types based on area frame survey
# data. Applied Vegetation Science, 22(4): 547-560 (see syllabus).
# (See also Lab 2)

# Load data
#----------

# Load training data from Norway
load("lab_data/Norway_sitka_training_data_2km")
summary(training_data)
training_data <- training_data[complete.cases(training_data),] # Remove missing values

# Load predictor raster stack (maps of predictors for prediction and plotting)
path <- "lab_data/RASTER/2km/"
r.list <- list.files(path, pattern="tif$", full.names=TRUE) # add all raster files in a given folder to a list
predictors <- stack(r.list)                                 # read and stack raster layers
names(predictors)                                           # Names of raster layers (maps of predictor variables)
names(predictors) <- sub("_2km","",names(predictors))       # Rename (get rid of "_2km") to match variable names in training_data
plot(predictors)                                            # Plots of raster layers - see Lab2 for more info on each layer


# Exercise 1: Variable exploration and model validation in MIAmaxent
#-------------------------------------------------------------------

# a) From the training data, prepare a data set for MIAmaxent with response variable
# in first column and predictor variables (EVs) in the second, third column, etc.
# and absence coded as NA (i.e. random background rather than true absence)
names(training_data)
maxent_data <- training_data[,c(3:ncol(training_data))]
maxent_data$presence[maxent_data$presence==0] <- NA  # recoding absences to NA (zero is interpreted as true absence)
head(maxent_data)

# b) Load MIAmaxent and plot frequency of observed presence (FOP) for the variables.
# Which variables show ecologically reaonable/interpretable patterns?
# What kinds of transformations could be suitable for these variables?
library(MIAmaxent)
par(mfrow=c(3,4))
for(i in 2:13) plotFOP(maxent_data,i)

# c) Do forward model selection with MIAmaxent and summarize and plot model results (as in Lab 4)
# How does the final maxent model compare to glms for Sitka spruce Picea sitchensis in previous labs?
# How can you evaluate the model predictions?
# How might climate change affect this species, according to the two models?
DV <- deriveVars(maxent_data[,c("presence","ar50_artype","bioclim_12","geo_norge123",
                                "Growing_season_length","swe_5","Topographic_Wetness_Index",
                                "Total_insolation")])
DVforEV <- selectDVforEV(DV$dvdata)
EV <- selectEV(DVforEV$dvdata)
ms.maxent <- EV$selectedmodel

# Response curve plots
par(mfcol=c(2,1))
plotResp(ms.maxent,DV$transformations,"Growing.season.length")
plotResp(ms.maxent,DV$transformations,"bioclim.12")

# Prediction plots
par(mfrow=c(1,1))
ps.maxent <- projectModel(ms.maxent,DV$transformations,predictors)
par(mfrow=c(1,2))
plot(ps.maxent$output,main="Maxent model")
points(training_data$x,training_data$y,pch=".",cex=2,col=c("blue","red")[training_data$presence+1]) # Presences and random absences on elevation map
plot(ps,main="Selected glm")  # Make sure ps contains predictions for Picea sitchensis (run relevant lines in Lab3 script if necessary)
points(training_data$x,training_data$y,pch=".",cex=2,col=c("blue","red")[training_data$presence+1]) # Presences and random absences on elevation map

# The model can be evaluated with 'testAUC' for independent data

# Climate change effects:
summary(ms)        # Positive effect: warmer temperatures (run relevant lines in Lab3 script to get ms)
summary(ms.maxent) # Difficult to interpret coefficients directly.
                   # Based on response curve plots: Unimodal response to growing season, more precip positive if not wet already.
# In sum: maxent model seems to capture some nonlinear responses to temp and precip that the glm misses,
# but are they realistic? FOP-plots suggest that trends at high values of bioclim_12 may be unreliable.



# Exercise 2: Model validation by cross validation with functions in 'sdm'
#-------------------------------------------------------------------------

# a) Preapare an sdm data object with the 'sdmData' function in pacakge 'sdm'.
# Use the formula of the selected glm for Picea sitchensis (ms, Lab3, exercise 1d)
library(sdm)
model.formula <- formula(ms) # storing the formula of the model to be tested, to save some typing in the next commands
sdm_data <- sdmData(model.formula,train=training_data)

# b) Data partitioning by subsampling.
# Fit the selected model in exercise 2a with the function 'sdm' in package 'sdm'.
# Use arguments methods='glm', replication='sub', test.percent=20 and n=10.
# If sdm throws an error about missing methods, try 'installAll()'.
# Plot the ROC with the function 'roc'. What does this say about model performance?
ms.sub <- sdm(model.formula,data=sdm_data,methods='glm',replications='sub',test.percent=20,n=10) # If 
roc(ms.sub)

# c) Cross validation.
# Repeat exercise 2b with the argument replication='cv' and 
# cv.folds=5 (instead of test.percent=20). Compare results to 4b: How / why do they differ?
ms1.cv <- sdm(model.formula,data=sdm_data,methods='glm',replications='cv',cv.folds=5,n=10)
roc(ms1.cv)

# d) Reflect on the pros and cons of validation with subsets and crossvalidation within
# a dataset (2b and c) and validation with a completely independent data set. Discuss with
# your neighbor.



# Exercise 3: Model ensembles - for the bold! The topic will be introduced on Monday
#-----------------------------------------------------------------------------------

# a) Fit ensemble DM to Picea sitchensis (same model formula as exercise 1a) with
# several methods in sdm, with no replication/cross-validation (to save time).
# If sdm throws an error about missing methods, try 'installAll()'.

me <- sdm(model.formula,data=sdm_data,methods=c('glm','rf'),replications='cv',cv.folds=5,n=1)

# b) Compare and interpret ROC plots for the different methods
roc(me)

# c) Make spatial predictions and plot maps (computer intensive)
p <- predict(me,predictors)
plot(p)
pe <- ensemble(me,newdata=predictors,filename="ens1.img",setting=list(method="weighted",stat="AUC",wtest=c("test.dep")),overwrite=TRUE)
plot(pe)

# d) Compare ensemble model predictions to predictions from maxent and selected glm (exercise 1c)
# Prediction plots
par(mfrow=c(1,3))
plot(ps.maxent$output,main="Maxent model")
points(training_data$x,training_data$y,pch=".",cex=2,col=c("blue","red")[training_data$presence+1]) # Presences and random absences on elevation map
plot(ps,main="Selected glm")  # Make sure ps contains predictions for Picea sitchensis (run relevant lines if necessary)
points(training_data$x,training_data$y,pch=".",cex=2,col=c("blue","red")[training_data$presence+1]) # Presences and random absences on elevation map
plot(pe,main="Ensemble glm+rf")
points(training_data$x,training_data$y,pch=".",cex=2,col=c("blue","red")[training_data$presence+1]) # Presences and random absences on elevation map




# Exercise 4: More species
#-------------------------

# Repeat exercises 1-3 for Hungarian lilac (Syringa jokasiea) and Wild boar (Sus scrofa) (see Lab 3)




# Exercise 5: Exploration
#------------------------

# Explore options for ensemble model fitting, prediction validation and evaluation
# in various packages, eg. 'sdm', 'dismo', 'biomod2' and 'AICcmodavg'

