# BIOS 5211/9211
# Olav Skarpaas, NHM, University of Oslo, Nov 2021

# Lab3. Distribution modelling with presence-absence data
#########################################################

# In this lab we will improve the modelling process illustrated briefly in Lab1,
# and detailed in the lecture by Rune Halvorsen earlier today. We will work with
# spatial data prepared as in Lab2, but with lower resolution.

# Load libraries
library(raster)
library(fields)    # for color ramps (e.g. tim.colors)
library(corrplot)  # useful correlation plot (reference: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html)

# Data are available in the 'lab_data' folder
# ===========================================
# The zipfile 'RASTER_2km_10km.zip' contains two folders with raster data:
# 2km: predictor rasters with 2km-resolution
# 10km: predictor rasters with 10km-resolution
# Unzip and put the two folder in 'lab_data/RASTER'


# Exercise 1: A model for Sitka spruce in Norway
# ==============================================

# Ecological model
# ----------------
# This first, important part of the modelling process is skipped here, but see
# introductory lecture day 1 (Olav Skarpaas) and Halvorsen 2012, Sommmerfeltia 34,
# and keep it in mind. Think of the purpose of modelling (ERM or SPM)
# and the relevance/importance of predictor variables.


# Data model: data structuring - picking up from Lab2
# ---------------------------------------------------

path <- "lab_data/RASTER/2km/" # Change this to 10km if you have problems with large files

# Load predictor raster stack (maps of predictors for prediction and plotting)
r.list <- list.files(path, pattern="tif$", full.names=TRUE) # add all raster files in a given folder to a list
predictors <- stack(r.list)                                 # read and stack raster layers
names(predictors)                                           # Names of raster layers (maps of predictor variables)
names(predictors) <- gsub("_2km","",names(predictors))      # Rename: get rid of "_2km" (or "_10km"), to match variable names in training data
names(predictors)                                           # Check renaming
plot(predictors)                                            # Plots of raster layers - see Lab2 for more info on each layer

# Load, inspect and clean training data (for regression analysis)
load("lab_data/Norway_sitka_training_data_2km")  # Change this to 10km if you have problems with large files
names(training_data)
summary(training_data)                                         # Note missing values
training_data <- training_data[complete.cases(training_data),] # Remove missing values
any(is.na(training_data))                                      # Check: no missing values left

# Plot map of training data, with elevation as background
par(mfrow=c(1,1))
plot(predictors$dem100,main="Presence/absence points on elevation map",col=terrain.colors(64))
points(training_data$x,training_data$y,col=c("blue","red")[training_data$presence+1]) # Presences and random absences on elevation map

# Biplots of all variables in training data set
plot(training_data)

# Continuous variables: check correlations
# Correlation plot for training data set, removing response (presence) and factor variables
str(training_data)
corrplot.mixed(cor(training_data[,-c(3,14,15)],use="complete.obs"),lower="number",upper="ellipse",outline=TRUE,tl.pos="lt",tl.cex=0.90,tl.col = "black",cl.cex = 0.90)

# Which predictors will you include/exclude, considering your
# ecological understanding of the species, the relationship between predictors
# and the quality of the data (ecological model + data model)?

# Categorical variables: check combinations
table(training_data$ar50_artype, training_data$geo_norge123)

# Zero observations for (combinations of) levels may cause trouble for
# model fitting and prediction (e.g. in the 'predict' function for rasters).
# Consider simplifying levels and/or using only one of the factor variables.

# For land cover (ar50) - how does resolution (e.g. 10km vs. 100m) affect
# classification of pixels with presences?


# Statistical model
# -----------------
# Logistic regression with GLM: A presence-absence DM for 
# Sitka spruce in Norway.

# a) Build a preliminary "full" regression model containing
# all the predictors you selected above, using the 'glm' function.
# (take a look at Lab1, or type '?glm' in the console, if you don't remember how to specify the model)

# b) Non-spatial model validation: Inspect model coefficients, tests of terms,
# residual deviance vs. null deviance. Is it a good model?

# c) Spatial model validation: Make predictions and plot prediction map
# Realistic predictions? Compare to species observations at gbif.org or artsdatabanken.no

# d) Model selection: Do stepwise backward model selection with AIC, using the 'step' function
# Which terms are removed from the full model by this procedure? Why?

# e) Plot and compare prediction maps for full model and selected model

# f) Look at plots and summaries of the final model and compare to global model (Lab1)
# Which predictors are important for this species?
# What is its response to environmental variables according to the national and global models?
# How / why do they differ?


# Model evaluation with TRAINING data
#------------------------------------

# g) Evaluate the selected model with the 'evaluate' function
# and ROC plot for the full training data in package 'dismo'.
# Is this a good model?


# Model evaluation with data partitioning
#----------------------------------------

# h) Set aside random subset of training_data (20%) for evaluation

# i) Refit the selected model WITHOUT THE EVALUATION DATA set aside in exercise 1h.

# j) Evaluate the refitted selected model from 3b with 'evaluate' function
# and ROC plot in package 'dismo'. Is the fit as good as for the training data (exercise 1g)?

# k) Repeat steps h-j. Do you get the same result? Why / why not?


# Exercise 2: Repeat Exercise 1 with other species (except 1f)
#=================================================

# Find and load training data for other species in the folder "lab_data":
# lilac: Hungarian lilac - Syringa josikaea
# boar: Wild boar - Sus scrofa

# Consider the following questions:
#  - which predictor variables are important?
#  - signs of overfitting or lack of explanatory/predictive power?
#  - are the final models plausible? why, why not?
#  - do the final models make good predictions? why, why not?



# Exercise 3: Independent study - do if time allows
#==============================

# a) If you have your own presence-absence data, apply the methods in this lab to your data.

# b) Download gbif data (see Lab2) for Norway or the World for a species you are interested in,
#    generate random absences (see Lab2), and deveolop a model with glm as in exercise 1.

# c) Take a look at tools for model validation in the sdm package.

# d) Prepare for tomorrow's lab by exploring the MIAmaxent package.
