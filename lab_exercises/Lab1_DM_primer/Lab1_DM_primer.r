# BIOS 5211/9211
# Olav Skarpaas & Eva Lieungh, NHM, University of Oslo, Nov 2021

# Lab1: A DM primer & R basics
##############################

# In this lab we will run quickly through the steps in distribution modelling
# from a ready-made training data set to model output and map predictions.
# WARNING! This exerceise neither demonstrates all aspects of good modelling practice,
# nor great models - it is just a quick and simple start.

# Preparations (Lab0):
# Download and install R: https://cran.r-project.org/
# Download and install R studio: https://www.rstudio.com/
# Create an R project (File menu)
# Download Lab material and put in your working directory:
# data: https://uio-my.sharepoint.com/:u:/g/personal/oskarpaa_uio_no/EZuM02rDipxAuWaFStbBKI0BQZGkxBbWLRi20XSabejyoA?e=WhKRBM
#   script: github
# Download and install packages:
# install.packages(c("raster", "fields"))


# A DM primer
# -----------

# Load packages
library(raster) # for raster data tools
library(fields) # for tim.colors and other color palettes

# Load training data
load("lab_data/DM_primer_training_data")

# View top and bottom of data
head(training_data)
tail(training_data)

# Load predictor variable data
temp <- raster("lab_data/DM_primer_predictor_map")

# Plot predictor map (temperature)
plot(temp)

# Plot observations on predictor map
presence_points <- training_data[training_data$presence==1, c("x","y")]
points(presence_points)
absence_points <- training_data[training_data$presence==0, c("x","y")]
points(absence_points,pch="+",col="red")

# Fit a first regression model to training data using GLM
m1 <- glm(presence~temp,family=binomial,data=training_data)

# Inspect model - what are the null and alternative hypotheses?
summary(m1)

# Make map of predictions - how does the species respond?
p1 <- predict(temp,m1,type="response")
plot(p1)

# Put plots on the same page, with headers and different colors for predicted probabilities of presence
par(mfrow=c(3,1),mar=c(2,2,4,5)+0.1) # make space for 3 rows of plots, specify plot margins
plot(temp,main="Temperature")
points(presence_points)
points(absence_points,pch="+",col="red")
plot(p1,col=tim.colors(64),main="Probability of presence - model 1")

# A better model, plotted below the first one
m2 <- glm(presence~temp+I(temp^2),family=binomial,data=training_data) # I() makes R treat temp^2 as a polynomial rather than an interaction term
p2 <- predict(temp,m2,type="response")
plot(p2,col=tim.colors(64),main="Probability of presence - model 2")

# Inspect model 2 and compare to model 1
summary(m2)

# Interpret and discuss the models. The species is Picea sitchensis (Sitka spruce, sitkagran).
# How does the species respond to the predictors?
# How do model predictions compare to the data?
# Search for Picea sitchensis on gbif.org. How do model predictions compare to all observations in GBIF?
# Are responses and predictions in accordance with ecological expectations for this species?
# What could be done to improve the models and predictions?


# R basics
#---------

# Confident users may skip this part of the lab, but please look through the script on your own.
# The following is a brief intro to R syntax, data objects and functions.
# Look at output in the console window and make sure you understand the effect of each command.
# Annotating lines with your own comments may be a good idea.
# Keep an eye on files and objects windows (right) as you go.
# Find more R basics and advanced stuff in RStudio Cheatsheets: https://www.rstudio.com/resources/cheatsheets/

# Command line exercises
1+1
2*2
3^2

# Create and inspect data objects: vectors
x <- 1:4
x
y <- c(5,7,3,4)
y
z <- seq(1,10,by=3)
z

# Access subsets of data with []
z[3]
z[3:4]
z[c(1,4)]
z[-1]
z[-c(1,4)]
z[z>5]

# Modify subsets of data with [] <-
z[1] <- NA
z
z[1:2] <- 5
z

# Create data objects: data.frames, lists
dat <- data.frame(x,y,z)
dat
dat <- data.frame(x_coord=x,y_coord=y,elevation=z)
dat
dat <- data.frame(x_coord=x,y_coord=y,elevation=z,ID=c("a","b"))
dat
dat.list <- list(x_coord=x,y_coord=y,elevation=z,ID=c("a","b"),dataset=dat)
dat.list

# Access and subset data with $, [], [[]]
dat$y_coord
dat[,2]
dat[,"y_coord"]
dat[,c("y_coord","elevation")]
dat[1,c("y_coord","elevation")]
dat[1:2,c("y_coord","elevation")]
dat[c(1,4),c("y_coord","elevation")]
dat[dat$y_coord>4,c("y_coord","elevation")]
dat.list$y_coord
dat.list[[2]]
dat.list[c(2,3)]

# Functions
print("Hello!")
exp(log(2))
mean(c(2,4))
f <- function(x) 2*x+5
f(1)
f(1:3)
f(dat$x_coord)

# Exporting (to e.g. QGIS) and importing
write.csv(dat,'datapoints.csv')
rm(dat)
dat <- read.csv('datapoints.csv')

# Regression
lm(elevation~y_coord,data=dat)
m <- lm(elevation~y_coord,data=dat)
summary(m)
m2 <- glm(elevation~y_coord,data=dat)
summary(m2)

# Plots
plot(x,y)
plot(dat)
plot(m)

# Packages: installing and loading (Packages tab)
install.packages("MIAmaxent")
library("MIAmaxent")

# Help: help files, examples, vignettes
?glm
example(glm)
?MIAmaxent
vignette(package="MIAmaxent")
vignette("a-modeling-example",package="MIAmaxent")


# Other stuff
# install.packages('installr')
# library('installr')
# updateR()
getwd()
setwd()
