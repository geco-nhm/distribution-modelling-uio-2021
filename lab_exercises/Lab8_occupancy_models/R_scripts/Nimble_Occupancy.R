# This file contains code adapted from the file R_BUGS_code_AHM_Vol_1_20170519.R
# available at https://www.mbr-pwrc.usgs.gov/pubanalysis/keryroylebook/, with
# permission of the authors.  It may have been modified so that examples originally
# written to work with JAGS or WinBUGS will work with NIMBLE.  More information
# about NIMBLE can be found at https://r-nimble.org, https://github.com/nimble-dev/nimble,
# and https://cran.r-project.org/web/packages/nimble/index.html.
#
# =========================================================================
#
#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#
#   Marc Kéry & J. Andy Royle
#
#   Created 2 Dec 2015 based on draft from 21 Oct 2015
#
# =========================================================================
#
#This code has been modified by Ryan Burner from code originally created by Jacob Levine and Perry de Valpine
#https://github.com/nimble-training
#
# Chapter 10.3 Simulation and analysis of the simplest possible site-occupancy model
# --------------------------------------------------------------------------

#load packages
library(nimble)
library(coda)
library(mcmcplots)
library(wiqid)
library(unmarked)
library(raster)


##############################
#
# MODEL 1 - SIMPLE OCCUPANCY MODEL WITH SIMULATED DATA
#
##############################

# Simulate some simple occupancy modeling data
# Choose sample sizes and prepare observed data array y
set.seed(24)                  # So we all get same data set
nSites <- 100  # Number of sites
nVisits <- 2    # Number of presence/absence measurements (visits/surveys) per site
y <- matrix(NA, nrow = nSites, ncol = nVisits) # to contain the obs. data

# Set real Parameter values for simulation
psi <- 0.8  # Probability of occupancy or presence
p <- 0.5    # Probability of detection

# Generate presence/absence data (the truth), for M sites with probability psi
z <- rbinom(n = nSites, size = 1, prob = psi)  # R has no Bernoulli

# Generate detection/non-detection data (i.e. presence/absence measurements)
# for each of two visits (columns)
for(v in 1:nVisits){
  y[,v] <- rbinom(n = nSites, size = 1, prob = z*p)
}

# Look at data
sum(z)          # True number of occupied sites
sum(apply(y, 1, max)) # Observed number of occupied sites
head(cbind(z=z, y))   # Truth, and then 2 measurements (visits/surveys), for first 6 sites

#how would we make a rough guess at p and psi from these data?
   #the model will do basically the same thing, but in a more rigorous way
#let's look at rows where the species was detected at any time
head(y[rowSums(y)>0,])
#let's see in what % of sites the species was ever detected
length(which(rowSums(y)>0))/nSites #'naive' estimate of psi
psi #true value

#now let's get an average probability of a detection in a survey
mean(rowSums(y[rowSums(y)>0,1:nVisits])/nVisits) #simple estimate of p
p #true value

######################
### QUESTION 1 #####
######################
##IS A psi ESTIMATE CALCULATED THIS WAY BIASED IN A PARTICULAR DIRECTION? WHAT ABOUT p?
  #Type your thoughts here:
    # https://docs.google.com/spreadsheets/d/176MnBGQafdaozYHhel5E-2H-yDpXh7HYTJzwelvsG_g/edit?usp=sharing 

######################
### QUESTION 2 #####
######################
##THESE EASY CALCULATIONS DO WORK REASONABLY WELL; WHAT IS THE USE OF THE MODEL THEN?
#Type your thoughts here:
# https://docs.google.com/spreadsheets/d/176MnBGQafdaozYHhel5E-2H-yDpXh7HYTJzwelvsG_g/edit?usp=sharing 


#now back to the nimble model....
# Bundle data and summarize data bundle
# nimble wants each input data object to be a named matrix, vector, or value in a list 
nim.data <- list(y = y)
nim.constants <- list(nS = nrow(y), nV = ncol(y))
nim.data

# Specify model in BUGS language
Section10p3_code <- nimbleCode({

    # Priors
  psi ~ dunif(0, 1)
  p ~ dunif(0, 1)
  
    # Likelihood
  for (i in 1:nS) {    # Loop over sites
    z[i] ~ dbern(psi)         # State/process model (the truth)
    for (j in 1:nV) { # Loop over replicate surveys
      y[i,j] ~ dbern(z[i]*p)  # Observation model (the detection)
    }
  }
})


# Supply initial values for MCMC for each parameter and data value that you will estimate
zst <- apply(y, 1, max)       # Ensure that z starts with a '1', whenever there were 1 or more 1 obs
# Initial values are also supplied as a list of named initial values
inits <- list(z = zst,p=0.5,psi=0.5)

# Parameters to be monitored (for mcmc to keep results)
params <- c("psi", "p")

# MCMC settings
ni <- 5000 # total iterations
nt <- 1    # thinning
nb <- 1000# burn-in
nc <- 3  # chains

# Call nimble and summarize posteriors
# This compiles (to C++, for speed) and runs the model all in one step
fm2 <- nimbleMCMC(
  code = Section10p3_code, # model object
  data = nim.data,  # data values
  constants = nim.constants,    # constants - index values and fixed values (e.g. number of sites) 
  inits = inits,      # initial values
  monitors = params,  # parameters for which to save mcmc chains 
  nburnin = nb,      # burn-in period
  nchains = nc,      # chains
  niter = ni,        # total iterations
  samplesAsCodaMCMC = TRUE
)

#output is an 'mcmc.list' object
class(fm2)

#check effective sample sizes of all parameters
bet=effectiveSize(fm2)
bet

## HOW DOES THIS COMPARE TO ACTUAL NUMBER OF SAMPLES?
## WHAT DOES THIS TELL US ABOUT MCMC?
(ni-nb)*nt*nc 

#check difference between chains for all parameters (potential scale reduction factor, psrf)
gel=gelman.diag(fm2,multivariate=FALSE)$psrf
gel

## HOW DOES THIS COMPARE TO THE IDEAL VALUES?

#see structure of chains
head(fm2)

#with few parameters, can use base R plot
plot(fm2)

#with more, try this
mcmcplot(fm2)

## IF YOU ONLY SAW THESE PLOTS WHAT WOULD YOU GUESS ABOUT YOUR ESS? PSRF? WHY?
## HOW WOULD YOU USE INFORMATION IN THESE PLOTS TO ADJUST YOUR MCMC? 

#let's see how our parameter estimates compare to the simulated values
p   # used for simulation
psi # used for simulation

res=round(summary(fm2)$quantiles,2) #from mcmc samples
res

#values for exercise 1....
#difference between psi median estimate and true psi
abs(res['psi','50%']-psi)
#inter-quartile range (IQR) of psi estimate
res['psi','75%']-res['psi','25%']


######################
### EXERCISE 1 #####
######################
#
# What is the lowest p that will allow you to get reasonably precise psi estimates?
#
# How does that change if you add an extra site visit (or two)?
# INSTRUCTIONS:
# Each person can pick corresponding boxes in the tables to fill in here:
# https://docs.google.com/spreadsheets/d/176MnBGQafdaozYHhel5E-2H-yDpXh7HYTJzwelvsG_g/edit?usp=sharing 
# Mark the boxes with your initials so we know that someone is working on each cell
# Then simulate data and run model
#
# [You could also simulate the effect of adding additional sites vs. surveys]
######################








##############################
#
# MODEL 2 - OCCUPANCY MODEL WITH SWISS BIRD SURVEY DATA (Willow Tit) AND COVARIATES
#
##############################

#this is an example modified from these tutorials: 
# https://bcss.org.my/tut/bayes-with-jags-a-tutorial-for-wildlife-researchers/occupancy-modelling/occupancy-with-site-covariates/
# https://bcss.org.my/tut/bayes-with-jags-a-tutorial-for-wildlife-researchers/occupancy-modelling/occupancy-with-survey-covariates/

#read in data from Willow Tit (Parus montanus) surveys in Switzerland
wt <- read.csv("http://bcss.org.my/data/willowtits.csv",comment="#")
head(wt)

#y.1 - y.3 are three surveys (pres/abs)
#rows are sites
#elev[ation] and forest [% cover] are covariates that may affect psi (presence)
#dur[ation of surveys], day [of the year], and length [of transects] may affect p (detection)

# Extract the detection data and convert to a matrix:
Y <- as.matrix(wt[, 1:3])
( n <- rowSums(!is.na(Y)) ) #this shows how many surveys per site

# Standardize the site covariates as before
elevS <- standardize(wt$elev) #elevation
elev2S <- standardize(wt$elev)^2 #elevation squared
forestS <- standardize(wt$forest) #forest cover

# Extract the survey covariates as matrices and standardize
day <- as.matrix(wt[, 12:14])  #survey day of the year
day2=day^2
dayS <- standardize(day)
day2S <- standardize(day2)
head(dayS)

dur <- as.matrix(wt[, 9:11]) #survey duration
durS <- standardize(dur)
head(durS)

# Replace missing covariate values with zero
dayS[is.na(dayS)] <- 0
day2S[is.na(day2S)] <- 0
durS[is.na(durS)] <- 0

#here we will use a uniform prior for the regression coefficients (could also use normal)
#for intercepts, p and psi must be between zero and one (probabilities)
#we can achieve this using a beta distribution
hist(rbeta(1000,1,1)) #here is a beta distribution
#but, because of the logit link function used for p and psi in the model...
#...we take a logit of that beta distribution for the intercept prior:
hist(logit(rbeta(1000,1,1))) #here is the logit of a beta distribution

# Write the model for NIMBLE
wt_surveyCovs.code <- nimbleCode({
  
  # Priors  
  psi0 ~ dbeta(1,1) #feeds into the next line...
  b0 <- logit(psi0) #...to make intercept prior
                    ### ADD FOREST BETA PRIOR ON THIS LINE
  bElev ~ dunif(-5, 5)   # elevation
  bElev2 ~ dunif(-5, 5)  # elevation^2
  
  p0 ~ dbeta(1, 1) #feeds into the next line...
  a0 <- logit(p0)  #...to make intercept prior
  aDur ~ dunif(-5, 5)    # duration
  aDay ~ dunif(-5, 5)    # day
  aDay2 ~ dunif(-5, 5)   # day^2
  
  # Likelihood
  for(i in 1:nSites) {
    # biological model  
    logit(psi[i]) <- b0 + #intercept
                          # ADD FOREST COVER BETA PARAMETER AND COVARIATE ON THIS LINE
      (bElev * ele[i]) +  #elevation
      (bElev2 * ele2[i]) #elevation squared
    
    z[i] ~ dbern(psi[i])  #here we estimate occupancy at each site
    
    # observation model
    for(j in 1:n[i]) {
      logit(p[i, j]) <- a0 + #intercept
        (aDur * dur[i, j]) + #duration
        (aDay * day[i, j]) +  #day
        (aDay2 * day2[i, j]) #day squared
      
      Y[i, j] ~ dbern(p[i, j] * z[i]) #here we estimate detection in each survey
    }
  }
  
  
  # Derived variable - Nimble lets us derive any variables we want
  #(easy way to capture posterior distributions of the derived quantities in the output)
  N <- sum(z[1:nSites]) #total number of sites occupied
})

# Prepare stuff for NIMBLE
y <- rowSums(Y, na.rm=TRUE)
z <- ifelse(y > 0, 1, NA) #initial values for z, so no sites with defections have a 0
nim.data <- list(Y = Y, z = z) # we won't bother with initial values for betas, just z
str(nim.data)
nim.const <- list(n = n, nSites = length(n),
                 fst = forestS, ele = elevS, ele2 = elev2S,
                 dur = durS, day = dayS, day2 = day2S)
str(nim.const)

#parameters to monitor ### ADD FOREST BETA COVARIATE TO THIS VECTOR ###
params <- c("p0", "a0", "aDur", "aDay", "aDay2",
            "psi0", "b0", "bElev",        "bElev2", "N") #ADD FOREST COVER BETA IN THIS GAP  

# MCMC settings
ni <- 10000 # total iterations
nt <- 2    # thinning
nb <- 1000# burn-in
nc <- 3  # chains

# Call nimble and summarize posteriors
# This compiles (to C++, for speed) and runs the model all in one step
# you'll get some warnings (esp about NAs) but not errors
fm3 <- nimbleMCMC(
  code = wt_surveyCovs.code, # model object
  data = nim.data,  # data values
  constants = nim.const,    # constants - index values and fixed values (e.g. number of sites) 
  monitors = params,  # parameters for which to save mcmc chains 
  nburnin = nb,      # burn-in period
  nchains = nc,      # chains
  niter = ni,        # total iterations
  samplesAsCodaMCMC = TRUE
)

#

#check effective sample sizes of all parameters
bet=effectiveSize(fm3)
hist(bet,50) #more parameters than last model; we'll look at histogram

#check difference between chains for all parameters (potential scale reduction factor, psrf)
gel=gelman.diag(fm3,multivariate=FALSE)$psrf
hist(gel,30)

#with more, try this
mcmcplot(fm3)

#view summary of parameter estimates
round(summary(fm3)$quantiles,3)

#let's plot estimates of p along a day of detection gradient
# Plot effect of day on detection at the average visit duration

#get series from min to max day
xx <- seq(min(wt$day.1,na.rm=TRUE), max(wt$day.3,na.rm=TRUE))
xx

# Need to scale xx the same way we scaled the day data
xxS <- standardize2match(xx, day)

# Do a plot with credible interval
a0 <- fm3$a0 #get mcmc chains for each parameter
aDay <- fm3$aDay
aDay2 <- fm3$aDay2

#now generate data using the model equation and each sample from the mcmc
toplot <- matrix(NA, length(xx), 3)
for(i in 1:length(xx)) { #only use 100 samples for speed
  logit.p.tmp <- a0 + aDay * xxS[i] + aDay2 * xxS[i]^2 #this is the equation for p
  p.tmp <- plogis(logit.p.tmp)
  toplot[i, 1] <- mean(p.tmp) #extract mean
  toplot[i, 2:3] <- hdi(p.tmp) #extract CI
}

plot(xx, toplot[, 1], ylim=range(toplot), type='n', las=1,
     xlab="Day", ylab="Detection probability")
polygon(x=c(xx, rev(xx)), y=c(toplot[, 2], rev(toplot[, 3])),
        col=adjustcolor('skyblue', 0.5), border=NA)
lines(xx, toplot[, 1], lwd=2, col='blue')
# Add intercept
abline(h=mean(fm3$p0))
abline(h=hdi(fm3$p0), lty=2)

#Note that p is estimated the best in the middle of the season
hist(day)





######################
### EXERCISE 2 #####
######################
#
# Make predictions of willow tit distributions across all of Switzerland
#
# INSTRUCTIONS:
# First, you'll need to estimate a covariate ('bFor') for the effect of forest cover
# We've already included the forest cover data ('fst') in the model data
# You just need to include it in three places:
#     1) The model prior for 'bFor' (the beta parameter for forest cover) 
#     2) The model likelihood function for psi - use these names: (bFor * fst[i])
#     3) The 'params' to monitor: 'bFor'
#   HINT: FOR HELP UPDATING THE CODE IN THIS WAY, SEE THE 'code hints' tab in the google sheet:
# https://docs.google.com/spreadsheets/d/176MnBGQafdaozYHhel5E-2H-yDpXh7HYTJzwelvsG_g/edit?usp=sharing 
# After you've rerun the model with forest covariates, paste your map (from below) into the google sheet
#
######################






#Now let's look at the effect of covariates on psi
# Convert to mcmcOutput
( mco1 <- mcmcOutput(fm3) )
summary(mco1)
sumry1 <- sumryList(mco1)
str(sumry1)


# Plot effect of elevation on occupancy at the average forest cover
range(wt$elev)
xx <- seq(250, 2750, length=101)
# Need to scale xx the same way we scaled the elev data
xxS <- standardize2match(xx, wt$elev)

# Do a simple plot using the point estimates of the coefficients
logit.psi <- with(sumry1$mean, b0 + bElev * xxS + bElev2 * xxS^2)
psi <- plogis(logit.psi)
plot(xx, psi, type='l') # small 'L', not 'one'

# Check on forest at high elevation
plot(wt$elev, wt$forest)
abline(h=35)

# Plot occupancy vs elevation with forest = 0, 35 and 80%
( forS <- standardize2match(c(0, 35, 80), wt$forest) )

# Extract the coefficients we need:
b0 <- mco1$b0
bFor <- mco1$bFor
bElev <- mco1$bElev
bElev2 <- mco1$bElev2

toplot0 <- matrix(NA, 101, 3)
for(i in 1:101) {
  logit.psi.tmp <- b0 + bFor * forS[1] + bElev * xxS[i] + bElev2 * xxS[i]^2
  psi.tmp <- plogis(logit.psi.tmp)
  toplot0[i, 1] <- mean(psi.tmp)
  toplot0[i, 2:3] <- hdi(psi.tmp)
}

which(xx == 2200)
toplot35 <- matrix(NA, 79, 3)
for(i in 1:79) {
  logit.psi.tmp <- b0 + bFor * forS[2] + bElev * xxS[i] + bElev2 * xxS[i]^2
  psi.tmp <- plogis(logit.psi.tmp)
  toplot35[i, 1] <- mean(psi.tmp)
  toplot35[i, 2:3] <- hdi(psi.tmp)
}

which(xx == 2000)
toplot80 <- matrix(NA, 71, 3)
for(i in 1:71) {
  logit.psi.tmp <- b0 + bFor * forS[3] + bElev * xxS[i] + bElev2 * xxS[i]^2
  psi.tmp <- plogis(logit.psi.tmp)
  toplot80[i, 1] <- mean(psi.tmp)
  toplot80[i, 2:3] <- hdi(psi.tmp)
}

plot(xx, toplot0[, 1], ylim=0:1, type='n', las=1,
     xlab="Elevation", ylab="Occupancy")
polygon(x=c(xx, rev(xx)), y=c(toplot0[, 2], rev(toplot0[, 3])),
        col=adjustcolor('yellow', 0.5), border=NA)
polygon(x=c(xx[1:79], rev(xx[1:79])), y=c(toplot35[, 2], rev(toplot35[, 3])),
        col=adjustcolor('skyblue', 0.5), border=NA)
polygon(x=c(xx[1:71], rev(xx[1:71])), y=c(toplot80[, 2], rev(toplot80[, 3])),
        col=adjustcolor('lightgreen', 0.5), border=NA)
lines(xx, toplot0[, 1], lwd=2, col='brown')
lines(xx[1:79], toplot35[, 1], lwd=2, col='blue')
lines(xx[1:71], toplot80[, 1], lwd=2, col='darkgreen')
legend('topleft', c("80% forest", "35% forest", "0% forest"),
       lwd=2, col=c('darkgreen', 'blue', 'brown'), bty='n')

# Species distribution map
# ========================
#download forest and elevation raster map of Switzerland
data(Switzerland)
str(Switzerland)

# Convert elevation and forest cover to rasters and plot
elevR <- rasterFromXYZ(Switzerland[,1:3])
forR <- rasterFromXYZ(Switzerland[,c(1,2,4)])
dev.new() #creates a new plot window; use this if it won't plot in your window
par(mfrow=1:2)
plot(elevR, main="Elevation", axes=FALSE, box=FALSE, col=terrain.colors(255))
plot(forR, main="Forest cover", axes=FALSE, box=FALSE)
dev.off() #removes the new plot window
par(mfrow=c(1, 1))

# Scale the elevation and forest covariates as we did for the willow tits data
# ('standardize2match' doesn't work for rasters.)
forRS <- (forR - mean(wt$forest)) / sd(wt$forest)
elevRS <- (elevR - mean(wt$elev)) / sd(wt$elev)

# Calculate occupancy for each pixel of the map
logit.psi <- with(sumry1$mean, b0 + bFor * forRS +
                    bElev * elevRS + bElev2 * elevRS^2)
psiR <- 1 / (1 + exp(-logit.psi))  # 'plogis' doesn't work for rasters

mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(psiR, col=mapPalette(100), main="Willow tit occupancy", axes=FALSE, box=FALSE)

# Overall mean occupancy in Switzerland
mean(values(psiR), na.rm=TRUE)



### For next time
## assessing model fit, and model selection
    #No time today....but don't forget these crucial steps in your own work!!


