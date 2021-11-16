##############################
# Prepare for HMSC exercises #
#     BIOS5211/9211 2021     #
##############################

## by Ryan Burner ryan.burner@outlook.com and Eva Lieungh
# If you make it through this code your computer and R environment should be ready for the exercises!
#Disclaimer: I am a Windows user; if you use Mac or other operating system we will need Google to troubleshoot system problems.

#Step 1
# install the Hmsc package. If there are issues with this version, or if you want to be fancy, you can download the probably more updated version from GitHub. Instructions at the bottom of this script!
install.packages("Hmsc")

#Step 2
#install a few other packages that we'll use (this may be an incomplete list)
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("lme4")
install.packages("MASS")
install.packages("Rcpp")
install.packages(c("abind", "ape", "corrplot", "plot.matrix", "Rlab", "usdm"))

#Step 3
#now make sure that all of the packages will load without errors (warnings are fine)
library(ggplot2)
library(gridExtra)
library(lme4)
library(MASS)
library(Hmsc)
library(abind)
library(ape)
library(corrplot)
library(plot.matrix)
library(Rcpp)
library(Rlab)
library(usdm)

#Step 4
#now figure out where your working directory is so you know where to put data files that we want to read in
getwd() # Is this your course/lab folder, containing the Lab 7 scripts and data?
#setwd('C:/.../yourcoursefolder') # Set the working directory to the correct folder if necessary!

#Step 5
#now let's simulate some data and run a bit of Hmsc test code to confirm that things are working.
library(Hmsc)
{
n = 50
x = rnorm(n)
beta1 = 0 #note: we are setting the intercept here for the simulated data
beta2 = 1 #note: we are setting the coefficient for the X covariate here
L = beta1 + beta2*x
y1 = L + rnorm(n, sd = 1)
Y = as.matrix(y1)
XData = data.frame(x = x)
m.normal = Hmsc(Y = Y, XData = XData, XFormula = ~ x)
nChains = 2
thin = 5
samples = 100
transient = 50
verbose = 150
m.normal = sampleMcmc(m.normal, thin = thin, samples = samples,
                      transient = transient, nChains = nChains,
                      verbose = verbose)
m.normal
}
#output should look like this:
#  "Hmsc object with 50 sampling units, 1 species, 2 covariates, 1 traits and 0 random levels
#  Posterior MCMC sampling with 2 chains each with 100 samples, thin 5 and transient 50" 

#if so, then you are all set! If not, ask your lab mate or teachers

#Step 6
#lastly, clean up the test objects
rm(list = ls())
#now you can continue to the first exercise in Lab7_2_HMSC_ex1.R! 


#################################################################

# instructions for Hmsc installation from GitHub

# first install Rtools from here: (this let's you install packages from other sources)
# https://cran.r-project.org/bin/windows/Rtools/
#Follow instructions for installation; can be a bit tricky to get the path sorted out so R knows where to find it, see link

# then install and load devtools package (to install Hmsc from github)
install.packages("devtools")
library(devtools)

# install Hmsc from github (see info here: https://github.com/hmsc-r/HMSC )
install_github("hmsc-r/HMSC")
##NOTE: sometimes installing Hmsc requires updating many packages (some of which may be loaded in your workspace)
#If this throws an error, try exiting R studio and opening a stand-alone R window, then running the install code above (and below)


