########################
#  WINDOWS USERS
########################
#Install Rtools
   ##YOU MAY NEED ADMINISTRATOR RIGHTS TO DO THIS
#use this link:  
   #  https://cran.r-project.org/bin/windows/Rtools/ 

#After install Rtools, you need to tell R where to find it:
write('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', file = "~/.Renviron", append = TRUE)

#Confirm it worked:
Sys.which("make")
#should get something like this:
## "C:\\rtools40\\usr\\bin\\make.exe"

###########################
# Mac and Linux users
###########################
#follow instructions here
# https://github.com/geco-nhm/distribution-modelling-uio-2021/blob/main/lab_exercises/Lab8_occupancy_models/R_scripts/Nimble_Occupancy.R

###########################
# Everyone
###########################

#now install any packages from below that you don't already have
#try to load them with the code below, then install any that are missing and reload them
install.packages(" PUT PACKAGE NAME HERE, ONE AT A TIME ")

#load packages - make sure you have all of these
library(nimble)
library(coda)
library(mcmcplots)
library(wiqid)
library(unmarked)
library(raster)


code <- nimbleCode({
  mu ~ dnorm(0, sd = 1000)
  sigma ~ dunif(0, 1000)
  for(i in 1:10) {
    x[i] ~ dnorm(mu, sd = sigma)
  }
})

data <- list(x = c(2, 5, 3, 4, 1, 0, 1, 3, 5, 3))
inits <- function() list(mu = rnorm(1,0,1), sigma = runif(1,0,10))
#this next step takes a few minutes...
mcmc.output <- nimbleMCMC(code, data = data, inits = inits,
                          monitors = c("mu", "sigma"), thin = 10,
                          niter = 20000, nburnin = 1000, nchains = 3,
                          summary = TRUE, WAIC = TRUE)
mcmc.output
