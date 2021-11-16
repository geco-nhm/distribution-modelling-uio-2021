####################################
#   HMSC lab, BIOS5211/9211 2021   #
#           EXERCISE 2             #
####################################
#by Ryan Burner, Luke Powell, Lukas Drag, and Eva Lieungh
#for questions after the lab contact Ryan at ryan.burner@outlook.com or Eva at eva.lieungh@nhm.uio.no
#most of this code is from Otso Ovakainen and Nerea Abrego's 2020 book:
# https://www.cambridge.org/us/academic/subjects/life-sciences/ecology-and-conservation/joint-species-distribution-modelling-applications-r?format=PB&isbn=9781108716789
# For more information and useful links, including course materials from a 2020 5-day course, see here: https://www2.helsinki.fi/en/researchgroups/statistical-ecology/hmsc 
####################################
#Data reference: 
#Miller, Jesse E.D., Damschen, Ellen I., Ives, Anthony R. (2018), Data from: Functional traits and community composition: a comparison among community-weighted means, weighted correlations, and multilevel models, Dryad, Dataset, https://doi.org/10.5061/dryad.7gj0s3b
####################################

# Before starting on this script, you should have gone through the 'Lab7_prep.R' script and Exercise 1 in Lab7_HMSC_1.R! If not, you need to load some libraries to continue.
# library(ggplot2)
# library(gridExtra)
# library(lme4)
# library(MASS)
# library(Hmsc)
# library(abind)
# library(ape)
# library(corrplot)
# library(plot.matrix)
# library(Rcpp)
# library(Rlab)
# library(usdm)

# You may want to clean up your R environment before starting on new data.
rm(list = ls())


#############################
#       Exercise 2          #
#   Plant case study w/     #
# multiple species, traits, #
#    and a phylogeny        #
#############################
#Data on plant occurrences documented in revisits to a low-elevation, non-serpentine subset of Robert Whittaker's historic vegetation plots in the Siskiyou Mountains.

data = read.csv('Data/whittaker_revisit_data.csv')
head(data)
str(data)

#now we need to do some formatting to get one row per site in X and Y...
#one column per species in Y...
#factors and factors
data$site = factor(data$site)
sites = levels(data$site)
data$species=factor(data$species)
species = levels(data$species)
n = length(sites)
ns = length(species)
Y = matrix(NA, nrow = n, ncol = ns)
env = rep(NA, n)
trait = rep(NA, ns)
#Question! What will the following for loop do?
for (i in 1:n){
  for (j in 1:ns){
    row = data$site==sites[i] & data$species==species[j]
    Y[i,j] = data[row,]$value
    env[i] = data[row,]$env
    trait[j] = data[row,]$trait
  }
}
colnames(Y) = species
XData = data.frame(TMG = env)
TrData = data.frame(CN = trait)
rownames(TrData)=colnames(Y)

#the env covariate is TMG (topographic moisture gradient)
#the trait is CN (Carbon:Nitrogen ratio, low C:N gives faster growth but less stress tolerance and vice versa)
#so we'd predict high CN to be associated with hotter and drier (high TMG)
head(XData)
Y[1:5,1:5] #head() doesn't work well on Y when so many columns

#let's look at the distribution of prevalence (occupancy) and abundance by species
Prev = colMeans(Y > 0)
Abund = colSums(Y)/colSums(Y > 0)
par(mfrow = c(1,2))
hist(Prev);hist(Abund) #most species are rare; this is very typical

#we want to add phylogeny to the model.
#we don't have a phylo (data class) tree, but we can use taxonomy to build a tree. (Another option could be to use timetree.org, where you can upload a species list and get a phylogenetic tree back! Eva has some code for formatting it to a Hmsc-readable tree if you want to try this)
taxonomy = read.csv('Data/plant_taxonomy.csv',stringsAsFactors = TRUE)
head(taxonomy)
plant.tree = as.phylo(~family/genus/species,
                      data = taxonomy, collapse = FALSE)
plant.tree$edge.length = rep(1, length(plant.tree$edge))
#this method assumes all equal branch lengths and so underestimates close relationships...
#..and overestimates distant ones. But it is something at least
plot(plant.tree, cex = 0.5)

#these are count data, so we'll use a lognormal poisson
#but, we'll also use a probit (pres/abs) model to compare.
#First we specify the formulas for the environmental covariate and trait
XFormula = ~TMG 
TrFormula = ~CN 
models = list()
#Then we specify two models with different distributions.
#for the probit model, we need to transform the Y matrix to 1s and 0s
Ypres=1*(Y > 0)
models[[1]] = Hmsc(Y=Ypres, 
                   XData = XData, XFormula = XFormula,
                   TrData = TrData, TrFormula = TrFormula, 
                   phyloTree = plant.tree,
                   distr = "probit")
models[[2]] = Hmsc(Y = Y, 
                   XData = XData, XFormula = XFormula,
                   TrData = TrData, TrFormula = TrFormula, 
                   phyloTree = plant.tree,
                   distr = "lognormal poisson")

models


#we'll use shorter chains on these more complex models for the sake of time.
#in real life you'd do the opposite (more complex = needs longer chains)!
#set mcmc settings for quick sampling
nChains = 2
thin = 3
samples = 200
transient = 500*thin
verbose = 500*thin

#fit models
for (i in 1:2){
  models[[i]] = sampleMcmc(models[[i]], thin = thin,
                           samples = samples, transient = transient,
                           nChains = nChains, verbose = verbose)
}

models

#Question! What do you want to check now to evaluate the models? (Spoilers below!)



#######################
#  model evaluation   #
#######################

#look at the model fit
preds = computePredictedValues(models[[1]])
MF = evaluateModelFit(hM = models[[1]], predY = preds)
MF$TjurR2
#Question! Why are there so many values? Why are they so small?

#Because there are many species, checking convergence looks a bit different than in last exercise.
mpost = convertToCodaObject(models[[1]], #note the [[1]] to select the model from the list
                            spNamesNumbers = c(TRUE,FALSE), 
                            covNamesNumbers = c(TRUE,FALSE))
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta,multivariate = FALSE)$psrf
#let's look at ESS and psrf of betas...
par(mfrow=c(1,2))
hist(ess.beta,30) #effective sample size should be close to full sample size ( e.g. 400 when 2 chains each with 200 samples)
hist(psrf.beta,30) # remember, psrf is what we would have to multiply the smaller chain by to make it equal the bigger one - should be close to 1
#Question! what do you think about the convergence of this model's mcmc sampling?
#Task! Run these same diagnostics for the second model (lognormal poisson) and compare the plots. Are they different?

#Task! We have fitted the same models with more samples and higher 'thin' value (~45min runtime). Read them in and run the same diagnostics again.
models=readRDS('ex2models.rds')

#Question! How many samples were taken for these models?



################################
#  ecological interpretation   #
################################

# Environment and trait responses 

#Is presence and abundance of most species correlated with the Topographic Moisture Gradient?
#get and plot the estimated Betas
#remember, Betas are the regression coefficients
postBeta = getPostEstimate(models[[1]], parName = 'Beta')
plotBeta(models[[1]], post = postBeta, param = 'Sign', 
         plotTree = TRUE, supportLevel = 0.95, split = 0.4,
         spNamesNumbers = c(FALSE,FALSE))
#Question! How do you interpret the figure?
postBeta = getPostEstimate(models[[2]], parName = 'Beta')
plotBeta(models[[2]], post = postBeta, param = 'Sign', 
         plotTree = TRUE, supportLevel = 0.95, split = 0.4,
         spNamesNumbers = c(FALSE,FALSE))
#Question! How and why are the plots different/similar between models 1 and 2?
#Question! What happens if you replace 'Sign' with 'Mean' in the plotting code above?
#Question! What environmental conditions do most species like?



#Are species with higher Carbon:Nitrogen ratios likely to respond differently to the moisture gradient TMG?
#here we are looking at the Gamma parameters for environmental covariates
postGamma = getPostEstimate(models[[1]], parName = 'Gamma')
plotGamma(models[[1]], post = postGamma, param = 'Mean', #or 'Mean'
          supportLevel = 0.95)
#Question! How do you interpret this figure?
#Task! Plot the Gammas for the second model as well. What are differences/similarities?



#let's plot gradients using each model
#the first plot for the probit model shows species richness across the TMG gradient
#the first plot for the lognormal Poisson model shows total plant abundance across the TMG gradient
#the second plot for both models shows mean community CN trait value across the TMG gradient
par(mfrow=c(2,2))
for (i in 1:2){
  m = models [[i]]
  Gradient = constructGradient(m, focalVariable = 'TMG')
  predY = predict(m, Gradient = Gradient, expected = TRUE)
  q = c(0.25,0.5,0.75)
  plotGradient(m, Gradient, pred = predY, measure = 'S',showData = TRUE, q = q)
  plotGradient(m, Gradient, pred = predY, measure = 'T', index = 2,showData = TRUE, q = q)
}

#ok, let's explore the plotgradient function a bit more
?plotGradient
#we have several options
#measure = 'S' plots species richness (for a probit model), or total count (logn Poisson) or summed response (Normal model)
#it does this by summing the columns of the predicted Y matrix and plotting those sums
par(mfrow=c(1,4))
plotGradient(m, Gradient, pred = predY, measure = 'S',showData = TRUE, q = q, showPosteriorSupport = FALSE )
# measure = 'T' plots the mean trait value in the community across the gradient
# index = <number> shows the number of the trait to plot, in this order (incl. intercept):
m$trNames
plotGradient(m, Gradient, pred = predY, measure = 'T',index = 2, showData = TRUE, q = q, showPosteriorSupport = FALSE)
# measure = 'Y' shows individual species responses
# index = <number> shows which species, following their order in the Y matrix
# let's plot two tree species, one evergreen (Abies concolor) and one deciduous (Acer macrophyllum)
colnames(Y)[1:10]
plotGradient(m, Gradient, pred = predY, measure = 'Y',index = 1, showData = TRUE, q = q, showPosteriorSupport = FALSE)
plotGradient(m, Gradient, pred = predY, measure = 'Y',index = 3, showData = TRUE, q = q, showPosteriorSupport = FALSE)
#NOTE: these plots don't fit the data very well....not sure what is going on here. Any ideas?

#############################################

#what about phylogenetic signal in Betas?
#remember that Hmsc looks only for 'residual' phylogenetic signal (i.e. the remaining signal not explained by the traits -- if you want to estimate 'total' phylo signal in niches, run a model without traits)
mpost = convertToCodaObject(models[[1]])
summary(mpost$Rho)$quant
plot(mpost$Rho)
#model 2
mpost = convertToCodaObject(models[[2]])
summary(mpost$Rho)$quant
plot(mpost$Rho)
#Is there a phylogenetic signal in this data set?

#we can also do Variance Partitioning to see what explains the most of the variation in the data
VP = computeVariancePartitioning(models [[1]],
                                 group = c(1,1),
                                 groupnames = 'TMG')
#VP$R2T shows how much variance in Betas is explained by the traits ($Beta)...
#and how much variance in the presence/abundance is explained by traits ($Y)
VP$R2T #quite a bit!

#now for model 2
VP = computeVariancePartitioning(models [[2]],
                                 group = c(1,1), groupnames = 'TMG')
VP$R2T


######################################################
#Congratulations, you made it to the end of the code!
#Did you know that most people have more arms than average?

#If you want more exercises (with spatial & temporal random effects, hurdle models and more), you can choose:
# Ryan's more advanced exercises (get files from Eva), 
# or one of the "Executable scripts" on the HMSC web page: https://www2.helsinki.fi/en/researchgroups/statistical-ecology/hmsc
