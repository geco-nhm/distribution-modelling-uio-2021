####################################
#   HMSC lab, BIOS5211/9211 2021   #
#           EXERCISE 1             #
####################################
#by Ryan Burner, Luke Powell, Lukas Drag, and Eva Lieungh
#for questions after the lab contact Ryan at ryan.burner@outlook.com or Eva at eva.lieungh@nhm.uio.no
#most of this code is from Otso Ovakainen and Nerea Abrego's 2020 book:
# https://www.cambridge.org/us/academic/subjects/life-sciences/ecology-and-conservation/joint-species-distribution-modelling-applications-r?format=PB&isbn=9781108716789
# For more information and useful links, including executable scripts and course materials from a 2020 5-day course, see here: https://www2.helsinki.fi/en/researchgroups/statistical-ecology/hmsc 
################################################################

# NB! Before starting on this script, you should have gone through the 'Lab7_prep.R' script! If not, you need to load some libraries to continue. (Un-comment the next lines with Ctrl+Shift+C !)
# library(ggplot2)
# library(gridExtra)
# library(lme4)
# library(MASS)
# library(Hmsc)
# library(abind)
# library(ape)
# library(corrplot)
# library(plot.matrix)
# library(Rlab)
# library(usdm)



############################
#      EXERCISE 1          #
# Fit and evaluate a basic #
# model with several types #
# of response variables    #
############################

#simulate response data (y)
#response to one covariate (x)
#you can imagine this is the body mass, occurrence and counts of one species, and its response to forest cover
set.seed(1) 
n = 50
x = rnorm(n)
beta1 = 0 #note: we are setting the intercept here for the simulated data
beta2 = 1 #note: we are setting the coefficient for the X covariate here
L = beta1 + beta2*x
y1 = L + rnorm(n, sd = 1) #normal
y2 = 1*((L + rnorm(n, sd = 1)) > 0) #pres/abs
y3 = rpois(n = n, lambda = exp(L + rnorm(n, sd = 1))) #count data

#plot the three response vars next to each other
par(mfrow = c(1,3))
#normal response variable (e.g. body mass)
plot(x,y1, main = "Normal")
#binary response (e.g. presence/absence)
plot(x,y2, main = "Probit")
#integer response (e.g. abundance counts)
plot(x,y3, main = "Lognormal Poisson")
#remember, we can think of the x covariate as forest cover. It has negative values, so think of it as scaled to make sense of the plots.

#fit a frequentist linear model to the normal response variable
df = data.frame(x, y1)
m.lm = lm(y1 ~ x, data = df )
# Question! Are the model intercept and x-coefficient similar to our simulated values? Look at a summary of the linear model:
summary(m.lm)

#now we are going to use the same, simulated data to make a simple Hmsc model. First format the same data for HMSC
Y = as.matrix(y1)
XData = data.frame(x = x)
#define simple Hmsc model
m.normal = Hmsc(Y = Y, 
								XData = XData, 
								XFormula = ~ x
								)

# the model is now defined, but not fitted. The parameters need to be estimated.
#set mcmc settings
nChains = 2 # number of mcmc chains
thin = 5 # number of MCMC steps between each recording of samples from the posterior
samples = 1000 # samples per chain
transient = 500*thin # burn-in length
verbose = 500*thin # frequency of progress reporting (messages appearing in console)

#fit the model
m.normal = sampleMcmc(m.normal, thin = thin, samples = samples,
                      transient = transient, nChains = nChains,
                      verbose = verbose) #while you are waiting, you could watch this video: https://www.youtube.com/watch?v=XV4yj4T4PBQ&ab_channel=pasky007
#Question! Look in the Console. What is the difference between 'transient' and 'sampling'?

#extract mcmc chains for each parameter so we can inspect them
mpost = convertToCodaObject(m.normal)

#summarize results
# Question! Are the model intercept and x-coefficient similar to our simulated values? Look at a summary of the Hmsc model:
summary(mpost$Beta) 
#note that a 'Beta' in Hmsc lingo = a regression coefficient
#remember that in Bayesian models we don't have to calculate 95% CIs from the SD.
#instead, we just look at the mcmc sample distributions

#to calculate model fit, we can compare the model predictions to the actual data
preds = computePredictedValues(m.normal)
MF = evaluateModelFit(hM = m.normal, predY = preds)
#this is our r2 value for the model
MF$R2

# Question! Can you write a plain-language description of the results of these two models (linear and Hmsc)? How does the imaginary species' body mass vary with forest cover?


#Ok, we have some results from two models, but...
#before we trust the HMSC model we need to know if:
# 1) the mcmc chains 'converged' on the same estimate for all parameters
# 2) the mcmc chains were long enough, and not too autocorrelated, to give us enough samples to describe the distribution

#take a look at the chains
plot(mpost$Beta) # you might need to make plot window bigger for this one
#Question! Is the chain convergence OK? Why, or why not? What would you do to improve convergence?

#our sample size is the number of samples in all chains combined
#but are all these samples independent? If they are correlated, we have less information than the number of samples would indicate
#check effective sample size (ESS)
effectiveSize(mpost$Beta) #sample size left, ESS right
#Question! Is there evidence of autocorrelation? How does ESS compare with raw sample size in this Hmsc model? What could you do to improve the ESS?

#To understand if the chains converged on similar values, we calculate the'potential scale reduction faction' (psrf)
#aka 'Gelman' diagnostics, it tells us what we would have to multiply the smaller chain by to make it equal to the bigger one
# 1 = perfect; under 1.05 or less is probably good (in larger models, most under 1.1 probably ok)
gelman.diag(mpost$Beta,multivariate = FALSE)$psrf

#good also to check model assumptions
#here for linear model
nres.lm = rstandard(m.lm) # regression deletion diagnostics
preds.lm = fitted.values(m.lm)
par(mfrow = c(1,2))
hist(nres.lm, las = 1) # should be normally distributed
plot(preds.lm, nres.lm, las = 1) # should not have clear trend
abline(a = 0, b = 0)
#Question! Does it look OK? Why or why not?

#now for the Hmsc model
preds.mean = apply(preds, FUN = mean, MARGIN = 1)
nres = scale(y1-preds.mean)
par(mfrow = c(1,2))
hist(nres)
plot(preds.mean, nres)
abline(a = 0, b = 0)

######################################################
# Progress check! How are you doing on time? If you are approaching halfway through the lab time, maybe skip (or go really fast through) the rest of this exercise and rather start Exercise 2 in Lab7_HMSC_2.R
# The rest of this exercise goes through mostly the same steps, but for the presence-absence and count data types.
######################################################

#now, let's consider different response variable types
#we just fit a model to normally distributed data above (called m.normal)...
#but for presence/absence, you want a probit distribution
# so we fit a new model with y2 data (instead of y1), same as before
Y = as.matrix(y2) #remember we simulated probit data above
m.probit = Hmsc(Y = Y, 
								XData = XData, 
								XFormula = ~ x,
                distr = "probit" # distribution defined here
								)

m.probit = sampleMcmc(m.probit, thin = thin, samples = samples,
                      transient = transient, nChains = nChains, verbose = verbose)

m.probit

#check convergence again
mpost = convertToCodaObject(m.probit)
#Gelman psrf
gelman.diag(mpost$Beta,multivariate = FALSE)$psrf
#ESS (effective sample size)
effectiveSize(mpost$Beta) 

#Question! Why is the Effective Sample Size smaller than in the normal model?

#check estimates
round(summary(mpost$Beta)$quantiles, 2)

#Question! How do the estimates look, compared to the previous two models? Why? 

#now lets look at fit
preds = computePredictedValues(m.probit)
evaluateModelFit(hM = m.probit, predY = preds)
#Look at the object that just printed. It has several values
#RMSE is mean square error (a measure of model fit)
#AUC is area under the curve (ranges from 0 to 1, 0.5 is as good as random model)
#TjurR2 is an r2 value used for pres/abs data (regular r2 isn't a good choice for binary models)

#now let's do all the same for count data
Y = as.matrix(y3)
m.lognormal.poisson = Hmsc(Y = Y, 
													 XData = XData, 
													 XFormula = ~ x,
                           distr = "lognormal poisson"
													 )
m.lognormal.poisson = sampleMcmc(m.lognormal.poisson,
                                 thin = thin, samples = samples,
                                 transient = transient,
                                 nChains = nChains,
                                 verbose = verbose)
m.lognormal.poisson

#Question! The following steps should be looking familiar now. Can you annotate each one and interpret them?
mpost = convertToCodaObject(m.lognormal.poisson)
effectiveSize(mpost$Beta)
gelman.diag(mpost$Beta,multivariate = FALSE)$psrf
round(summary(mpost$Beta)$quantiles,2)
preds = computePredictedValues(m.lognormal.poisson,
                               expected = FALSE)
evaluateModelFit(hM = m.lognormal.poisson, predY = preds)
#now we have even more measures of fit -- briefly, those starting with 'O.' refer to pres/abs prediction
#and those starting with 'C.' refer to numeric count prediction

#now we'll make model predictions across a gradient of Y values for each of the 3 models
#grey points are raw data, ribbons are 95% CIs
#these ribbons provide estimates of the true value of Y at each value of X
par(mfrow = c(1,3))
for (i in 1:3){
  m = switch (i, m.normal,m.probit, m.lognormal.poisson)
  Gradient = constructGradient(m, focalVariable = "x")
  predY = predict(m, Gradient = Gradient, expected = TRUE)
  plotGradient(m, Gradient, pred = predY, measure = "Y",showPosteriorSupport = FALSE,
               index = 1, showData = TRUE, main = c("Normal",
                                                    "Probit",
               																		  "Lognormal Poisson")[i])
}

# Question! How do the modelled responses compare to the data we plotted at the beginning?


######################################################
#Congratulations, you made it to the end of Exercise 1!
#Did you know that sloths are faster in water than on the ground?

#Now have a short break and then do Exercise 2!