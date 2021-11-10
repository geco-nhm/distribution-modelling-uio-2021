## ----setup, include = FALSE----------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ---- fig.show='hold', fig.width=5.5, fig.height=5.5---------------------------------------------------------------------
library(raster)
EV1 <- raster(list.files(system.file("extdata", "EV_continuous", 
                                     package="MIAmaxent"), full.names=TRUE)[1])
PO <- read.csv(system.file("extdata", "occurrence_PO.csv", package="MIAmaxent"))
plot(EV1, legend=FALSE)
points(PO$POINT_X, PO$POINT_Y, pch = 20, cex = 0.5, col = 'blue')


## ------------------------------------------------------------------------------------------------------------------------
library(MIAmaxent)
grasslandPO <- readData(
  occurrence=system.file("extdata", "occurrence_PO.csv", package="MIAmaxent"), 
  contEV=system.file("extdata", "EV_continuous", package="MIAmaxent"),
  catEV=system.file("extdata", "EV_categorical", package="MIAmaxent"),
  maxbkg=20000)


## ------------------------------------------------------------------------------------------------------------------------
str(grasslandPO)
sum(grasslandPO$RV == 1, na.rm = TRUE)
sum(is.na(grasslandPO$RV))


## ---- fig.show='hold', fig.width=5.5, fig.height=4-----------------------------------------------------------------------
teraspifFOP <- plotFOP(grasslandPO, "teraspif")
terslpdgFOP <- plotFOP(grasslandPO, "terslpdg")


## ---- fig.width=5.5, fig.height=4----------------------------------------------------------------------------------------
terslpdgFOP <- plotFOP(grasslandPO, "terslpdg", span = 0.75, intervals = 20)
terslpdgFOP


## ---- fig.width=5.5, fig.height=4----------------------------------------------------------------------------------------
geobergFOP <- plotFOP(grasslandPO, 10)


## ------------------------------------------------------------------------------------------------------------------------
geobergFOP


## ---- warning=FALSE------------------------------------------------------------------------------------------------------
grasslandDVs <- deriveVars(grasslandPO, 
                           transformtype = c("L","M","D","HF","HR","T","B"))


## ------------------------------------------------------------------------------------------------------------------------
summary(grasslandDVs$dvdata)
head(summary(grasslandDVs$transformations))
length(grasslandDVs$transformations)


## ---- fig.show='hold', fig.width=2.75------------------------------------------------------------------------------------
plot(grasslandPO$terslpdg, grasslandDVs$dvdata$terslpdg$terslpdg_D2, pch=20, 
     ylab="terslpdg_D2")
plot(grasslandPO$terslpdg, grasslandDVs$dvdata$terslpdg$terslpdg_M, pch=20,
     ylab="terslpdg_M")


## ---- warning=FALSE------------------------------------------------------------------------------------------------------
grasslandDVselect <- selectDVforEV(grasslandDVs$dvdata, alpha = 0.001, quiet = TRUE)


## ------------------------------------------------------------------------------------------------------------------------
summary(grasslandDVs$dvdata)
sum(sapply(grasslandDVs$dvdata[-1], length))
summary(grasslandDVselect$dvdata)
sum(sapply(grasslandDVselect$dvdata[-1], length))


## ------------------------------------------------------------------------------------------------------------------------
grasslandDVselect$selection$terdem


## ------------------------------------------------------------------------------------------------------------------------
grasslandEVselect <- selectEV(grasslandDVselect$dvdata, alpha = 0.001, 
                              interaction = TRUE)


## ------------------------------------------------------------------------------------------------------------------------
summary(grasslandDVselect$dvdata)
length(grasslandDVselect$dvdata[-1])
summary(grasslandEVselect$dvdata)
length(grasslandEVselect$dvdata[-1])


## ------------------------------------------------------------------------------------------------------------------------
grasslandEVselect$selectedmodel$formula


## ------------------------------------------------------------------------------------------------------------------------
grasslandEVselect$selection[!duplicated(grasslandEVselect$selection$round), ]


## ---- fig.width=4, fig.height=4------------------------------------------------------------------------------------------
plot(grasslandEVselect$selection$round, grasslandEVselect$selection$Dsq, 
     xlab="round", ylab="Dsq")


## ------------------------------------------------------------------------------------------------------------------------
grasslandmodel <- chooseModel(grasslandDVselect$dvdata, 
                              formula("~ prbygall + geoberg + lcucor1 + 
                                      tertpi09 + geolmja1"))


## ---- fig.width=5.5, fig.height=4----------------------------------------------------------------------------------------
plotResp(grasslandmodel, grasslandDVs$transformations, "prbygall")


## ---- fig.show='hold', fig.width=5.5, fig.height=4-----------------------------------------------------------------------
prbygallFOP <- plotFOP(grasslandPO, "prbygall")
plotResp(grasslandmodel, grasslandDVs$transformations, "prbygall")


## ---- fig.show='hold', fig.width=5.5, fig.height=4-----------------------------------------------------------------------
geolmja1FOP <- plotFOP(grasslandPO, "geolmja1")
plotResp(grasslandmodel, grasslandDVs$transformations, "geolmja1")


## ---- fig.width=5.5, fig.height=4----------------------------------------------------------------------------------------
plotResp2(grasslandmodel, grasslandDVs$transformations, "prbygall")


## ------------------------------------------------------------------------------------------------------------------------
calculateFTVA(grasslandEVselect, formula("~ prbygall + geoberg + lcucor1 + 
                                      tertpi09 + geolmja1"))


## ----fig.show='hold', fig.width=5.5, fig.height=5.5----------------------------------------------------------------------
EVfiles <- c(list.files(system.file("extdata", "EV_continuous", package="MIAmaxent"), 
             full.names=TRUE),
             list.files(system.file("extdata", "EV_categorical", package="MIAmaxent"), 
             full.names=TRUE))
EVstack <- raster::stack(EVfiles)
names(EVstack) <- gsub(".asc", "", basename(EVfiles))
grasslandPreds <- projectModel(model = grasslandmodel,
                               transformations = grasslandDVs$transformations,
                               data = EVstack)


## ----fig.show='hold', fig.width=5.5, fig.height=5.5----------------------------------------------------------------------
plot(log2(grasslandPreds$output + 1))


## ------------------------------------------------------------------------------------------------------------------------
grasslandPreds


## ------------------------------------------------------------------------------------------------------------------------
grasslandPA <- readData(
  occurrence = system.file("extdata", "occurrence_PA.csv", package="MIAmaxent"), 
  contEV = system.file("extdata", "EV_continuous", package="MIAmaxent"),
  catEV = system.file("extdata", "EV_categorical", package="MIAmaxent"),
  PA = TRUE, XY = TRUE)
head(grasslandPA)
tail(grasslandPA)


## ---- fig.show='hold', fig.width=5.5, fig.height=5.5---------------------------------------------------------------------
plot(log2(grasslandPreds$output + 1))
presences <- grasslandPA[grasslandPA$RV==1, ]
absences <- grasslandPA[grasslandPA$RV==0, ]
points(presences$x, presences$y, pch = 20, cex = 0.5, col = 'red')
points(absences$x, absences$y, pch = 20, cex = 0.5, col = 'blue')
legend('topleft', c('presence', 'absence'), col = c('red', 'blue'), 
       pch = c(20, 20))


## ---- fig.width=4, fig.height=4------------------------------------------------------------------------------------------
testAUC(model = grasslandmodel, transformations = grasslandDVs$transformations,
        data = grasslandPA)


## ---- fig.show='hold', fig.width=5.5, fig.height=4-----------------------------------------------------------------------
grasslandPA <- readData(
  occurrence = system.file("extdata", "occurrence_PA.csv", package="MIAmaxent"), 
  contEV = system.file("extdata", "EV_continuous", package="MIAmaxent"),
  catEV = system.file("extdata", "EV_categorical", package="MIAmaxent"),
  PA = TRUE, XY = FALSE)
str(grasslandPA)

plotFOP(grasslandPA, "teraspif")
plotFOP(grasslandPA, "terslpdg")
plotFOP(grasslandPA, 10)


## ---- warning = FALSE----------------------------------------------------------------------------------------------------
PA.grasslandDVs <- deriveVars(grasslandPA, algorithm = "LR")


## ---- fig.show='hold', fig.width=5.5, fig.height=4-----------------------------------------------------------------------
PA.grasslandDVselect <- selectDVforEV(PA.grasslandDVs$dvdata, alpha = 0.001, 
                                      algorithm = "LR", quiet = TRUE) 

PA.grasslandEVselect <- selectEV(PA.grasslandDVselect$dvdata, alpha = 0.001, algorithm = "LR")
PA.grasslandEVselect$selection
PA.grasslandEVselect$selectedmodel

plotResp(PA.grasslandEVselect$selectedmodel, PA.grasslandDVs$transformations, "terdem")
plotFOP(grasslandPA, "terdem")


## ---- fig.show='hold', fig.width=5.5, fig.height=4, warning = FALSE------------------------------------------------------
plotResp(PA.grasslandEVselect$selectedmodel, PA.grasslandDVs$transformations, "prbygall")
plotFOP(grasslandPA, "prbygall")


## ------------------------------------------------------------------------------------------------------------------------
calculateFTVA(PA.grasslandEVselect)


## ----fig.show='hold', fig.width=5.5, fig.height=5.5----------------------------------------------------------------------
PA.grasslandPreds <- projectModel(model = PA.grasslandEVselect$selectedmodel,
                               transformations = PA.grasslandDVs$transformations,
                               data = EVstack)


## ------------------------------------------------------------------------------------------------------------------------
# logistic regression model
testAUC(PA.grasslandEVselect$selectedmodel, PA.grasslandDVs$transformations, 
        grasslandPA, plot = FALSE) 
# maxent model
testAUC(grasslandmodel, grasslandDVs$transformations, grasslandPA, plot = FALSE) 


## ------------------------------------------------------------------------------------------------------------------------
# logistic regression model
testAUC(PA.grasslandEVselect$selectedmodel, PA.grasslandDVs$transformations, 
        grasslandPO, plot = FALSE)
# maxent model
testAUC(grasslandmodel, grasslandDVs$transformations, grasslandPO, plot = FALSE)

