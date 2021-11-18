# Dataprep functions

generate.training.data <- function(sp,predictors,generate.absences=FALSE,absence.value=0,n.abs=NA,factors=NA)
# The arguments of this function are:
# sp: species obs returned by rgbif::occ_search
# predictors: predictor raster stack, e.g. generated from Horvath et al. 2019 (see Lab2)
# generate.absences: should random absences be generated? TRUE or FALSE
# absence.value: value for coding of absences, depending on intended use of training data (e.g. 0 for glm (default), NA for MIAmaxent)
# n.abs: maximum number of absences to generate (default NA, which will give the same number of absences as presences, see below)
# factors: names of variables (character vector) to be treated as factors (or unique parts of the variable names)
{
  sp.data <- data.frame(sp$data)

  # Convert lat-long coordinates to coordinate system (crs) of environmental raster data
  occ_points <- data.frame(x=sp.data$decimalLongitude,y=sp.data$decimalLatitude)
  occ_points <- SpatialPoints(occ_points,proj4string=CRS("+proj=longlat +datum=WGS84"))
  occ_UTM33 <- spTransform(occ_points,crs(predictors))
  sp.data$x <- occ_UTM33$x
  sp.data$y <- occ_UTM33$y

  # Point data like 'sp.data' can be used directly as input to MIAmaxent (Labs 4-5),
  # along with environmental rasters.
  # For glm (Lab 3), we need absence data in addition to presences, and we
  # need to combine presences and absences with environmental data in a
  # data set suitable for the glm function. This can be done as follows.
  
  # Make background raster
  bg.raster <- predictors[[1]]/predictors[[1]] # raster with values=1 in all cells with values of predictors[[1]] (e.g. all of mainland Norway)
  
  # Rasterize occurrences
  occ_ras <- rasterize(occ_UTM33,bg.raster,fun='count',background=0) # raster with counts of occurrences in each cell of bg.raster
  occ_ras <- mask(occ_ras,bg.raster)                                 # filtering occurrence raster with bg.raster
  occ_ras[occ_ras>0] <- 1                                            # reducing counts>0 to 1 (presence)
  
  # Take the occurrence cells as presences
  presences <- which(raster::values(occ_ras)==1)

  if(generate.absences)
  {
    # If required, generate absence data by sampling from the cells without occurrence observations.
    # NB! This rests on risky assumptions, but we need the absences
    # for logistic regression with glm (Lab 3). We make this absence data set for
    # educational purposes, but ideally one would want a data set with true, verified presences
    # and absences for logistic regression (and also for validation of 'presence-only' methods,
    # such as maxent).
    absences <- which(raster::values(occ_ras)==0)
    if(is.na(n.abs)) n.abs <- length(presences)
    if(n.abs > length(absences)) n.abs <- length(absences)
    absences_sample <- sample(absences,size=n.abs) # sample of same number of absence cells as presence cells
 
    # Combine presences, absences and environmental data
    selected <- c(presences,absences_sample)
  }
  else selected <- presences
  
  occ_ras[occ_ras==0] <- absence.value   # recoding absences to desired absence value
  
  xy <- raster::coordinates(occ_ras)[selected, ]
  training_data <- data.frame(presence=raster::values(occ_ras)[selected],
                              xy,
                              extract(predictors,xy))

  # Convert discrete environmental predictors to factor variables
  if(!is.na(factors[1]))
  {
    for(i in 1:length(factors))
    {
      j <- grep(factors[i],names(training_data))     # find column index with name containing the i'th factor
      training_data[,j] <- factor(training_data[,j]) # convert to factor
    }
  }
  
  # Return training data
  return(training_data)
}
