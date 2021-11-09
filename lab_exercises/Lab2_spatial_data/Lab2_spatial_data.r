# BIOs 5211/9211
# Olav Skarpaas, Peter Horvath (environmental data) and Dag Endresen (GBIF data)
# NHM, University of Oslo, Nov 2021


# Lab2. Data download, visualization and basic spatial analysis
###############################################################

# In this lab we will go through some tools for working
# with spatial data in R to prepare for distribution modelling.

# Data for Norway
#----------------
# The environmental data set for Norway in this lab is a subset of
# the data compiled and described in Horvath et al. 2019. Distribution 
# modelling of vegetation types based on area frame survey
# data. Applied Vegetation Science, 22(4): 547-560 (see syllabus).

# We will first import and work with the spatial environmental data
# layers, and then combine this with species occurrences to prepare for
# distribution modelling.

# Libraries: spatial data tools
# You may need to install these using 'install.packages()' or the Packages tab
# library(readxl) # excel file tools
library(sp)     # vector tools
library(raster) # raster tools
library(rgdal)  # spatial data tools (geospatial data abstraction library)
library(fields) # color ramps (e.g. tim.colors)

# if you have NOT downloaded the spatial data for this exercise yet, then do so, by following this link: https://uio-my.sharepoint.com/:f:/g/personal/peterhor_uio_no/EjbrdH5bzjVLmWyz7JVVxhkB9vluk4RznDJbxmY54hKVsw?e=0luLHx

# save the files and copy file path
# Set file path for data - modify path to files on your computer
path <- "lab_data/RASTER/" 
# path <- "E:/Project_1_FINAL_layers/TIFF_BIOS5211/"



### Environmental predictor variables: rasters

# Load environmental predictor maps for Norway
# Coordinate Reference System of the raster files is EPSG:32633 - WGS 84 / UTM zone 33N
# Topographic 
DEM <- raster(paste(path,"dem100.tif",sep=""))                    # Digital elevation model (Kartverket)
ASP <- raster(paste(path,"Aspect.tif",sep=""))                    # Aspect (Derived from DEM --- Direction? unit -> radians (?? rad = 180?))
TWI <- raster(paste(path,"Topographic_Wetness_Index.tif",sep="")) # Topographic Wetness Index (Derived from DEM --- unitless - index)
SOL <- raster(paste(path,"Total_insolation.tif",sep=""))          # Potential Incoming Solar Radiation (Derived from DEM --- SAGA, )
# Geological
GEO <- raster(paste(path,"geo_norge123.tif",sep=""))              # Bedrock nutrient status (three classes, SOURCE)
# Climatic
TMP <- raster(paste(path,"bioclim_1.tif",sep=""))                 # Annual mean temperature (Bioclim)
TMP2 <- raster(paste(path,"bioclim_10.tif",sep=""))               # Mean Temperature of Warmest Quarter (Bioclim)
PRC <- raster(paste(path,"bioclim_12.tif",sep=""))                # Annual Precipitation (Bioclim)
GSL <- raster(paste(path,"Growing_season_length.tif",sep=""))     # Growing season length (MET - SeNorge 2, CDO algorithm - derived from temperature)
SW5 <- raster(paste(path,"swe_5.tif",sep=""))                     # Snow water equivalent in May (MET - SeNorge 1.1.1)
SW1 <- raster(paste(path,"swe_1.tif",sep=""))                     # Snow water equivalent in January (MET - SeNorge 1.1.1)
# Land cover
LC    <- raster(paste(path,"ar50_artype.tif",sep=""))             # Land cover classification, based on AR50 (Nibio) 
# 10	Developed area. Residential area, town center, city, transport, industrial area and alike.
# 20	Agricultural area. Cropland, cultivated soil and cultivated pastures
# 30	Forest. Afforested area
# 50	Barren land. Land areas with natural vegetation cover that is not forest.
# 60	Bog and Fen. Land areas with characteristics of marshes
# 70	Glacier. Ice and snow that does not melt during the summer.
# 81	Freshwater. Rivers and lakes
# 82	Ocean.
# 99	Not mapped.

# Plot maps to inspect data
# These are high-resolution maps, so each takes time to render,
par(mfrow=c(4,3)) # split plotting area into 4 rows and 3 columns
plot(DEM,main="Elevation (m)")
plot(ASP,main="Aspect (RAD)")
plot(TWI,main="Topographic Wetness Index")
plot(SOL,main="Solar Radiation")
plot(GEO,main="Bedrock (nutrient class)", breaks = c(0,1,2,3), col = terrain.colors(4))
plot(TMP,main="Temperature (C, annual mean)")
plot(TMP2,main="Temperature (C, warmest quarter)")
plot(PRC,main="Precipitation (mm, annual)")
plot(SW5,main="Snow water equivalent (mm, May)")
plot(SW1,main="Snow water equivalent (mm, January)")
plot(GSL,main="Growing season length (days)")
plot(LC,main="Land Cover AR50")
par(mfrow=c(1,1)) # set plotting area back to normal

# Inspect distributions of values
hist(DEM, maxpixels=1000) # continuous
hist(GEO, maxpixels=1000) # categorical

# Calculate terrain parameters from DEM on your own - instead of loading pre-calculated
# this will take extra time (2 min on my PC)
start_time <- Sys.time() 
TERR <- terrain(DEM,c("slope","aspect","TPI", "TRI"),  unit='degrees')
end_time <- Sys.time()
end_time - start_time
plot(TERR)

# You can build a raster stack of predictors loaded individually as follows:
# predictors <- stack(DEM,TERR,TWI,SOL,as.factor(GEO),as.factor(SUB),as.factor(LC),TMP,RCP85)

# We use an alternative method for loading several files at once,
# which keeps file names as variable names:
r.list <- list.files(path, pattern="tif$", full.names=TRUE) # add all raster files in a given folder to a list
predictors <- stack(r.list)                                 # read and stack raster layers
names(predictors)
# If required, the rasters can be renamed: names(r.stack) <- c("land_cover",...)

# Convert raster layers with factor variables to factors
predictors[["ar50_artype"]] <- as.factor(predictors[["ar50_artype"]])
predictors[["geo_norge123"]] <- as.factor(predictors[["geo_norge123"]])

# Plot individual raster in the stack of predictor data
plot(predictors[[2]])
plot(predictors[["bioclim_1"]])
# or plot all 
plot(predictors)
# subset layers from rasterstack
LC <- subset(predictors, 1)
# find and save current extent of a layer
extent.LC <- extent(LC)
# modify extent to southern Norway
ymax(extent.LC) <- 7100000
xmax(extent.LC) <- 400000

# Limit to Oslo and surrounding area extent 
# This extent is great for fast testing of scripts
extent.OSL <- extent(220000, 280000, 6620000, 6680000)

# Check plotting
plot(predictors[[1]])
plot(extent.LC, add=TRUE, col="red")
plot(extent.OSL, add=TRUE, col="blue")

# Crop out from raster based on extent
predictors.OSLO <- crop(predictors, extent.OSL)

# Check for correlation between variables (testing only on small extent for practical reasons)
p_cor <- layerStats(predictors.OSLO, 'pearson', na.rm=T)
pairs(p_cor)
# or for example checking correlation like this (only 4 test layers)
pairs(predictors.OSLO[[2:5]])
# second and third are strongly correlated
# you can choose to drop one of them
predictors.OSLO <- dropLayer(predictors.OSLO,3)

# Adjust resolution of rasters
res(predictors.OSLO)
predictors.OSLO.2km <- aggregate(predictors.OSLO, fact=20) # from original resolution 100m down to 2km
res(predictors.OSLO.2km)

# Saving your raster back to disc
path <- "lab_data/" 
writeRaster(predictors.OSLO.2km, filename=paste0(path, names(predictors.OSLO.2km), "_2km"), bylayer=TRUE, format="GTiff", overwrite=TRUE )


### Response variable: species occurrences

# Download GBIF occurrence data
library(rgbif)      # global biodiversity data (GBIF) tools
library(maps)       # Provides functions that let us plot the maps
library(mapdata)    # Contains the hi-resolution points that mark out the countries.require(maps)
require(mapproj)
require(maptools) # mapping tools for spatial objects
# require(rgdal) # provides the Geospatial Data Abstraction Library
# require(raster) # spatial raster data management, works well with dismo

# The approach below gives quick access to data sets of limited size, suitable
# for our purpose in this course.
# However, for scientific publications you should use asynchronous downloading,
# as demonstrated here by Anders Finstad in a workshop at the Nordic Oikos conference:
# https://gbif-europe.github.io/nordic_oikos_2018_r/s3_gbif_demo/3.x_async_download_gbif.Rmd
# This allows downloading larger data sets, and citation of a download with a single doi.

key <- name_backbone(name="Picea sitchensis", kingdom="Plantae")$speciesKey
sp <- occ_search(taxonKey=key, hasCoordinate=TRUE, country="NO", limit=1000) 
gbif_citation(sp)                # Overview of data sources with references, including doi. See https://www.gbif.org/tool/81747/rgbif#citations

# Convert lat-long coordinates to coordinate system of environmental raster data
occ_points <- data.frame(x=sp$data$decimalLongitude,y=sp$data$decimalLatitude)
occ_points <- SpatialPoints(occ_points,proj4string=CRS("+proj=longlat +datum=WGS84"))
occ_UTM33 <- spTransform(occ_points,CRS("+proj=utm +zone=33 ellps=GRS80 +units=m"))
sp$data$x <- occ_UTM33$x
sp$data$y <- occ_UTM33$y

# Point data like 'sp' can be used directly as input to MIAmaxent (Labs 4-5),
# along with environmental rasters.
# For glm (Lab 3), we need absence data in addition to presences, and we
# need to combine presences and absences with environmental data in a
# data set suitable for the glm function. This can be done as follows.

# Make background raster for Norway
norway <- predictors[["dem100"]]/predictors[["dem100"]]           # raster with values=1 in mainland Norway, based on elevation raster
plot(norway)

# Rasterize occurrences
occ_ras <- rasterize(occ_UTM33,norway,fun='count',background=0)   # raster with counts of occurrences in each cell of norway
plot(occ_ras)
occ_ras <- occ_ras*norway                                         # filtering occurrence raster with mainland raster
occ_ras[occ_ras>0] <- 1                                           # reducing counts>0 to 1 (presence)
plot(occ_ras)
zoom(occ_ras)                                                     # zoom by clicking twice to western Norway to see the rasterized presence pixels
plot(occ_ras, xlim = c(-74000,10000), ylim = c(6610000, 6680000)) # zoom to western Norway to see the rasterized presence pixels
summary(occ_ras)
table(values(occ_ras))

## Similar plots as above, with different color scales
# plot(occ_ras,col=colorRampPalette(c("blue","red"))(2))
# plot(occ_ras,col=colorRampPalette(c("blue","red"))(2), xlim = c(-74000,10000), ylim = c(6610000, 6680000) )

# Take the occurrence cells as presences
presences <- which(values(occ_ras)==1)

# Then generate absence data by sampling from the cells without occurrence observations.
# NB! This rests on risky assumptions, but we need the absences
# for logistic regression with glm (Lab 3). We make this absence data set for
# educational purposes, but ideally one would want a data set with true, verified presences
# and absences for logistic regression (and also for validation of 'presence-only' methods,
# such as maxent).
absences <- which(values(occ_ras)==0)
absences_sample <- sample(absences,size=length(presences)) # sample of same number of absence cells as presence cells

# Combine presences, absences and environmental data
selected <- c(presences,absences_sample)
xy <- coordinates(occ_ras)[selected, ]
training_data <- data.frame(xy,presence=values(occ_ras)[selected],
                  extract(predictors,xy))
head(training_data)
tail(training_data)
plot(training_data[,c("x","y")],col=c("blue","red")[training_data$presence+1])

# Convert discrete environmental predictors to factor variables
training_data$ar50_artype <- factor(training_data$ar50_artype)
training_data$geo_norge123 <- factor(training_data$geo_norge123)

# Save data (you'll overwrite any existing files with these names if you run these lines)
save(training_data,file="lab_data/Norway_sitka_training_data")
save(predictors,file="lab_data/Norway_predictor_raster_stack")




# Global data - do if you have time
#------------
# The code below is how training data for the DM primer (Lab1) were prepared.
# Go through the steps make sure you understand what the functions are doing. Look up with '?'
# Then change the species and/or predictors to something of your own choice (i.e. manipulate the code below),
# make a new data set and repeat the exercies in Lab1.


# Environmental predictor variable data: download WorldClim data, extract temperature
env <- raster::getData("worldclim",var="bio",res=10)  # Biologically relevant climate variables: http://www.worldclim.org/bioclim
temp <- env[[1]]; names(temp) <- "temp"       # Selecting annual mean temperature and naming it 'temp'

# Response variable data: download GBIF presence data
species <- "Picea sitchensis"
n <- 100
key <- name_backbone(name=species, kingdom="Plantae")$speciesKey
occ <- occ_search(taxonKey=key, return="data", hasCoordinate=TRUE, limit=n) 

# Plot observations on predictor map
plot(temp,main=paste("Temperature and","GBIF occurrences of ",species))
occ_points <- data.frame(x=occ$data$decimalLongitude,y=occ$data$decimalLatitude)
points(occ_points)

# Generate random absence points across MAP EXTENT and add to map
abs_points <- spsample(as(env@extent, 'SpatialPolygons'),n=n, type="random")  # simple random points, including sea
points(abs_points,pch="+",col="red")

# Generate random absence points across PREDICTOR DATA (takes some time) and make new map
abs_points <- spsample(as(temp, 'SpatialPolygons'),n=n, type="random")
plot(temp,main=paste("Temperature and","GBIF occurrences of ",species))
points(occ_points)
points(abs_points,pch="+",col="red")

# Extract data from rasters and polygons at specific points (e.g. presences and absences of a species)
presence_data <- data.frame(occ_points,temp=extract(temp,occ_points))
absence_data <- data.frame(abs_points,temp=extract(temp,abs_points))
training_data <- rbind(presence_data,absence_data)
training_data$presence <- c(rep(1,nrow(presence_data)),rep(0,nrow(absence_data)))
head(training_data)
tail(training_data)

# Save primer training data and predictor map (overwrites any existing files with these names)
# path <- "lab_data/"
# save(training_data,file=paste(path,"DM_primer_training_data",sep=""))
# writeRaster(predictor_map,file=paste(path,"DM_primer_predictor_map",sep=""))



# When you're done, here's a link for inspiration and further independent study:
# GIS with R (shift+click) https://pakillo.github.io/GISwithR/
