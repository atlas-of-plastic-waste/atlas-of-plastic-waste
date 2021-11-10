############################################################################
## Mapping Plastic Waste Sites 
## Matt Gordon (m.gordon@yale.edu) and Anna Papp (ap3907@columbia.edu)
## Last modified: 11/10/21, AP
############################################################################

############################################################################
# PART 0 - Setup
############################################################################

# clear
rm(list = ls())

# set working directory for project
setwd("/Users/annapapp/Desktop/atlas-of-plastic-waste/gee")

### libraries ###
#install.packages("rgee")
#install.packages("googledrive")
#install.packages("googleCloudStorageR")
#install.packages("leafem")
#install.packages("raster")
#install.packages("geojsonio")
library(rgee)
library(googledrive)
#library(googleCloudStorageR)
library(leafem)
library(geojsonio)
library(sf)
library(raster)
#ee_install(py_env = "rgee") # only need to do this once 
ee_Initialize(user = 'annapapp94@gmail.com', drive=TRUE) # if need to add GCS, need to change

### geometries
testarea <- ee$Geometry$Polygon(list(c(100.30806104215795, 5.7343542479935),c(100.25312940153295, 5.351631511547638),c(100.54976026090795, 3.741777202509949),c(101.24189893278295, 2.8753122233284727),c(101.89009229215795, 2.436332193522161),c(102.28560010465795, 3.1715298579348636),c(100.67060987028295, 6.138667549099112)))
apwtable <- ee$FeatureCollection("users/annnapappp/atlas_of_plastic")
training19 <- ee$FeatureCollection("projects/atlasplasticwaste/assets/training19")

### constants ###
# date filters 
datefilter19 <- ee$Filter$date('2019-01-01', '2019-12-31')

# smoothing radius for Sentinel-1
smoothingradius <- 50

# bands for classification 
#bands <- c('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12','VV', 'VH')
bands <- c('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12')

# image visualization parameters 
sentinel2VisParams <- list(min = 0.0,max = 3000,bands = c("B4", "B3", "B2"))
classificationVisParams <- list(min = 0.0,max = 1.0, bands=c("classification"))
sentinel1VisParams <- list(min=-15, max=0, bands=c("VH"))

############################################################################
# PART 1 - Load Dump Sites
############################################################################

# function for setting class and label of dumpsites 
setclass <- function(feature){
  feature$set(list(class=1, label=dumpsites))
}

# filter imported shapefile - separate file for each year
apw19 <- ee$FeatureCollection(apwtable)$
  filter(ee$Filter$eq("year", "2019"))

############################################################################
# PART 2 - Sentinel-2 Imagery
############################################################################

# cloud function helper 
getQABits <- function(image, qa) {
  # Convert decimal (character) to decimal (little endian)
  qa <- sum(2^(which(rev(unlist(strsplit(as.character(qa), "")) == 1))-1))
  # Return a single band image of the extracted QA bits, giving the qa value.
  image$bitwiseAnd(qa)$lt(1)
}

# sentinel cloud mask 
maskS2clouds <- function(image) {
  qa <- image$select('QA60')
  mask <- getQABits(qa, "110000000000")
  image$updateMask(mask)
}

# get sentinel-2 level2 image for each year, median for now
sentinel219 <- ee$ImageCollection('COPERNICUS/S2_SR')$
  filterBounds(testarea)$
  filter(datefilter19)$
  filter(ee$Filter$lt('CLOUDY_PIXEL_PERCENTAGE',20))$
  map(maskS2clouds)$
  median()$
  clip(testarea)

############################################################################
# PART 3 - SAR Imagery
############################################################################

# get Sentinel-1 data for entire time period 
s1 <- ee$ImageCollection('COPERNICUS/S1_GRD')$
  filterDate('2018-01-01', '2020-12-31')$
  filterMetadata('resolution_meters', 'equals' , 10)$
  filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VV'))$
  filter(ee$Filter$listContains('transmitterReceiverPolarisation', 'VH'))$
  filter(ee$Filter$eq('instrumentMode', 'IW'))$
  filterBounds(testarea)

# separate ascending and descending orbit images into distinct collections
s1ascend <- s1$filter(ee$Filter$eq('orbitProperties_pass', 'ASCENDING'))
s1descend <- s1$filter(ee$Filter$eq('orbitProperties_pass', 'DESCENDING'))

# mean vh ascending / descending, apply filter to reduce speckle
#vvAsc19 <- s1ascend$filter(datefilter19)$select('VV')$mean()$focal_mean(smoothingradius, 'circle', 'meters')
#vvDesc19 <- s1ascend$filter(datefilter19)$select('VV')$mean()$focal_mean(smoothingradius, 'circle', 'meters')
#vhAsc19 <- s1ascend$filter(datefilter19)$select('VH')$mean()$focal_mean(smoothingradius, 'circle', 'meters')
#vhDesc19 <- s1descend$filter(datefilter19)$select('VH')$mean()$focal_mean(smoothingradius, 'circle', 'meters')
vvAsc19 <- s1ascend$filter(datefilter19)$select('VV')$mean()
vvDesc19 <- s1ascend$filter(datefilter19)$select('VV')$mean()
vhAsc19 <- s1ascend$filter(datefilter19)$select('VH')$mean()
vhDesc19 <- s1descend$filter(datefilter19)$select('VH')$mean()

# combine into 1 image for each year 
sentinel119 <- ee$Image$cat(vvAsc19,vhAsc19)

############################################################################
# PART 4 - Classification
############################################################################

# combine images from Sentinel 2 and Sentinel 1
sentinel19 <- sentinel219$addBands(sentinel119)$select(bands)

# sample training points from the sentinel imagery 
sampleSet19 <- sentinel19$sampleRegions(collection = training19,
                                        properties=list("class"),
                                        scale = 30,
                                        tileScale = 16)

# create a random forest classifier with 200 trees and get model training on the sampleSet
RF_LCZ19 <- ee$Classifier$smileRandomForest(200)$train(sampleSet19, 'class', bands)

#classify the map and visualize it
LCZmap19 <- sentinel19$classify(RF_LCZ19)
classified <- LCZmap19
#classified <- LCZmap19$updateMask(LCZmap19$eq(1))

############################################################################
# PART 5 - Map
############################################################################

# add sites to map (do not yet include classification)
Map$centerObject(testarea, 8)
Map$addLayer(eeObject = apw19,
             visParams = list(color= 'FF0000'),
             name = 'Plastic Dumps - 2019') +
Map$addLayer(eeObject=sentinel219,
             visParams=sentinel2VisParams,
             name='Sentinel-2 Composite 2019') +
Map$addLayer(eeObject=sentinel119, 
             visParams=sentinel1VisParams,
             name='Sentinel-1 Ascending') + 
  Map$addLayer(eeObject=classified, 
               visParams=classificationVisParams,
               name='Classification') 

# to do, export?
#task <- ee$batch$Export$image$toDrive(
#  image = wastesites,
#  description = 'classification',
#  scale = 30,
#  region = testarea,
#  maxPixels = 100000000
#)
#task$start()
#ee_monitoring(task)

#exported_stats <- ee_drive_to_local(task = task,dsn = "exported_stats.csv")
#read.csv(exported_stats)


