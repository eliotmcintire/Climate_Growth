rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn);library(raster)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "F:/MBgrowth"
} else {
  workPath <- "J:/MBgrowth"
}


inputData <- read.csv(file.path(workPath, "MBdataSimplified.csv"), header = TRUE,
                      stringsAsFactors = FALSE) %>%
  data.table
masterTable <- read.csv(file.path(workPath, "masterTable.csv"),
                        header = TRUE, stringsAsFactors = FALSE) %>%
  data.table

masterTable <- masterTable[PlotID %in% unique(inputData$PlotID),]
source("~/GitHub/landwebNRV/landwebNRV/R/UTMtoLongLat.R")

locationData <- unique(masterTable[,.(PlotID, Easting, Northing, Zone = 14)], by = "PlotID")
locationData <- UTMtoLongLat(UTMTable = locationData)
locationData <- locationData$Transformed
demraster <- raster(file.path(workPath, "MBDEM.tif"))
locationGeoPoint <- SpatialPointsDataFrame(coords = data.frame(locationData[,.(Longitude, Latitude)]),
                                           data = data.frame(locationData[,.(PlotID)]),
                                           proj4string = CRS("+proj=longlat"),
                                           match.ID = TRUE)

locationGeoPoint <- spTransform(locationGeoPoint, CRSobj = crs(demraster))
locationGeoPoint@data$elevations <- extract(demraster, locationGeoPoint)
locationData <- setkey(locationData, PlotID)[setkey(data.table(locationGeoPoint@data), PlotID),
                                                nomatch = 0]

locationData <- locationData[!is.na(elevations),.(Name = PlotID, ID = PlotID,
                                                  Latitude, Longitude,
                                                  Elevation = elevations)]

write.csv(locationData,
          file.path(workPath, "StudyAreaClimates_BiomSIM", "plotlocations.csv"),
          row.names = FALSE)




