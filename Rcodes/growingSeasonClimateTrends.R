rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn); library(raster)
library(maptools); library(rgeos)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPathMB <- "~/Associates/Yong Luo/Climate_Growth/MBgrowth"
  workPathAB <- "~/Associates/Yong Luo/Climate_Growth/ABgrowth"
} else {
  workPathMB <- "J:/MBgrowth"
}

canadamap <- readRDS(file.path(workPathMB,
                               "StudyAreaClimates_BiomSIM",
                               "canadamap.rds"))
provinces <- c("Quebec", "Nova Scotia", "Saskatchewan", "Alberta", "Newfoundland and Labrador", 
               "British Columbia", "New Brunswick", "Prince Edward Island", 
                "Manitoba", "Ontario")
canadamap <- canadamap[canadamap@data$NAME %in% provinces,]
ecozones <- readRDS(file.path(workPathMB,"ecozones.rds"))
ecozoneses <- c("Boreal Cordillera", "Boreal PLain",  "Boreal Shield",
                "Taiga Shield", "Taiga Plain", "Taiga Cordillera",
                "Hudson Plain")
boreal <- ecozones[ecozones@data$ZONE_NAME %in% ecozoneses,]
boreal <- spTransform(boreal, CRSobj = crs(canadamap))
boreal <- gIntersection(boreal, canadamap)

 df<- data.frame(id = "boreal")
 row.names(df) <- 1
 boreal <- SpatialPolygonsDataFrame(boreal, data = df)
# data collecting
# monthlyData <- read.csv(file.path(workPathMB, "StudyAreaClimates_BiomSIM",
#                                 "background_monthlyclimates.csv"),
#                       header = TRUE, stringsAsFactors = FALSE) %>%
#   data.table
# monthlyData <- monthlyData[,.(ID, Latitude, Longitude,
#                           Year = as.numeric(lapply(strsplit(Year.Month, "/", fixed = TRUE), function(x){x[[1]]})),
#                           Month = as.numeric(lapply(strsplit(Year.Month, "/", fixed = TRUE), function(x){x[[2]]})),
#                           Temp = Mean.TMean...C.,
#                           Prep = Total.Precipitation..mm.,
#                           PET = PET..mm.)]
# monthlyData <- monthlyData[Year>=1984 & Year <= 2011 & Month >=5 & Month <= 10,]
# monthlyData[,CMI:=Prep-PET]
# Climates <- monthlyData[, .(meanTemp = mean(Temp), totalPrep = sum(Prep),
#                             totalPET = sum(PET), totalCMI = sum(CMI)),
#                         by = c("ID", "Latitude", "Longitude", "Year")]
# 
# climateTrends <- data.table(ID = numeric(), Latitude = numeric(),
#                             Longitude = numeric(), TempRate = numeric(),
#                             CMIRate = numeric(), PrepRate = numeric())
# IDs <- unique(Climates$ID)
# for(id in IDs){
#   climateTrendsAdd <- data.table(ID = id, Latitude = Climates[ID == id,]$Latitude[1],
#                                  Longitude = Climates[ID == id,]$Longitude[1],
#                                  TempRate = as.numeric(coef(lm(meanTemp~Year, Climates[ID == id,]))[2]),
#                                  CMIRate = as.numeric(coef(lm(totalCMI~Year, Climates[ID == id,]))[2]),
#                                  PrepRate = as.numeric(coef(lm(totalPrep~Year, Climates[ID == id,]))[2]))
#   
#   climateTrends <- rbind(climateTrends, climateTrendsAdd)
#   }
# xyzdata_a <- data.frame(climateTrends[,.(Longitude, Latitude, CMIRate)])
# studyAreaRaster_CMITrend <- rasterFromXYZ(xyz = xyzdata_a, crs = crs("+proj=longlat"))
# studyAreaRaster_CMITrend <- projectRaster(studyAreaRaster_CMITrend,
#                                  res = c(0.05, 0.05),
#                                  crs = crs(studyAreaRaster_CMITrend))
# studyAreaRaster_CMITrend <- projectRaster(studyAreaRaster_CMITrend, crs = crs(canadamap))
# MBmap <- canadamap[canadamap@data$NAME == "Manitoba",]
# studyAreaRaster_CMITrend <- mask(studyAreaRaster_CMITrend, MBmap)
# studyAreaPoly <- studyAreaRaster_CMITrend
# theValues <- getValues(studyAreaPoly)
# theValues[!is.na(theValues)] <- 1
# studyAreaPoly[] <- theValues
# studyAreaPoly <- rasterToPolygons(studyAreaPoly, dissolve = TRUE)


locations <- read.csv(file.path(workPathMB, "masterTable.csv"),
                      header = TRUE,
                      stringsAsFactors = FALSE) %>%
  data.table
workingData <- read.csv(file.path(workPathMB, "MBdataSimplified.csv"),
                        header = TRUE,
                        stringsAsFactors = FALSE) %>%
  data.table
locations <- locations[PlotID %in% unique(workingData$PlotID),]

locations <- locations[,.(PlotID, Easting, Northing)] %>%
  unique(., by = "PlotID")
MBlocation <- SpatialPoints(locations[,.(Easting, Northing)], proj4string = CRS("+proj=utm +zone=14 datum=NAD83"))
MBlocation <- spTransform(MBlocation, crs(canadamap))
MBlocation <- as.data.frame(MBlocation@coords)
ABlocation <- read.csv(file.path(workPathAB, "plotlocations.csv"), header=TRUE,
                       stringsAsFactors = FALSE) %>%
  data.table %>%
  unique(., by = "Groupnumber")
ABlocation <- ABlocation[,.(Groupnumber, Easting, Northing)]
ABlocation <- SpatialPointsDataFrame(ABlocation[,.(Easting, Northing)], data = ABlocation,
                            proj4string = CRS("+proj=utm +zone=11 datum=NAD83"),
                            match.ID = TRUE)
ABlocation <- spTransform(ABlocation, crs(canadamap))
ABlocation1 <- intersect(ABlocation, boreal)

ABlocation <- as.data.frame(ABlocation@coords)


canadamapall <- fortify(canadamap, region = "NAME") %>%
  data.table
canadamapall[,fill:=1]
borealall <- fortify(boreal, region = "id") %>%
  data.table

# studyareaall <- fortify(studyAreaPoly, region = "CMIRate") 


figure1_a1 <- ggplot(data = canadamapall, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group, fill = as.factor(fill)))+
  scale_fill_manual(values = "gray")+
  geom_polygon(data = borealall, aes(x = long, y = lat, group = group),
               fill = "green")+
  geom_path(aes(group = group), col = "white", size = 1)+
  # geom_path(data = studyareaall, aes(x = long, y = lat, group = group),
  #           col = "red", size = 1.5)+
  theme_bw()+
  geom_point(data = MBlocation,
             aes(x = MBlocation$Easting, y = MBlocation$Northing), col = "blue",
             pch = 1, size = 2)+
  geom_point(data = ABlocation,
             aes(x = ABlocation$Easting, y = ABlocation$Northing), col = "red",
             pch = 1, size = 2)+
  annotate("text", label = "a", x = min(canadamapall$long), 
            y = max(canadamapall$lat), size = 13)+
  annotate("text", label = "BC", x = -1800000, y = 2000000, size = 5)+
  annotate("text", label = "AB", x = -1200000, y = 1700000, size = 5)+
  annotate("text", label = "SK", x = -680000, y = 1500000, size = 5)+
  annotate("text", label = "MB", x = -80000, y = 1500000, size = 5)+
  annotate("text", label = "ON", x = 600000, y = 1200000, size = 5)+
  annotate("text", label = "QC", x = 1500000, y = 1700000, size = 5)+
  annotate("text", label = "USA", x = -680000, y = 900000, size = 5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "in"),
        panel.border = element_rect(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        panel.margin = unit(0, "in"))


studyAreaRaster_CMITrend_p <- rasterToPoints(studyAreaRaster_CMITrend)
studyAreaRaster_CMITrend_p <- data.frame(studyAreaRaster_CMITrend_p)
colnames(studyAreaRaster_CMITrend_p) <- c("long", "lat", "CMITrends")
# 
# ABSKmapdata <- fortify(ABSKmap, region = "NAME")
studyAreaRaster_CMITrend_p$CMITrends_New <- cut(x = studyAreaRaster_CMITrend_p$CMITrends, 
                                                round(seq(-3, 6, length = 10)),
                                                include.lowes = TRUE)



figure1_a <- ggplot(data = studyAreaRaster_CMITrend_p, aes(x = long, y = lat))+
  geom_raster(data = studyAreaRaster_CMITrend_p, aes(fill = as.factor(CMITrends_New)))+
   scale_fill_manual(name = "Temporal trend of \ngrowing season \nclimate moisture index \n(1984 - 2011) (mm per year)",
                     values = rev(terrain.colors(10, alpha = 1)))+

  theme_bw()+
  annotate("text", label = "b", x = min(studyAreaRaster_CMITrend_p$long),
           y = max(studyAreaRaster_CMITrend_p$lat), size = 13)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        
        # plot.margin = unit(c(0, 0, 0, 0), "in"),
        plot.background = element_rect(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        # legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.9, 0.25),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 13),
        panel.margin = unit(0.1, "in"))


xyzdata_b <- data.frame(climateTrends[,.(Longitude, Latitude, TempRate)])
studyAreaRaster_TempTrend <- rasterFromXYZ(xyz = xyzdata_b, crs = crs("+proj=longlat"))
studyAreaRaster_TempTrend <- projectRaster(studyAreaRaster_TempTrend,
                                          res = c(0.05, 0.05),
                                          crs = crs(studyAreaRaster_TempTrend))
studyAreaRaster_TempTrend <- projectRaster(studyAreaRaster_TempTrend, crs = crs(canadamap))
studyAreaRaster_TempTrend <- mask(studyAreaRaster_TempTrend, MBmap)

studyAreaRaster_TempTrend_p <- rasterToPoints(studyAreaRaster_TempTrend)
studyAreaRaster_TempTrend_p <- data.frame(studyAreaRaster_TempTrend_p)
colnames(studyAreaRaster_TempTrend_p) <- c("long", "lat", "TempTrends")
# 
# ABSKmapdata <- fortify(ABSKmap, region = "NAME")
studyAreaRaster_TempTrend_p$TempTrends_New <- cut(x = studyAreaRaster_TempTrend_p$TempTrends, 
                                                c(seq(-0.03, 0.06, by = 0.01), 0.08),
                                                include.lowes = TRUE)


figure1_b <- ggplot(data = studyAreaRaster_TempTrend_p, aes(x = long, y = lat))+
  geom_raster(data = studyAreaRaster_TempTrend_p, aes(fill = as.factor(TempTrends_New)))+
  scale_fill_manual(name = "Temporal trend of \ngrowing season \nclimate moisture index \n(1984 - 2011) (mm per year)",
                    values = rev(terrain.colors(10, alpha = 1)))+
  theme_bw()+
  annotate("text", label = "c", x = min(studyAreaRaster_TempTrend_p$long),
           y = max(studyAreaRaster_TempTrend_p$lat), size = 13)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        
        # plot.margin = unit(c(0, 0, 0, 0), "in"),
        plot.background = element_rect(colour = "black"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        # legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.position = c(0.85, 0.25),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 13),
        panel.margin = unit(0.1, "in"))
# figure1a <-  figure1_a+
#   annotation_custom(grob = figure1_a1, xmin = -600000, xmax = 100000,
#                     ymin = 1880000, ymax = 2500000)


library(gridExtra)
plotlayout <- rbind(c(1, 1), c(1, 1), c(2, 3), c(2, 3), c(2, 3))
dev(4)
a <- grid.arrange(figure1_a1, figure1_a, figure1_b, layout_matrix = plotlayout)

orgPath <- getwd()
setwd(file.path(workPathMB))
ggsave(file = "Figure 1.png", a, width = 9, height = 10)
setwd(orgPath)








