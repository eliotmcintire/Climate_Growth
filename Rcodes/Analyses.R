rm(list=ls())
library(data.table);library(ggplot2); library(dplyr); library(nlme)
library(SpaDES); library(MuMIn)
load("~/LandWeb/Data/MB/MBPSP.RData")
names(allPSP) <- c("PlotID", "FMU", "TWP", "RGE", "PlotNumber",
                   "StandType", "StructureType", "VegeType",
                   "SStructure", "Moist", "YearEstabish", "Easting",
                   "Northing", "TreeNumber", "Species", "Distance",
                   "Angle", "DBH", "Height", "Status", "Class",
                   "TreeAge", "Year", "PlotSize")
orgData <- data.frame(allPSP) %>% data.table
print(length(unique(allPSP$PlotID))) # 425
# plot level selection
# 1. Natureal regeneration
allPSP <- allPSP[StandType == "NAT",]
print(length(unique(allPSP$PlotID))) # 271
# 2. minimum 3 times
allPSP[,measureTime:=length(unique(Year)), by = PlotID]
allPSP <- allPSP[measureTime>=3,]
print(length(unique(allPSP$PlotID))) # 172
# 3. plot size bigger than 500m^2
allPSP[,plotsizelength := length(unique(PlotSize)), by = PlotID]
allPSP[plotsizelength == 2,.(maxDistance=max(Distance)), by = PlotID] 
# PlotID maxDistance
# 1: 61-295       12.55
#*unique(allPSP[plotsizelength == 2,]$PlotSize)
# [1] 500.34 500.31
allPSP[plotsizelength == 2, PlotSize:=500.34]
#*unique(allPSP$PlotSize)
# [1] 500.34  50.01 250.17
allPSP <- allPSP[PlotSize>=500,]
print(length(unique(allPSP$PlotID))) # 169
set(allPSP, , c("measureTime", "plotsizelength"), NULL)
# 4. FA information
allPSP[, ':='(baseTreeAge= TreeAge-Year+min(Year), baseYear=min(Year)), by = PlotID]
allPSP[,baseFA:=round(mean(baseTreeAge, na.rm = TRUE)), by = PlotID]
allPSP[,FA:=Year-baseYear+baseFA]
#* range(allPSP$FA)
# [1]   5 168
allPSP[,':='(baseTreeAge = NULL, baseYear = NULL, baseFA = NULL)]

# CIindex <- read.csv("F:/MBgrowth/MBdata.csv", header = TRUE,
#                     stringsAsFactors = FALSE) %>%
#   data.table
# CIindex[,uniTreeID:=paste(PlotID, "_", TreeNumber, sep = "")]
# CIindex <- CIindex[,.(uniTreeID, Year, Hegyi, IntraHegyiRatio)]
# write.csv(CIindex, file.path("F:/MBgrowth/competitionIndex.csv"),
#           row.names = FALSE)

####$$$$%%%%%%%%%% tree level selection
# 1. select all alive trees
allPSP[,':='(uniTreeID = paste(PlotID, "_", TreeNumber, sep = ""))]
allPSP <- allPSP[Status == 1 | Status == 0,]

# 2. measured at least three times
allPSP[,':='(maxLength = length(DBH)),
       by = uniTreeID]
allPSP <- allPSP[maxLength >2,]
allPSP <- allPSP[order(uniTreeID, Year),]
allPSP <- allPSP[,c("FinYear", "FinDBH",  "tempuniTreeID",
                     "FinFA") 
                  := shift(x = allPSP[,.(Year, DBH, uniTreeID,
                                         FA)],
                           n = 1, fill = NA, type = "lead")]
allPSP <- allPSP[tempuniTreeID == uniTreeID,][,tempuniTreeID := NULL]
allPSP <- allPSP[Year != FinYear,]

allPSP <- setkey(allPSP, uniTreeID, Year)[setkey(CIindex, uniTreeID, Year), nomatch = 0]
write.csv(allPSP, file.path("F:/MBgrowth/MBdata.csv"),
          row.names = FALSE)




