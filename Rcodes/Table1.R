rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn); library(gridExtra)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Associates/Yong Luo/Climate_Growth/MBgrowth"
} else {
  workPath <- "J:/MBgrowth"
}
DBHCut <- 0
analysesData <- read.csv(file.path(workPath, "MBdataSimplified.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>%
  data.table

analysesData[,':='(TreeMLength = length(IniDBH)),
             by = uniTreeID]

analysesData[, ':='(minPlotMY = min(IniYear),
                    PlotMLength = length(unique(IniYear))), by = PlotID]
minFirstDBH <- analysesData[IniYear==minPlotMY, .(uniTreeID, FirstDBH = IniDBH)]
analysesData <- dplyr::left_join(analysesData, minFirstDBH, by = "uniTreeID") %>% data.table
analysesData[,Status:="IngrowthOrDead"]


analysesData[TreeMLength == PlotMLength & FirstDBH >= DBHCut, Status:="Survive"]
set(analysesData, , c("TreeMLength", "PlotMLength"), NULL)
minDBHdata <- analysesData[Status == "Survive",.(PlotID, IniYear, IniDBH, IniHeight)][
  ,.(minDBH=min(IniDBH), minHeight = min(IniHeight)), by = c("PlotID", "IniYear")]
analysesData <- setkey(analysesData, PlotID, IniYear)[setkey(minDBHdata, PlotID, IniYear), 
                                                      nomatch = 0]
analysesData[,':='(maxDBH=max(IniDBH), maxHeight = max(IniHeight)), by = c("PlotID", "IniYear")]
analysesData[,DominanceIndex:=100*(IniDBH-minDBH)/(maxDBH-minDBH)]
analysesData[, minBAGR:=round(min(BAGR), 3), by = uniTreeID]
analysesData <- analysesData[Status == "Survive" & minBAGR > 0,]


climateinputs <- read.csv(file.path(workPath, "plotClimates.csv"), header = TRUE,
                          stringsAsFactors = FALSE) %>%
  data.table
nrow(analysesData) #72486
analysesData <- setkey(analysesData, PlotID, IniYear, FinYear)[setkey(climateinputs, PlotID, IniYear, FinYear),
                                                               nomatch = 0]

nrow(analysesData) #72486

specieses <- c("JP", "BS", "TA")
allspeciesdata <- analysesData[Species %in% specieses, ]
for(indispecies in specieses){
  speciesdata <- allspeciesdata[Species == indispecies,.(PlotID, BAGR, IniDBH, Year, Hegyi, DominanceIndex,
                                                         ATA, GSTA, NONGSTA,
                                                         APA, GSPA, NONGSPA,
                                                         ACMIA, GSCMIA, NONGSCMIA,
                                                         ACO2A, GSCO2A, NONGSCO2A)]
  tempcolnames <- names(speciesdata)[2:18]
  speciesdata <- reshape(data = speciesdata, varying = tempcolnames,
                         v.names = "Values",
                         # times = tempcolnames, 
                         timevar = "Variable",
                         direction = "long")
  speciesdata <- speciesdata[,.(mean = round(mean(Values), 2), sd = round(sd(Values), 2),
                                min = round(min(Values), 2),
                                max = round(max(Values), 2)), by = Variable]
  speciesdata <- speciesdata[,.(Variable, mean = paste(mean, "±", sd), 
                                range = paste(min, "to", max))]
  names(speciesdata)[2:3] <- paste(indispecies, "_", names(speciesdata)[2:3], sep = "")
  if(indispecies == specieses[1]){
    Table1Output <- speciesdata
  } else {
    Table1Output <- setkey(Table1Output, Variable)[setkey(speciesdata, Variable), nomatch = 0]
  }
}
Table1Output$Variable <- tempcolnames
write.csv(Table1Output, file.path(workPath, "table1.csv"), row.names = FALSE)

for(indispecies in specieses){
  speciesdata <- allspeciesdata[Species == indispecies,]
  cat("Number of plot:", length(unique(speciesdata$PlotID)), "\n")
  cat("Number of tree:", length(unique(speciesdata$uniTreeID)), "\n")
  cat("Number of observation:", nrow(speciesdata), "\n")
}
