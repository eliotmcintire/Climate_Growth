rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Associates/Yong Luo/Climate_Growth/MBgrowth"
} else {
  workPath <- "J:/MBgrowth"
}
analysesComponent <- "BAGR" # or "BiomassGR
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
analysesDataOrg <- data.table(data.frame(analysesData))
analysesData1 <- analysesData[Status == "Survive",]
analysesData <- analysesDataOrg[Status == "Survive" & minBAGR > 0,]


speciesTable <- unique(analysesData[,.(Species, uniTreeID)], by = "uniTreeID")[
  , .(NumberofTree=length(uniTreeID)), by = c("Species")] 

speciesTable[,totNumbofTree:=sum(NumberofTree), by = Species]
speciesTable <- speciesTable[order(-totNumbofTree),.(Species,  NumberofTree)]
set(analysesData, ,c("minPlotMY", "FirstDBH", "Status", "minDBH",
                     "minHeight", "maxDBH", "maxHeight", "minBAGR"),NULL)
climateinputs <- read.csv(file.path(workPath, "plotClimates.csv"), header = TRUE,
                          stringsAsFactors = FALSE) %>%
  data.table
nrow(analysesData) #72486
analysesData <- setkey(analysesData, PlotID, IniYear, FinYear)[setkey(climateinputs, PlotID, IniYear, FinYear),
                                                               nomatch = 0]

nrow(analysesData) #72486



OverallModels <- list()
OverallResults <- data.table(Species = character(),
                             Climate = character(),
                             R2 = numeric(),
                             Variable = character(),
                             Value = numeric(), Std.Error = numeric(),
                             DF = numeric(), t.value = numeric(),
                             p.value = numeric())
optim <- lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000)

testSpecieses <- c("JP", "BS", "TA")

climates <- c("ATA", "GSTA", "NONGSTA",
              "APA", "GSPA", "NONGSPA",
              "ACMIA", "GSCMIA", "NONGSCMIA",
              "ACO2A", "GSCO2A", "NONGSCO2A")

i <- 1
for(testspecies in testSpecieses){
  speciesData <- analysesData[Species == testspecies,]
  for(indiclimate in climates){
    if(i > 1){ set(speciesData, ,"climate", NULL)}
    setnames(speciesData, indiclimate, "climate")
    speciesData[,':='(logY = log(BAGR), 
                      logDBHctd = log(IniDBH)-mean(log(IniDBH)), 
                      climatectd = climate-mean(climate),
                      logHctd = log(Hegyi)-mean(log(Hegyi)),
                      Dominancectd = DominanceIndex - mean(DominanceIndex))]
    theModel <- lme(logY~logDBHctd+climatectd+logHctd+Dominancectd+
                      logDBHctd:climatectd+logDBHctd:logHctd+logDBHctd:Dominancectd+
                      climatectd:logHctd+climatectd:Dominancectd+logHctd:Dominancectd+
                      logDBHctd:climatectd:logHctd+logDBHctd:climatectd:Dominancectd+
                      climatectd:logHctd:Dominancectd,
                    random = ~1|PlotID/uniTreeID,
                    data = speciesData,
                    control = optim)
    OverallModels[[i]] <- theModel
    names(OverallModels)[i] <- paste(testspecies, "_", indiclimate, sep = "")
    i <- i+1
    coeff <- data.frame(summary(theModel)$tTable)
    coeff$Variable <- row.names(coeff)
    coeff <- data.table(coeff)[,':='(Species = testspecies, 
                                     Climate = indiclimate,
                                     R2 = as.numeric(r.squaredGLMM(theModel)[2]))]
    OverallResults <- rbind(OverallResults, coeff[,.(Species, Climate, R2,
                                                     Variable,
                                                     Value, Std.Error, DF, t.value, p.value)])
  }
  
}
rm(i, coeff, optim, testspecies, testSpecieses, theModel, analysesComponent,
   analysesData1, analysesDataOrg, DBHCut, minDBHdata, minFirstDBH, speciesData)

save.image(file.path(workPath, "BAGR_ClimateResults.RData"))
