rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "F:/MBgrowth"
} else {
  workPath <- "J:/MBgrowth"
}
analysesComponent <- "BAGR" # or "BiomassGR

analysesData <- read.csv(file.path(workPath, "MBdataSimplified.csv"), header = TRUE,
                         stringsAsFactors = FALSE) %>%
  data.table
analysesData[,':='(IntraHegyi = Hegyi*IntraHegyiRatio,
                   InterHegyi = Hegyi*(1-IntraHegyiRatio))]
if (analysesComponent == "BAGR"){
  analysesData[,':='(logY = log(BAGR), 
                     logDBH = log(IniDBH)-mean(log(IniDBH)), 
                     Yearctd = GSCMIA-mean(GSCMIA),
                     logIntraHegyi = log(IntraHegyi+0.1)-mean(log(IntraHegyi+0.1)),
                     logInterHegyi = log(InterHegyi+0.1)-mean(log(InterHegyi+0.1)))]
} else if(analysesComponent == "BiomassGR"){
  analysesData[,':='(logY = log(BiomassGR), 
                     logDBH = log(IniDBH)-mean(log(IniDBH)), 
                     logIntraHegyi = log(IntraHegyi+0.1)-mean(log(IntraHegyi+0.1)),
                     logInterHegyi = log(InterHegyi+0.1)-mean(log(InterHegyi+0.1)))]
}

OverallModels <- list()
OverallResults <- data.table(Species = character(), Dominance = numeric(),
                             R2 = numeric(),
                             Model = character(), NofObservation = numeric(),
                             NofTree = numeric(), Variable = character(),
                             Value = numeric(), Std.Error = numeric(),
                             DF = numeric(), t.value = numeric(),
                             p.value = numeric())
optim <- lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000)
testSpecieses <- c("JP", "BS", "TA", "BF")
i <- 1
for(testspecies in testSpecieses){
  theModel <- lme(logY~logDBH+Yearctd+logIntraHegyi+logInterHegyi+
                    logDBH:Yearctd+logDBH:logIntraHegyi+logDBH:logInterHegyi+
                    Yearctd:logIntraHegyi+Yearctd:logInterHegyi+
                    logIntraHegyi:logInterHegyi,
                  random = ~1|PlotID,
                  data = analysesData[Species == testspecies],
                  control = optim)
  OverallModels[[i]] <- theModel
  names(OverallModels)[i] <- paste(testspecies, "_all", sep = "")
  i <- i+1
  coeff <- data.frame(summary(theModel)$tTable)
  coeff$Variable <- row.names(coeff)
  coeff <- data.table(coeff)[,':='(Species = testspecies, Dominance = 0,
                                   R2 = as.numeric(r.squaredGLMM(theModel)[2]),
                                   Model = paste("Overall_", testspecies, sep = ""),
                                   NofObservation = nrow(analysesData[Species == testspecies,]),
                                   NofTree = length(unique(analysesData[Species == testspecies,]$uniTreeID)))]
  OverallResults <- rbind(OverallResults, coeff[,.(Species, Dominance, R2, Model,
                                                   NofObservation, NofTree, Variable,
                                                   Value, Std.Error, DF, t.value, p.value)])
}
rm(i, coeff, dominanceAssigner, optim, testspecies, testSpecieses, theModel)
print(AIC(OverallModels[[1]])) # 39252.37 for ATA
# 39242.47 for GSTA

# 38237.03 for ACMIA
# 38187.82 for GSCMIA



# minObsForTree = 3 # the trees was observed at least 2 times 

rmdead <- F
if(rmdead){
  analysesData[,maxdead:=max(Dead), by = uniTreeID]
  analysesData <- analysesData[maxdead == 0, ][,maxdead:=NULL]
}
rmrecruitment <- F
if(rmrecruitment){
  analysesData[,maxrec:=max(Regen), by = uniTreeID]
  analysesData <- analysesData[maxrec == 0, ][,maxrec:=NULL]
}
rm(rmdead, rmrecruitment)
# analysesData[,treeMeasureLength:=length(IniDBH), by = uniTreeID]
# analysesData <- analysesData[treeMeasureLength>=2,][,treeMeasureLength:=NULL]


# remove both regen and dead trees
# analysesData <- analysesData[Regen==0 & Dead==0,]

optim <- lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000)

allresults <- data.table(Species = character(), Dominance = numeric(),
                         R2 = numeric(),
                         Model = character(), NofObservation = numeric(),
                         NofTree = numeric(), Variable = character(),
                         Value = numeric(), Std.Error = numeric(),
                         DF = numeric(), t.value = numeric(),
                         p.value = numeric())
DomianceAllModels <- list()
YearResiduals <- data.table(Species = character(), Dominance = numeric(),
                            uniTreeID = character(),
                            Year = numeric(), OrigBAGR = numeric(),
                            FittedBAGR = numeric(), Residuals = numeric())
targetSpecies <- c("JP", "BS", "TA", "BF")
# DominanceLevels <- c("Suppressed", "ModerateSuppressed", "Codominant", "Dominant")
DominanceLevels <- c(paste("D", 1:4, sep = ""))
relativeImportanceOutput <- data.table(Species = character(), Dominance = numeric(),
                                       Variable = character(), Fvalue = numeric(),
                                       Pvalue = numeric())
i <- 1
for(indiSpecies in targetSpecies){
  analysesData_Species <- analysesData[Species == indiSpecies,]
  for(dominance in DominanceLevels){
    analysesData_Species_Dominance <- analysesData_Species[DominanceClass == dominance,]
    theModelDBHYear <- lme(logY~logDBH+Yearctd+logIntraHegyi+logInterHegyi+
                             logDBH:Yearctd+logDBH:logIntraHegyi+logDBH:logInterHegyi+
                             Yearctd:logIntraHegyi+Yearctd:logInterHegyi+
                             logIntraHegyi:logInterHegyi,
                           random = ~1|PlotID,
                           data = analysesData_Species_Dominance,
                           control = optim)
    indiResult <- data.frame(summary(theModelDBHYear)$tTable)
    indiResult$Variable <- as.character(row.names(indiResult))
    indiResult <- data.table(indiResult)[,':='(R2 = as.numeric(r.squaredGLMM(theModelDBHYear)[2]))]
    
    indiResult[,':='(Model = paste(indiSpecies, "_", dominance, sep = ""),
                     Species = indiSpecies,
                     Dominance = as.numeric(gsub("D", "", dominance)),
                     NofObservation = nrow(analysesData_Species_Dominance),
                     NofTree = length(unique(analysesData_Species_Dominance$uniTreeID)))]
    
    allresults <- rbind(allresults, indiResult[,.(Species, Dominance, R2,
                                                  Model, NofObservation, NofTree,
                                                  Variable, Value, Std.Error, DF,
                                                  t.value, p.value)])
    relativeImportance <- data.frame(anova.lme(theModelDBHYear))
    relativeImportance$Variable <- row.names(relativeImportance)
    relativeImportance <- data.table(relativeImportance)[,':='(Species = indiSpecies,
                                                               Dominance = as.numeric(gsub("D", "", dominance)))]
    relativeImportanceOutput <- rbind(relativeImportanceOutput,
                                      relativeImportance[,.(Species = indiSpecies,
                                                            Dominance = as.numeric(gsub("D", "", dominance)),
                                                            Variable,
                                                            Fvalue = F.value,
                                                            Pvalue = p.value)])
    DomianceAllModels[[i]] <- theModelDBHYear
    names(DomianceAllModels)[i] <- paste(indiSpecies, "_", dominance, sep = "")
    i <- i+1
  }
}
allresults <- rbind(allresults, OverallResults)
rm(analysesData_Species, analysesData_Species_Dominance, dominance, DominanceLevels,
   i, indiResult, indiSpecies, optim, targetSpecies, YearResiduals)
save.image(file.path(workPath,paste(analysesComponent, "GSCMIA_mainResults.RData")))
