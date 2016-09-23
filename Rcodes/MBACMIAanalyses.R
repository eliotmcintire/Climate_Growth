rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "F:/MBgrowth"
} else {
  workPath <- "J:/MBgrowth"
}


inputData <- read.csv(file.path(workPath, "MBdataSimplified.csv"), header = TRUE,
                      stringsAsFactors = FALSE) %>%
  data.table



dominanceLevel <- 3
DBHCutoff <- 5 # 5 cm is the cut off point for Initial DBH
superDominance <- 5 # first biggest 5 trees as super dominant status
minObsForTree <- 3
dominanceMethod <- "quantileDBH" # or "quantile"

source(file.path(workPath,"MBanalysesDataGenerator.R"))
# dominanceLevel <- 4
analysesData <- MBanalysesDataGenerator(inputData = inputData,
                                        DBHCutoff = DBHCutoff, # 5 cm is the cut off point for Initial DBH
                                        dominanceLevel = dominanceLevel, # 4 dominant levels
                                        dominanceMethod = dominanceMethod,
                                        superDominance = superDominance, # first biggest 5 trees as super dominant status
                                        minObsForTree = 3) # the trees was observed at least 2 times 

analysesData <- analysesData[BAGR>0,]
analysesData[,treeMeasureLength:=length(IniDBH), by = uniTreeID]
analysesData <- analysesData[treeMeasureLength>=2,][,treeMeasureLength:=NULL]

analysesData[,':='(minPlotYear=min(IniYear), maxPlotYear=max(FinYear)), by = PlotID]
analysesData[,':='(minTreeYear=min(IniYear), maxTreeYear=max(FinYear)), by = uniTreeID]
analysesData[minPlotYear==minTreeYear, Regen:=0]
analysesData[minPlotYear!=minTreeYear, Regen:=1]
analysesData[maxPlotYear==maxTreeYear, Dead:=0]
analysesData[maxPlotYear!=maxTreeYear, Dead:=1]

speciesTable <- unique(analysesData[,.(Species, uniTreeID, Regen, Dead)], by = "uniTreeID")[, ':='(NumberofTree=length(uniTreeID),
                                                                                                   NumberofRegen = sum(Regen),
                                                                                                   NumberofDead = sum(Dead)), by = Species] %>%
  unique(., by = "Species")
speciesTable <- speciesTable[order(-NumberofTree),]
set(analysesData, ,c("minPlotYear", "maxPlotYear", "minTreeYear", "maxTreeYear"),NULL)

# remove both regen and dead trees
# analysesData <- analysesData[Regen==0 & Dead==0,]

analysesData[,':='(logBAGR = log(BAGR), 
                   logDBH = log(IniDBH)-mean(log(IniDBH)), 
                   logACMIA = log(ACMIA+10)-mean(log(ACMIA+10)),
                   logHegyiIntra = log(Hegyi*IntraHegyiRatio+0.001)-
                     mean(log(Hegyi*IntraHegyiRatio+0.001)),
                   logHegyiInter = log(Hegyi*(1-IntraHegyiRatio)+0.001)-
                     mean(log(Hegyi*(1-IntraHegyiRatio)+0.001)))]
optim <- lmeControl(opt="optim", maxIter=10000, msMaxIter = 10000)


allresults <- data.table(Species = character(), Dominance = numeric(),
                         ModelSeq = character(), R2 = numeric(),
                         Model = character(), NofObservation = numeric(),
                         NofTree = numeric(), Variable = character(),
                         Value = numeric(), Std.Error = numeric(),
                         DF = numeric(), t.value = numeric(),
                         p.value = numeric())
allModels <- list()
ACMIAResiduals <- data.table(Species = character(), Dominance = numeric(),
                            uniTreeID = character(),
                            ACMIA = numeric(), OrigBAGR = numeric(),
                            FittedBAGR = numeric(), Residuals = numeric())
targetSpecies <- speciesTable$Species[1:3]
# DominanceLevels <- c("Suppressed", "ModerateSuppressed", "Codominant", "Dominant")
DominanceLevels <- c(paste("D", 1:(dominanceLevel+1), sep = ""))
relativeImportanceOutput <- data.table(Species = character(), Dominance = numeric(),
                                       Variable = character(), Fvalue = numeric(),
                                       Pvalue = numeric())
i <- 1
for(indiSpecies in targetSpecies){
  analysesData_Species <- analysesData[Species == indiSpecies,]
  for(dominance in DominanceLevels){
    analysesData_Species_Dominance <- analysesData_Species[DominanceClass == dominance,]
    theModelDBHACMIA <- lme(logBAGR~logDBH+logHegyiIntra+logHegyiInter+logACMIA+
                             logDBH:logHegyiIntra+logDBH:logHegyiInter+logDBH:logACMIA+
                             logHegyiIntra:logHegyiInter+logHegyiIntra:logACMIA+
                             logHegyiInter:logACMIA,
                           random = ~1|PlotID/uniTreeID,
                           data = analysesData_Species_Dominance,
                           control = optim)
    indiResult <- data.frame(summary(theModelDBHACMIA)$tTable)
    indiResult$Variable <- as.character(row.names(indiResult))
    indiResult <- data.table(indiResult)[,':='(ModelSeq = "Simultaneous",
                                               R2 = as.numeric(r.squaredGLMM(theModelDBHACMIA)[2]))]
    
    indiResult[,':='(Model = paste(indiSpecies, "_", dominance, sep = ""),
                     Species = indiSpecies,
                     Dominance = as.numeric(gsub("D", "", dominance)),
                     NofObservation = nrow(analysesData_Species_Dominance),
                     NofTree = length(unique(analysesData_Species_Dominance$uniTreeID)))]
    
    allresults <- rbind(allresults, indiResult[,.(Species, Dominance, ModelSeq, R2,
                                                  Model, NofObservation, NofTree,
                                                  Variable, Value, Std.Error, DF,
                                                  t.value, p.value)])
    relativeImportance <- data.frame(anova.lme(theModelDBHACMIA))
    relativeImportance$Variable <- row.names(relativeImportance)
    relativeImportance <- data.table(relativeImportance)[,':='(Species = indiSpecies,
                                                               Dominance = as.numeric(gsub("D", "", dominance)))]
    relativeImportanceOutput <- rbind(relativeImportanceOutput,
                                      relativeImportance[,.(Species = indiSpecies,
                                                            Dominance = as.numeric(gsub("D", "", dominance)),
                                                            Variable,
                                                            Fvalue = F.value,
                                                            Pvalue = p.value)])
    allModels[[i]] <- theModelDBHACMIA
    names(allModels)[i] <- paste(indiSpecies, "_", dominance, sep = "")
    i <- i+1
  }
}
# 
figureImportance <- ggplot(data = relativeImportanceOutput[Variable == "logDBH" | Variable == "logHegyiIntra" |
                                                             Variable == "logHegyiIntra" | Variable == "logHegyiInter" |
                                                             Variable == "logACMIA",], 
                           aes(x = as.factor(Dominance), y = Fvalue))+
  geom_line(aes(col = Variable, group = Variable), size = 1.5)+
  scale_color_manual(values = c("red", "green", "blue", "black"))+
  scale_y_log10(name = "Relative importance",
                breaks = c(0, 1, 10, 100, 1000, 10000))+
  scale_x_discrete(name = "none", labels = c("Severe suppressed", "Suppressed",
                                             "Codominant", "Dominant",
                                             "Super dominant"))+
  # geom_errorbar(aes(x = Dominance, ymin = Value-1.98*Std.Error, ymax = Value+1.98*Std.Error, col = ModelSeq))+
  # geom_hline(yintercept = 0, col = "blue")+
  facet_wrap(~Species)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        legend.position = c(0.5, 0.2))



# & ModelSeq != "Simultaneous"
figureACMIAEffectChange <- ggplot(data = allresults[Variable == "logACMIA", ], 
                                 aes(x = as.factor(Dominance), y = Value))+
  geom_point(aes(group = Species))+
  scale_y_continuous(name = "ACMIA effect")+
  scale_x_discrete(name = "none", labels = c( "Suppressed",
                                             "Codominant", "Dominant",
                                             "Super dominant"))+
  geom_errorbar(aes(x = Dominance, ymin = Value-1.98*Std.Error, ymax = Value+1.98*Std.Error, col = ModelSeq))+
  geom_hline(yintercept = 0, col = "blue")+
  facet_wrap(~Species)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        legend.position = "none")

figureIntraCIEffectChange <- ggplot(data = allresults[Variable == "logHegyiIntra", ], 
                                    aes(x = as.factor(Dominance), y = Value))+
  geom_point(aes(group = Species))+
  scale_y_continuous(name = "Intra-CI effect")+
  scale_x_discrete(name = "none", labels = c("Severe suppressed", "Suppressed",
                                             "Codominant", "Dominant",
                                             "Super dominant"))+
  geom_errorbar(aes(x = Dominance, ymin = Value-1.98*Std.Error, ymax = Value+1.98*Std.Error, col = ModelSeq))+
  geom_hline(yintercept = 0, col = "blue")+
  facet_wrap(~Species)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        legend.position = "none")

figureInterCIEffectChange <- ggplot(data = allresults[Variable == "logHegyiInter", ], 
                                    aes(x = as.factor(Dominance), y = Value))+
  geom_point(aes(group = Species))+
  scale_y_continuous(name = "Inter-CI effect")+
  scale_x_discrete(name = "none", labels = c("Severe suppressed", "Suppressed",
                                             "Codominant", "Dominant",
                                             "Super dominant"))+
  geom_errorbar(aes(x = Dominance, ymin = Value-1.98*Std.Error, ymax = Value+1.98*Std.Error, col = ModelSeq))+
  geom_hline(yintercept = 0, col = "blue")+
  facet_wrap(~Species)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        legend.position = "none")

figureIntraCIACMIAEffectChange <- ggplot(data = allresults[Variable == "logHegyiIntra:logACMIA", ], 
                                        aes(x = as.factor(Dominance), y = Value))+
  geom_point(aes(group = Species))+
  scale_y_continuous(name = "Intra-CI*ACMIA effect")+
  scale_x_discrete(name = "none", labels = c("Suppressed",
                                             "Codominant", "Dominant",
                                             "Super dominant"))+
  geom_errorbar(aes(x = Dominance, ymin = Value-1.98*Std.Error, ymax = Value+1.98*Std.Error, col = ModelSeq))+
  geom_hline(yintercept = 0, col = "blue")+
  facet_wrap(~Species)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        legend.position = "none")

figureInterCIACMIAEffectChange <- ggplot(data = allresults[Variable == "logHegyiInter:logACMIA", ], 
                                        aes(x = as.factor(Dominance), y = Value))+
  geom_point(aes(group = Species))+
  scale_y_continuous(name = "Inter-CI*ACMIA effect")+
  scale_x_discrete(name = "none", labels = c("Severe suppressed", "Suppressed",
                                             "Codominant", "Dominant",
                                             "Super dominant"))+
  geom_errorbar(aes(x = Dominance, ymin = Value-1.98*Std.Error, ymax = Value+1.98*Std.Error, col = ModelSeq))+
  geom_hline(yintercept = 0, col = "blue")+
  facet_wrap(~Species)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        legend.position = "none")

