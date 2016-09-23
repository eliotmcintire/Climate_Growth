rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "F:/MBgrowth"
} else {
  workPath <- "J:/MBgrowth"
}
analysesComponent <- "BAGR"
load(file.path(workPath, paste(analysesComponent, "_mainResults.RData")))
studySpecies <- c("JP", "TA", "BS")
output <- data.table(Species = character(), Dominance = numeric(),
                     Variable = character(), Value = numeric())
for(indispecies in studySpecies){
  themodel <- OverallModels[[paste(indispecies,"_all", sep = "")]]
  speciesData <- data.frame(analysesData[Species == indispecies,])
  # for DBH
  predictData <- speciesData %>% data.table
  predictData[,':='(Yearctd = 0, logIntraHctd = 0, logInterHctd = 0,
                    DBHbyYearctd = 0, DBHbyIntraHctd = 0, DBHbyInterctd = 0,
                    YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
  DBHsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
  # for Year
  predictData <- speciesData %>% data.table
  predictData[,':='(logDBHctd = 0, logIntraHctd = 0, logInterHctd = 0,
                    DBHbyYearctd = 0, DBHbyIntraHctd = 0, DBHbyInterctd = 0,
                    YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
  
  Yearsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
  # for intraCI
  predictData <- speciesData %>% data.table
  predictData[,':='(Yearctd = 0, logDBHctd = 0, logInterHctd = 0,
                    DBHbyYearctd = 0, DBHbyIntraHctd = 0, DBHbyInterctd = 0,
                    YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
  
  IntraCIsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
  # for interCI
  predictData <- speciesData %>% data.table
  predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                    DBHbyYearctd = 0, DBHbyIntraHctd = 0, DBHbyInterctd = 0,
                    YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
  InterCIsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
  # for DBHbyYearctd
  predictData <- speciesData %>% data.table
  predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                    logInterHctd = 0, DBHbyIntraHctd = 0, DBHbyInterctd = 0,
                    YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
  DBHbyYearsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
  # for DBHbyIntraHctd
  predictData <- speciesData %>% data.table
  predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                    logInterHctd = 0, DBHbyYearctd = 0, DBHbyInterctd = 0,
                    YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
  DBHbyIntraHsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
  # for DBHbyInterctd
  predictData <- speciesData %>% data.table
  predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                    logInterHctd = 0, DBHbyYearctd = 0, DBHbyIntraHctd = 0,
                    YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
  DBHbyInterHsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
  # for YearbyIntraHctd
  predictData <- speciesData %>% data.table
  predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                    logInterHctd = 0, DBHbyYearctd = 0, DBHbyIntraHctd = 0,
                    DBHbyInterctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
  YearbyIntraHsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
  # for YearbyInterHctd
  predictData <- speciesData %>% data.table
  predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                    logInterHctd = 0, DBHbyYearctd = 0, DBHbyIntraHctd = 0,
                    DBHbyInterctd = 0, YearbyIntraHctd = 0, InterHbyIntraHctd = 0)]
  YearbyInterHsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
  
  
  output <- rbind(output, data.table(Species = indispecies,
                                     Dominance = 0,
                                     Variable = c("DBH", "Year", "IntraCI", "InterCI",
                                                  "DBHYear", "DBHIntraCI", "DBHInterCI",
                                                  "YearIntraCI", "YearInterCI"),
                                     Value = c(DBHsenstivity, Yearsenstivity,
                                               IntraCIsenstivity, InterCIsenstivity,
                                               DBHbyYearsenstivity, DBHbyIntraHsenstivity,
                                               DBHbyInterHsenstivity, YearbyIntraHsenstivity,
                                               YearbyInterHsenstivity)))
  rm(DBHsenstivity, Yearsenstivity, IntraCIsenstivity, InterCIsenstivity, themodel)
  for(i in 1:4){
    speciesDominanceData <- data.frame(data.table(speciesData)[DominanceClass == paste("D", i, sep = ""),])
    themodel <- DomianceAllModels[[paste(indispecies,"_D", i, sep = "")]]
    predictData <- speciesDominanceData %>% data.table
    predictData[,':='(Yearctd = 0, logIntraHctd = 0, logInterHctd = 0,
                      DBHbyYearctd = 0, DBHbyIntraHctd = 0, DBHbyInterctd = 0,
                      YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
    DBHsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
    # for Year
    predictData <- speciesDominanceData %>% data.table
    predictData[,':='(logDBHctd = 0, logIntraHctd = 0, logInterHctd = 0,
                      DBHbyYearctd = 0, DBHbyIntraHctd = 0, DBHbyInterctd = 0,
                      YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
    
    Yearsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
    # for intraCI
    predictData <- speciesDominanceData %>% data.table
    predictData[,':='(Yearctd = 0, logDBHctd = 0, logInterHctd = 0,
                      DBHbyYearctd = 0, DBHbyIntraHctd = 0, DBHbyInterctd = 0,
                      YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
    
    IntraCIsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
    # for interCI
    predictData <- speciesDominanceData %>% data.table
    predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                      DBHbyYearctd = 0, DBHbyIntraHctd = 0, DBHbyInterctd = 0,
                      YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
    InterCIsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
    # for DBHbyYearctd
    predictData <- speciesDominanceData %>% data.table
    predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                      logInterHctd = 0, DBHbyIntraHctd = 0, DBHbyInterctd = 0,
                      YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
    DBHbyYearsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
    # for DBHbyIntraHctd
    predictData <- speciesDominanceData %>% data.table
    predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                      logInterHctd = 0, DBHbyYearctd = 0, DBHbyInterctd = 0,
                      YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
    DBHbyIntraHsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
    # for DBHbyInterctd
    predictData <- speciesDominanceData %>% data.table
    predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                      logInterHctd = 0, DBHbyYearctd = 0, DBHbyIntraHctd = 0,
                      YearbyIntraHctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
    DBHbyInterHsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
    # for YearbyIntraHctd
    predictData <- speciesDominanceData %>% data.table
    predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                      logInterHctd = 0, DBHbyYearctd = 0, DBHbyIntraHctd = 0,
                      DBHbyInterctd = 0, YearbyInterHctd = 0, InterHbyIntraHctd = 0)]
    YearbyIntraHsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
    # for YearbyInterHctd
    predictData <- speciesDominanceData %>% data.table
    predictData[,':='(Yearctd = 0, logIntraHctd = 0, logDBHctd = 0,
                      logInterHctd = 0, DBHbyYearctd = 0, DBHbyIntraHctd = 0,
                      DBHbyInterctd = 0, YearbyIntraHctd = 0, InterHbyIntraHctd = 0)]
    YearbyInterHsenstivity <- sd(predict(themodel, newdata = predictData, level = 0))
    
    
    output <- rbind(output, data.table(Species = indispecies,
                                       Dominance = i,
                                       Variable = c("DBH", "Year", "IntraCI", "InterCI",
                                                    "DBHYear", "DBHIntraCI", "DBHInterCI",
                                                    "YearIntraCI", "YearInterCI"),
                                       Value = c(DBHsenstivity, Yearsenstivity,
                                                 IntraCIsenstivity, InterCIsenstivity,
                                                 DBHbyYearsenstivity, DBHbyIntraHsenstivity,
                                                 DBHbyInterHsenstivity, YearbyIntraHsenstivity,
                                                 YearbyInterHsenstivity)))
    rm(DBHsenstivity, Yearsenstivity, IntraCIsenstivity, InterCIsenstivity, themodel)
  }
}

figureImportance <- ggplot(data = output, 
                           aes(x = as.factor(Dominance), y = Value))+
  geom_line(aes(col = Variable, group = Variable), size = 1.5)+
  # scale_color_manual(values = c("red", "green", "blue", "black"))+
  scale_y_continuous(name = "Relative importance",
                breaks = seq(0, 1, length = 5))+
  scale_x_discrete(name = "none", labels = c("Overall", "Severe suppressed", 
                                             "Suppressed", "Codominant",
                                             "Dominant"))+
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

signalTable <- output[Variable == "Year",.(Signal = sum(Value)), by = c("Species", "Dominance")]
signalTableAdd <- output[Variable == "InterCI" | Variable == "IntraCI",.(Noise = sum(Value)), by = c("Species", "Dominance")]
signalTable <- setkey(signalTable, Species, Dominance)[setkey(signalTableAdd, Species, Dominance),
                                                       nomatch = 0]
signalTable[,Ratio:=Signal/Noise]
rm(signalTableAdd)

figureSignalNoise <- ggplot(data = signalTable, 
                            aes(x = as.factor(Dominance), y = Ratio))+
  geom_point(aes(col = Species), size = 1.5)+
  scale_color_manual(values = c("red", "green", "blue"))+
  scale_y_continuous(name = "Signal-to-noise ratio")+
  scale_x_discrete(name = "none", labels = c("Overall", "Severe suppressed", 
                                             "Suppressed", "Codominant",
                                             "Dominant"))+
  # geom_errorbar(aes(x = Dominance, ymin = Value-1.98*Std.Error, ymax = Value+1.98*Std.Error, col = ModelSeq))+
  # geom_hline(yintercept = 0, col = "blue")+
  # facet_wrap(~Species)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        legend.position = c(0.95, 0.2))





