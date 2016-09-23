rm(list = ls())
library(data.table); library(ggplot2); library(SpaDES)
library(nlme); library(dplyr);library(MuMIn)
if(as.character(Sys.info()[6]) == "yonluo"){
  workPath <- "~/Associates/Yong Luo/Climate_Growth/MBgrowth"
} else {
  workPath <- "J:/MBgrowth"
}
load(file.path(workPath, paste("BAGR_ClimateResults.RData", sep = "")))
rm(climateinputs, climates, indiclimate, speciesTable)
tempdata <- data.table(Dominance = seq(0, 100, by = 10))
maineffectTable <- data.table(Species = OverallResults[Variable == "climatectd",]$Species,
                              Climate = OverallResults[Variable == "climatectd",]$Climate,
                              maineffect = OverallResults[Variable == "climatectd",]$Value,
                              maineffect_SE = OverallResults[Variable == "climatectd",]$Std.Error)
interactionTable <- data.table(Species = OverallResults[Variable == "climatectd:Dominancectd",]$Species,
                               Climate = OverallResults[Variable == "climatectd:Dominancectd",]$Climate,
                               interactioneffect = OverallResults[Variable == "climatectd:Dominancectd",]$Value,
                               Pvalue = OverallResults[Variable == "climatectd:Dominancectd",]$p.value)
coeffTable <- setkey(maineffectTable, Species, Climate)[setkey(interactionTable, Species, Climate),
                                                        nomatch = 0]
climateWithDominanceTable <- setkey(coeffTable[,k:=1], k)[setkey(tempdata[,k:=1], k), nomatch = NA, allow.cartesian = TRUE][,k:=NULL]
speciesdata <- analysesData[Species %in% c("JP", "BS", "TA"),]
speciesdata <- speciesdata[,.(meandomi=mean(DominanceIndex)), by = Species]
climateWithDominanceTable <- setkey(climateWithDominanceTable, Species)[setkey(speciesdata, Species),
                                                                        nomatch = 0]
climateWithDominanceTable[,climateEffect:=(Dominance-meandomi)*interactioneffect+maineffect]
climateWithDominanceTable[, linetype:=1]
climateWithDominanceTable[Pvalue>=0.05, linetype:=2]
climateWithDominanceTable$Species <- factor(climateWithDominanceTable$Species, levels = c("JP", "TA", "BS"),
                                            labels = c("Jack pine", "Trembling aspen", "Black spruce"))
climateWithDominanceTable$linetype <- factor(climateWithDominanceTable$linetype, levels = c(1, 2),
                                             labels = c(1, 2))
climateWithDominanceTable[,newDomincneIndex:=Dominance]
climateWithDominanceTable[Climate %in% c("GSTA", "GSPA", "GSCMIA", "GSCO2A"), 
                          newDomincneIndex:=newDomincneIndex+110]
climateWithDominanceTable[Climate %in% c("NONGSTA", "NONGSPA", "NONGSCMIA", "NONGSCO2A"), 
                          newDomincneIndex:=newDomincneIndex+220]

rm(coeffTable, interactionTable, maineffectTable, speciesdata, tempdata)
maineffectTable <- unique(climateWithDominanceTable[,.(Species, Climate, maineffect, maineffect_SE, meandomi)],
                          by = c("Species", "Climate"))
setnames(maineffectTable, c("maineffect", "meandomi"), c("climateEffect", "Dominance"))
maineffectTable[, newDomincneIndex:=Dominance]
maineffectTable[Climate %in% c("GSTA", "GSPA", "GSCMIA", "GSCO2A"), 
                          newDomincneIndex:=newDomincneIndex+110]
maineffectTable[Climate %in% c("NONGSTA", "NONGSPA", "NONGSCMIA", "NONGSCO2A"), 
                          newDomincneIndex:=newDomincneIndex+220]


climates <- read.csv(file.path(workPath, "plotClimates.csv"), header = TRUE,
                          stringsAsFactors = FALSE) %>%
  data.table
climates[,Year:=(FinYear+IniYear)/2]
longcol <- names(climates)[4:15]
climate_longform <- reshape(data = climates, varying = longcol, v.names = "Value",
                            times = longcol, timevar = "DependentVariable", 
                            direction = "long")
climate_longform[,Yearctd:=Year-mean(Year)]

temperatureNames <- c("ATA", "GSTA", "NONGSTA")
precipitationNames <- c("APA", "GSPA", "NONGSPA")
CMINames <- c("ACMIA", "GSCMIA", "NONGSCMIA")
CO2Names <- c("ACO2A", "GSCO2A", "NONGSCO2A")
for(indiclimategroup in c("temperatureNames", "precipitationNames", "CMINames", "CO2Names")){
  climateData <- climate_longform[DependentVariable %in% get(indiclimategroup),]
  climateModel <- lme(Value ~ DependentVariable/Yearctd, random =~(DependentVariable-1)|PlotID,
                      data = climateData)
  climateData$predValue <- predict(climateModel, newdata = climateData, level = 0) 
  coeff <- data.frame(summary(climateModel)$tTable)
  coeff <- coeff[row.names(coeff) %in% paste("DependentVariable", get(indiclimategroup), ":Yearctd", sep = ""),c(1,5)] %>%
    data.table
  coeff[,':='(DependentVariable = get(indiclimategroup),
              Value = round(coeff$Value, 2), linetype=1)]
  coeff[p.value>=0.05, linetype:=2]
  climateData <- setkey(climateData, DependentVariable)[setkey(coeff[,.(DependentVariable, linetype)], DependentVariable),
                                                        nomatch = 0]
  if(indiclimategroup == "temperatureNames"){ 
    allClimateData <- climateData
    allcoeff <- coeff
  } else {
      allClimateData <- rbind(allClimateData, climateData)
      allcoeff <- rbind(allcoeff, coeff)
    }
}
allClimateData$linetype <- factor(allClimateData$linetype, levels = c(1, 2))

print(allcoeff[DependentVariable %in% temperatureNames,.(DependentVariable, Value)])
textYdif <- 0.3
textYcentral <- -1.1
Figure_a <- ggplot(data = allClimateData[DependentVariable %in% temperatureNames],
                            aes(x = Year, y = Value))+
  geom_point(aes(group = interaction(PlotID, DependentVariable), col = DependentVariable),
            alpha = 0.1)+
  geom_line(aes(x = Year, y = predValue, col = DependentVariable, linetype = linetype), size = 1)+
  scale_linetype(guide = "none")+
  scale_color_manual(name = "Season", 
                     values = c("red", "green", "blue"),
                     label = c("Whole year", "Growing season", "Non-growing season"))+
  scale_y_continuous(name = expression(paste("Temperature anomaly (", degree, "C)")), limits = c(-1.5, 1.5),
                     breaks = seq(-1.5, 1.5, by = 0.5))+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  annotate("text", x = c(1985, 1999, rep(2002, 3)), 
           y = c(1.5, textYcentral, textYcentral+textYdif, textYcentral, textYcentral-textYdif),
           label = c("a", "Slope:", as.character(allcoeff[DependentVariable %in% temperatureNames,]$Value)),
           size = c(10, rep(5, 4)), col = c("black", "black", "red", "green", "blue"))+
  annotate("text", x = 2007, y = textYcentral, label = paste("~degree~C~year^{-1}"),
           size = 5, parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = c(0.3, 0.82),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))



# the a
Figure_b <- ggplot(data = climateWithDominanceTable[Climate %in% temperatureNames,], 
                aes(x = Dominance, y = climateEffect))+
  geom_segment(x = 0, xend = 100, y = 0, yend = 0, 
               linetype = 2, size = 1, colour = "gray")+
  geom_line(aes(group = interaction(Species, Climate), col = Species, linetype = linetype), size = 1)+
  geom_point(data = maineffectTable[Climate %in% temperatureNames,],
             aes(x = Dominance, y = climateEffect, col = Species), size = 2)+
  geom_errorbar(data = maineffectTable[Climate %in% temperatureNames,],
                aes(x = Dominance, ymin = climateEffect-1.98*maineffect_SE, 
                    ymax = climateEffect+1.98*maineffect_SE, col = Species))+
  facet_wrap(~Climate)+
  geom_text(data = data.frame(Dominance = 0, climateEffect = 3.2, label = "b", Climate = "ATA"),
            aes(x = Dominance, y = climateEffect, label = label), size = 10)+
  geom_text(data = data.frame(Dominance = 50, climateEffect = 3.0, 
                              label = c("Whole year", "Growing season", "Non-growing season"),
                              Climate = c("ATA", "GSTA", "NONGSTA")),
            aes(x = Dominance, y = climateEffect, label = label), size = 5)+
  scale_color_manual(name = "Species", 
                     values = c("red", "green", "blue"),
                     labels = c("Jack pine", "Trembling aspen", "Black spruce"))+
  scale_linetype(guide = "none")+
  scale_y_continuous(name = "Temperature effect",
                     limits = c(-0.4, 3.2), 
                     breaks = round(seq(-0.4, 3.2, by = 0.6),1))+
  scale_x_continuous(name = "Dominance index", limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = c(0.8, 0.6),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))



PrepFigure <- "FALSE"
if(PrepFigure){
  print(allcoeff[DependentVariable %in% precipitationNames,.(DependentVariable, Value)])
  precipitationFigure_c <- ggplot(data = allClimateData[DependentVariable %in% precipitationNames],
                                  aes(x = Year, y = Value))+
    geom_point(aes(group = interaction(PlotID, DependentVariable), col = DependentVariable),
               alpha = 0.1)+
    geom_line(aes(x = Year, y = predValue, col = DependentVariable, linetype = linetype), size = 1)+
    scale_linetype(guide = "none")+
    scale_color_manual(name = "Season", 
                       values = c("red", "green", "blue"),
                       label = c(expression(paste("Whole year (slope: 3.65 mm ", year^{-1},")")),
                                 expression(paste("Growing season (slope: 2.86 mm ", year^{-1},")")),
                                 expression(paste("Non-growing season (slope: 1.04 mm ", year^{-1},")"))))+
    scale_y_continuous(name = paste("Precipitation anomaly (mm)"), limits = c(-100, 130),
                       breaks = seq(-100, 130, by = 35))+
    scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
    annotate("text", x = 1985, y = 130, label = "c", size = 10)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.title.x = element_blank(),
          legend.position = c(0.8, 0.2),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15))
  
  
  precipitationFigure_d <- ggplot(data = climateWithDominanceTable[Climate %in% precipitationNames,], 
                                  aes(x = newDomincneIndex, y = climateEffect))+
    geom_segment(data = data.frame(x = c(0, 110, 220), xend = c(100, 210, 320)),
                 aes(x = x, xend = xend, y = 0, yend = 0), 
                 linetype = 2, size = 1, colour = "gray")+
    geom_line(aes(group = interaction(Species, Climate), col = Species, linetype = linetype), size = 1)+
    geom_point(data = maineffectTable[Climate %in% precipitationNames,],
               aes(x = newDomincneIndex, y = climateEffect, col = Species), size = 2)+
    geom_errorbar(data = maineffectTable[Climate %in% precipitationNames,],
                  aes(x = newDomincneIndex, ymin = climateEffect-1.98*maineffect_SE, 
                      ymax = climateEffect+1.98*maineffect_SE, col = Species))+
    scale_color_manual(name = "Species", 
                       values = c("red", "green", "blue"),
                       labels = c("Jack pine", "Trembling aspen", "Black spruce"))+
    scale_linetype(guide = "none")+
    scale_y_continuous(name = "Effect of precipitation on growth",
                       limits = c(-0.025, 0.015), 
                       breaks = round(seq(-0.02, 0.01, by = 0.1),2))+
    scale_x_continuous(name = "Dominance index", limits = c(0, 320),
                       breaks = c(seq(0, 100, by = 20), seq(110, 210, by = 20),
                                  seq(220, 320, 20)))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          legend.position = c(0.8, 0.2),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 15))
}


print(allcoeff[DependentVariable %in% CMINames,.(DependentVariable, Value)])
textYdif <- 30
textYcentral <- -100
Figure_c <- ggplot(data = allClimateData[DependentVariable %in% CMINames,],
                              aes(x = Year, y = Value))+
  geom_point(aes(group = interaction(PlotID, DependentVariable), col = DependentVariable),
             alpha = 0.1)+
  geom_line(aes(x = Year, y = predValue, col = DependentVariable, linetype = linetype), size = 1)+
  scale_linetype(guide = "none")+
  scale_color_manual(name = "Season", 
                     values = c("red", "green", "blue"),
                     label = c(expression(paste("Whole year (slope: 4.18 mm ", year^{-1},")")),
                               expression(paste("Growing season (slope: 4.18 mm ", year^{-1},")")),
                               expression(paste("Non-growing season (slope: 0.31 mm ", year^{-1},")"))))+
  scale_y_continuous(name = paste("CMI anomaly (mm)"),
                     limits = c(-140, 130),
                     breaks = seq(-140, 130, by = 40))+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  annotate("text", x = c(1985, 1999, rep(2002, 3)), 
           y = c(130, textYcentral, textYcentral+textYdif, textYcentral, textYcentral-textYdif),
           label = c("c", "Slope:", as.character(allcoeff[DependentVariable %in% CMINames,]$Value)),
           size = c(10, rep(5, 4)), col = c("black", "black", "red", "green", "blue"))+
  annotate("text", x = 2007, y = textYcentral, label = paste("~mm~year^{-1}"),
           size = 5, parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))


CMINames1 <- CMINames[1:2]
climateWithDominanceTable[Climate %in% CMINames1, climateEffect:=climateEffect*100]
maineffectTable[Climate %in% CMINames1,':='(climateEffect = climateEffect*100,
                                            maineffect_SE = maineffect_SE*100)]
Figure_d <- ggplot(data = climateWithDominanceTable[Climate %in% CMINames1,], 
                                aes(x = Dominance, y = climateEffect))+
  geom_segment(data = data.frame(Climate = CMINames1, x = 0, xend = 100, y = 0, yend = 0),
               aes(x = x, xend = xend, y = y, yend = yend), 
               linetype = 2, size = 1, colour = "gray")+
  geom_line(aes(group = interaction(Species, Climate), col = Species, linetype = linetype), size = 1)+
  facet_wrap(~Climate)+
  geom_point(data = maineffectTable[Climate %in% CMINames1,],
             aes(x = Dominance, y = climateEffect, col = Species), size = 2)+
  geom_errorbar(data = maineffectTable[Climate %in% CMINames1,],
                aes(x = Dominance, ymin = climateEffect-1.98*maineffect_SE, 
                    ymax = climateEffect+1.98*maineffect_SE, col = Species))+
  scale_color_manual(name = "Species", 
                     values = c("red", "green", "blue"),
                     labels = c("Jack pine", "Trembling aspen", "Black spruce"))+
  scale_linetype(guide = "none")+
  scale_y_continuous(name = expression(paste("CMI effect (",10^{-2}, ")")),
                     limits = c(-0.9, 1), 
                     breaks = round(seq(-0.9, 1, by = 0.3),1))+
  scale_x_continuous(name = "Dominance index", limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  geom_text(data = data.frame(Dominance = 0, climateEffect = 1, label = "d", Climate = c("ACMIA", "NONCMIA")),
            aes(x = Dominance, y = climateEffect, label = label), size = 10, col = c("black", "white"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_blank(),
        legend.position = "none")


print(allcoeff[DependentVariable %in% CO2Names,.(DependentVariable, Value)])
textYdif <- 4
textYcentral <- -12
Figure_e <- ggplot(data = allClimateData[DependentVariable %in% CO2Names,],
                    aes(x = Year, y = Value))+
  geom_point(aes(group = interaction(PlotID, DependentVariable), col = DependentVariable),
             alpha = 0.1)+
  geom_line(aes(x = Year, y = predValue, col = DependentVariable, linetype = linetype), size = 1)+
  scale_linetype(guide = "none")+
  scale_color_manual(name = "Season", 
                     values = c("red", "green", "blue"),
                     label = c(expression(paste("Whole year (slope: 1.80 ppm ", year^{-1},")")),
                               expression(paste("Growing season (slope: 1.80 ppm ", year^{-1},")")),
                               expression(paste("Non-growing season (slope: 1.79 ppm ", year^{-1},")"))))+
  scale_y_continuous(name = expression(paste(CO[2], " anomaly (ppm)")), limits = c(-20, 22),
                     breaks = seq(-20, 20, by = 8))+
  scale_x_continuous(name = "Year", limits = c(1985, 2010), breaks = seq(1985, 2010, by = 5))+
  annotate("text", x = c(1985, 1999, rep(2002, 3)), 
           y = c(22, textYcentral, textYcentral+textYdif, textYcentral, textYcentral-textYdif),
           label = c("e", "Slope:", "1.80", "1.80", "1.79"),
           size = c(10, rep(5, 4)), col = c("black", "black", "red", "green", "blue"))+
  annotate("text", x = 2007, y = textYcentral, label = paste("~ppm~year^{-1}"),
           size = 5, parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = "none",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 15))
climateWithDominanceTable[Climate %in% CO2Names,climateEffect:=climateEffect*100]
maineffectTable[Climate %in% CO2Names,':='(climateEffect = climateEffect*100, maineffect_SE = maineffect_SE*100)]

Figure_f <- ggplot(data = climateWithDominanceTable[Climate %in% CO2Names,], 
                      aes(x = Dominance, y = climateEffect))+
  geom_segment(x = 0, xend = 100, y = 0, yend = 0, 
               linetype = 2, size = 1, colour = "gray")+
  geom_line(aes(group = interaction(Species, Climate), col = Species, linetype = linetype), size = 1)+
  facet_wrap(~Climate, nrow = 1)+
  geom_point(data = maineffectTable[Climate %in% CO2Names,],
             aes(x = Dominance, y = climateEffect, col = Species), size = 2)+
  geom_errorbar(data = maineffectTable[Climate %in% CO2Names,],
                aes(x = Dominance, ymin = climateEffect-1.98*maineffect_SE, 
                    ymax = climateEffect+1.98*maineffect_SE, col = Species))+
  geom_text(data = data.frame(Dominance = 0, climateEffect = 3.6, label = "f", Climate = "ACO2A"),
            aes(x = Dominance, y = climateEffect, label = label), size = 10)+
  scale_color_manual(name = "Species", 
                     values = c("red", "green", "blue"),
                     labels = c("Jack pine", "Trembling aspen", "Black spruce"))+
  scale_linetype(guide = "none")+
  scale_y_continuous(name = expression(paste(CO[2], " effect (", 10^{-2}, ")")),
                     limits = c(-4.8, 3.6), 
                     breaks = round(seq(-4.8, 3.6, by = 1.2),2))+
  scale_x_continuous(name = "Dominance index", limits = c(0, 100),
                     breaks = seq(0, 100, by = 20))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = 1),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = "none")



Figa_Grob <- ggplotGrob(Figure_a)
Figb_Grob <- ggplotGrob(Figure_b)
Figc_Grob <- ggplotGrob(Figure_c)
Figd_Grob <- ggplotGrob(Figure_d)
Fige_Grob <- ggplotGrob(Figure_e)
Figf_Grob <- ggplotGrob(Figure_f)

# height
Figa_Grob$heights <- Fige_Grob$heights
Figb_Grob$heights <- Figf_Grob$heights
Figc_Grob$heights <- Fige_Grob$heights
Figd_Grob$heights <- Figf_Grob$heights

# width
Figa_Grob$widths <- Figc_Grob$widths
Fige_Grob$widths <- Figc_Grob$widths


Figb_Grob$widths <- Figf_Grob$widths
Figd_Grob$widths <- Figf_Grob$widths


dev(4)
clearPlot()
plotlayout <- rbind(c(1, 2, 2), c(3, 4, 4), c(5, 6, 6))
c <- grid.arrange(Figa_Grob, Figb_Grob, Figc_Grob,Figd_Grob, Fige_Grob, Figf_Grob,
                  layout_matrix = plotlayout)

ggsave(file = file.path(workPath, "ClimateTrandsAndEffects.png"), c,
       width = 18, height = 10)
setwd(orgPath)

