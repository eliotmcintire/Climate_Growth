#' to generate the data for analysis for MBdata
#'
#' @param inputData data.table 
#' 
#' @param DBHCutoff numeric, specify the cut off point for trees into analyses
#' 
#' @param dominanceLevel numeric, specify the dominance levels for trees in a stand
#'                               
#' 
#' @param dominanceMethod character, specify the method that classify the dominant layer
#'                        "quantile" and "equalBA"
#' 
#' @param superDominance numeric, specify the first N biggest trees used for super dominance 
#'                       status, super dominant trees are the trees keep within the N biggest status
#'                       and alive throught the monitoring period
#' @param minObsForTree numeric, specify the minimum observations for each trees
#'
#' @return a data table
#' 
#' @importFrom data.table data.table ':='
#' @importFrom dplyr left_join '%>%' 
#'
#' @note no note
#'
#' @seealso no
#'
#' @export
#' @docType methods
#' @rdname MBanalysesDataGenerator
#'
#' @author Yong Luo
#'
setGeneric("MBanalysesDataGenerator",
           function(inputData,
                    DBHCutoff,
                    dominanceLevel,
                    dominanceMethod,
                    superDominance,
                    minObsForTree) {
             standardGeneric("MBanalysesDataGenerator")
           })

#' @export
#' @rdname MBanalysesDataGenerator
setMethod("MBanalysesDataGenerator",
          signature = c(inputData = "data.table",
                        DBHCutoff = "numeric",
                        dominanceLevel = "numeric",
                        dominanceMethod = "character",
                        superDominance = "numeric",
                        minObsForTree = "numeric"),
          definition = function(inputData,
                                DBHCutoff,
                                dominanceLevel,
                                dominanceMethod,
                                superDominance,
                                minObsForTree){
            if(dominanceMethod == "equalBA"){
              allData <- data.frame(inputData) %>% data.table # a copy of whole data
            }   
            inputDataRef <- inputData[IniDBH >= DBHCutoff,]
            inputDataRef <- inputDataRef[,.(referenceYear=min(IniYear)), by = uniTreeID]
            inputData <- setkey(inputData, uniTreeID)[setkey(inputDataRef, uniTreeID), nomatch = 0]
            inputData <- inputData[IniYear>=referenceYear,]
            inputData[,NumberofTree:=length(IniDBH), by = c("PlotID", "IniYear")]
            inputData[,minNumberofTree:=min(NumberofTree), by = PlotID]
            inputData <- inputData[minNumberofTree>=30,][,':='(NumberofTree = NULL,
                                                               minNumberofTree = NULL)]
            superDomiTrees <- inputData[order(PlotID, IniYear,  -IniDBH, -IniHeight),]
            superDomiTrees[,DBHrank:=1:length(IniDBH), by = c("PlotID", "IniYear")]
            superDomiTrees <- superDomiTrees[DBHrank <= superDominance,]
            superDomiTrees[,':='(plotMeasureLength=length(unique(IniYear))), by = PlotID]
            superDomiTrees[,':='(treeMeasureLength=length(unique(IniYear))), by = uniTreeID]
            superDomiTrees <- superDomiTrees[plotMeasureLength == treeMeasureLength,.(uniTreeID)] %>%
              unique(.,by = "uniTreeID")
            analysesTrees <- inputData[IniYear == referenceYear,.(PlotID, uniTreeID, IniYear)]
            #length(unique(analysesTrees$uniTreeID)) == nrow(analysesTrees)
            inputData[, referenceYear:=NULL]
            dataForDominance <- unique(analysesTrees[,.(PlotID, IniYear)],
                                       by = c("PlotID", "IniYear"))
            if(dominanceMethod == "quantileDBH"){
              dataForDominance <- setkey(inputData, PlotID, IniYear)[setkey(dataForDominance,
                                                                            PlotID, IniYear),
                                                                     nomatch = 0]
              dataForDominance[,NumberofTree:=length(IniDBH), by = c("PlotID", "IniYear")]
              
              plotlevels <- unique(dataForDominance$PlotID)
              for(indiplot in plotlevels){
                indiplotdata <- dataForDominance[PlotID == indiplot,]
                indiplotdata <- indiplotdata[order(IniYear, IniDBH, IniHeight),]
                indiplotdata <- indiplotdata[,DBHrank:=1:length(Species), by = IniYear]
                # indiplotdata[NumberofTree <= (dominanceLevel+1), DominanceClass:=paste("D", dominanceLevel, sep = "")]
                indiplotdata[, DominanceClass:=as.character(cut(DBHrank,
                                                                round(c(0, 0.5, 0.9, 1)*length(Species)),
                                                                labels = c(paste("D", 1:dominanceLevel, sep = "")),
                                                                include.lowest = TRUE)),
                             by = IniYear]
                indiplotdata <- indiplotdata[,.(uniTreeID, IniYear, DominanceClass)]
                if(indiplot == plotlevels[1]){
                  dominanceTable <- indiplotdata
                } else {
                  dominanceTable <- rbind(dominanceTable, indiplotdata)
                }
                rm(indiplotdata)
              }
              
            } else if(dominanceMethod == "quantileHeight"){
              dataForDominance <- setkey(inputData, PlotID, IniYear)[setkey(dataForDominance,
                                                                            PlotID, IniYear),
                                                                     nomatch = 0]
              dataForDominance[,NumberofTree:=length(IniDBH), by = c("PlotID", "IniYear")]
              
              plotlevels <- unique(dataForDominance$PlotID)
              for(indiplot in plotlevels){
                indiplotdata <- dataForDominance[PlotID == indiplot,]
                indiplotdata <- indiplotdata[order(IniYear, IniHeight, IniDBH),]
                indiplotdata <- indiplotdata[,heightrank:=1:length(Species), by = IniYear]
                # indiplotdata[NumberofTree <= (dominanceLevel+1), DominanceClass:=paste("D", dominanceLevel, sep = "")]
                indiplotdata[, DominanceClass:=as.character(cut(heightrank,
                                                                round(c(0, 0.3, 0.6, 0.9, 1)*length(Species)),
                                                                labels = c(paste("D", 1:dominanceLevel, sep = "")),
                                                                include.lowest = TRUE)),
                             by = IniYear]
                indiplotdata <- indiplotdata[,.(uniTreeID, IniYear, DominanceClass)]
                if(indiplot == plotlevels[1]){
                  dominanceTable <- indiplotdata
                } else {
                  dominanceTable <- rbind(dominanceTable, indiplotdata)
                }
                rm(indiplotdata)
              }
              
            } else if(dominanceMethod == "equalBA"){
              dataForDominance <- setkey(inputData, PlotID, IniYear)[setkey(dataForDominance,
                                                                            PlotID, IniYear),
                                                                     nomatch = 0]
              # include all the trees for a given plot and initial year
              dataForDominance[,NumberofTree:=length(IniDBH), by = c("PlotID", "IniYear")]
              dataForDominance[,BA:=3.1415926*((IniDBH/2)^2)]
              dataForDominance[,totalBA:=sum(BA), by = c("PlotID", "IniYear")]
              dataForDominance <- dataForDominance[order(PlotID, IniYear, BA),]
              dataForDominance[,BArank:=1:length(IniDBH), by = c("PlotID", "IniYear")]
              dominanceTable <- data.table(uniTreeID = character(), IniYear = numeric(),
                                           DominanceClass = character())
              plotlevels <- unique(inputData$PlotID)
              for(indiplot in plotlevels){
                plotdata <- dataForDominance[PlotID == indiplot,]
                yearlevels <- unique(plotdata$IniYear)
                for(year in yearlevels){
                  plotyeardata <- plotdata[IniYear == year,]
                  cutpoint <- rep(nrow(plotyeardata),dominanceLevel) 
                  BAsum <- plotyeardata$BA[1]
                  for(i in 2:nrow(plotyeardata)){
                    prevBAsum <- BAsum
                    BAsum <- BAsum+plotyeardata$BA[i]
                    for(j in 1:dominanceLevel){
                      if(prevBAsum <= plotyeardata$totalBA[1]*j/dominanceLevel &
                         BAsum > plotyeardata$totalBA[1]*j/dominanceLevel){
                        cutpoint[j] <- i
                      }
                    }
                  }
                  for(k in 1:length(cutpoint)){
                    if(length(unique(cutpoint)) < length(cutpoint)){
                      cutpointTemp <- data.table(cutpoint=cutpoint, temcol = 1:length(cutpoint))
                      cutpointTemp[,Freqcy:=length(temcol), by = cutpoint]
                      cutpointonce <- cutpointTemp[Freqcy==1,]$cutpoint
                      cutpointmulti <- unique(cutpointTemp[Freqcy>1,], by = c("cutpoint", "Freqcy"))
                      for(cutpointind in unique(cutpointmulti$cutpoint)){
                        for(freqcyind in 1:(cutpointmulti[cutpoint==cutpointind,]$Freqcy-1)){
                          cutpointmultiexp <- cutpointmulti[cutpoint==cutpointind,][
                            ,cutpoint:=cutpoint-freqcyind]
                          cutpointmulti <- rbind(cutpointmulti, cutpointmultiexp)
                        }
                      }
                      cutpoint <- sort(c(cutpointonce, cutpointmulti$cutpoint))
                      rm(cutpointTemp,cutpointonce, cutpointmulti, cutpointind, freqcyind)
                    }
                  }
                  plotyeardata[,DominanceClass := as.character(cut(BArank, breaks = c(1, cutpoint),
                                                                   labels = c(paste("D", 1:dominanceLevel, sep = "")),
                                                                   include.lowest = TRUE))]
                  dominanceTable <- rbind(dominanceTable,plotyeardata[,.(uniTreeID, IniYear, DominanceClass)])
                  
                }
              }
            } else if(dominanceMethod == "equalBiomass"){
              dataForDominance <- setkey(inputData, PlotID, IniYear)[setkey(dataForDominance,
                                                                            PlotID, IniYear),
                                                                     nomatch = 0]
              # include all the trees for a given plot and initial year
              dataForDominance[,NumberofTree:=length(IniDBH), by = c("PlotID", "IniYear")]
              dataForDominance[,BA:=3.1415926*((IniDBH/2)^2)]
              dataForDominance[,totalBA:=sum(BA), by = c("PlotID", "IniYear")]
              dataForDominance <- dataForDominance[order(PlotID, IniYear, BA),]
              dataForDominance[,BArank:=1:length(IniDBH), by = c("PlotID", "IniYear")]
              dominanceTable <- data.table(uniTreeID = character(), IniYear = numeric(),
                                           DominanceClass = character())
              plotlevels <- unique(inputData$PlotID)
              for(indiplot in plotlevels){
                plotdata <- dataForDominance[PlotID == indiplot,]
                yearlevels <- unique(plotdata$IniYear)
                for(year in yearlevels){
                  plotyeardata <- plotdata[IniYear == year,]
                  cutpoint <- rep(nrow(plotyeardata),dominanceLevel) 
                  BAsum <- plotyeardata$BA[1]
                  for(i in 2:nrow(plotyeardata)){
                    prevBAsum <- BAsum
                    BAsum <- BAsum+plotyeardata$BA[i]
                    for(j in 1:dominanceLevel){
                      if(prevBAsum <= plotyeardata$totalBA[1]*j/dominanceLevel &
                         BAsum > plotyeardata$totalBA[1]*j/dominanceLevel){
                        cutpoint[j] <- i
                      }
                    }
                  }
                  for(k in 1:length(cutpoint)){
                    if(length(unique(cutpoint)) < length(cutpoint)){
                      cutpointTemp <- data.table(cutpoint=cutpoint, temcol = 1:length(cutpoint))
                      cutpointTemp[,Freqcy:=length(temcol), by = cutpoint]
                      cutpointonce <- cutpointTemp[Freqcy==1,]$cutpoint
                      cutpointmulti <- unique(cutpointTemp[Freqcy>1,], by = c("cutpoint", "Freqcy"))
                      for(cutpointind in unique(cutpointmulti$cutpoint)){
                        for(freqcyind in 1:(cutpointmulti[cutpoint==cutpointind,]$Freqcy-1)){
                          cutpointmultiexp <- cutpointmulti[cutpoint==cutpointind,][
                            ,cutpoint:=cutpoint-freqcyind]
                          cutpointmulti <- rbind(cutpointmulti, cutpointmultiexp)
                        }
                      }
                      cutpoint <- sort(c(cutpointonce, cutpointmulti$cutpoint))
                      rm(cutpointTemp,cutpointonce, cutpointmulti, cutpointind, freqcyind)
                    }
                  }
                  plotyeardata[,DominanceClass := as.character(cut(BArank, breaks = c(1, cutpoint),
                                                                   labels = c(paste("D", 1:dominanceLevel, sep = "")),
                                                                   include.lowest = TRUE))]
                  dominanceTable <- rbind(dominanceTable,plotyeardata[,.(uniTreeID, IniYear, DominanceClass)])
                  
                }
              }
            }
            
            treeWithDominance <- setkey(analysesTrees,
                                        uniTreeID, IniYear)[setkey(dominanceTable,
                                                                   uniTreeID, IniYear),
                                                            nomatch = 0][,.(uniTreeID, DominanceClass)]
            
            analysesData <- setkey(inputData, uniTreeID)[setkey(treeWithDominance, uniTreeID),
                                                         nomatch = 0]
            
            analysesDataForSupperDomi <- inputData[uniTreeID %in% superDomiTrees$uniTreeID, ][
              , DominanceClass:=paste("D", dominanceLevel+1, sep = "")]
            analysesDataForSupperDomi[,uniTreeID:=paste(uniTreeID, "_S", sep = "")]
            rm(superDomiTrees)
            analysesData <- rbind(analysesData, analysesDataForSupperDomi) 
            analysesData <- analysesData[FinYear-IniYear == 5,]
            analysesData[,treeMeasureLength:=length(IniYear), by = uniTreeID]
            analysesData <- analysesData[treeMeasureLength >= (minObsForTree-1),][,treeMeasureLength:=NULL]
            analysesData[,':='(Year = (IniYear+FinYear)/2, FA = (IniFA+FinFA)/2,
                               BAGR = (3.1415926*((FinDBH/2)^2-(IniDBH/2)^2))/5)]
            return(analysesData)
          })

