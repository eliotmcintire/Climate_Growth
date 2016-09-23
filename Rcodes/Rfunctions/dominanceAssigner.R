#' to assign dominance information for the dataset for each census and each plot, therefore, the dominance
#' level can be change for a given tree. The dominance index is calculated for ranked trees from the smallest to biggest tree (BA based or Biomass based).
#' 
#' 
#' 
#' 
#' 
#'
#'
#' @param inputData data.table 
#' 
#' 
#' @param dominanceMethod character, specify the method that classify the dominant layer,
#'                        "BA" and "Biomass"
#' 
#' @param dominanceCut numeric, from 0 to 1, specify the cut point for dividing the trees into dominant level 
#'                     based on the percentage of total BA or total Biomass
#'                                                   
#' 
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
#' @rdname dominanceAssigner
#'
#' @author Yong Luo
#'
setGeneric("dominanceAssigner",
           function(inputData,
                    dominanceCut,
                    dominanceMethod) {
             standardGeneric("dominanceAssigner")
           })
#' @export
#' @rdname dominanceAssigner
setMethod("dominanceAssigner",
          signature = c(inputData = "data.table",
                        dominanceCut = "numeric",
                        dominanceMethod = "character"),
          definition = function(inputData,
                                dominanceCut,
                                dominanceMethod){
            dominanceCut <- dominanceCut[dominanceCut!=1 & dominanceCut != 0] %>% unique %>% sort
            if(dominanceMethod == "BA"){
              inputData[,':='(totalBA = sum(IniBA)),
                        by = c("PlotID", "IniYear")]
              dominanceTable <- data.table(uniTreeID = character(), IniYear = numeric(),
                                           DominanceClass = character())
              plotlevels <- unique(inputData$PlotID)
              for(indiplot in plotlevels){
                plotdata <- inputData[PlotID == indiplot,]
                yearlevels <- unique(plotdata$IniYear)
                for(year in yearlevels){
                  plotyeardata <- plotdata[IniYear == year,]
                  cutpoint <- rep(nrow(plotyeardata),length(dominanceCut)+1) 
                  plotyeardata <- plotyeardata[order(IniBA, IniHeight),][,BArank:=1:nrow(plotyeardata)]
                  BAsum <- plotyeardata$IniBA[1]
                  k <- 1
                  for(i in 2:nrow(plotyeardata)){
                    prevBAsum <- BAsum
                    BAsum <- BAsum+plotyeardata$IniBA[i]
                    if(k <= length(dominanceCut)){
                      if(prevBAsum <= plotyeardata$totalBA[1]*dominanceCut[k] &
                         BAsum > plotyeardata$totalBA[1]*dominanceCut[k]){ 
                        cutpoint[k] <- i
                        k <- k+1
                      }
                    }
                  }
                  rm(k)
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
                                                                   labels = c(1:length(cutpoint)),
                                                                   include.lowest = TRUE))]
                  dominanceTable <- rbind(dominanceTable,plotyeardata[,.(uniTreeID, IniYear, DominanceClass)])
                }
              }
              inputData[,totalBA := NULL]
            } else if(dominanceMethod == "Biomass"){
              inputData[,':='(totalBiomass = sum(IniBiomass)),
                        by = c("PlotID", "IniYear")]
              dominanceTable <- data.table(uniTreeID = character(), IniYear = numeric(),
                                           DominanceClass = character())
              plotlevels <- unique(inputData$PlotID)
              for(indiplot in plotlevels){
                plotdata <- inputData[PlotID == indiplot,]
                yearlevels <- unique(plotdata$IniYear)
                for(year in yearlevels){
                  plotyeardata <- plotdata[IniYear == year,]
                  cutpoint <- rep(nrow(plotyeardata),length(dominanceCut)+1) 
                  plotyeardata <- plotyeardata[order(IniBiomass, IniHeight),][,Biomassrank:=1:nrow(plotyeardata)]
                  Biomasssum <- plotyeardata$IniBiomass[1]
                  k <- 1
                  for(i in 2:nrow(plotyeardata)){
                    prevBiomasssum <- Biomasssum
                    Biomasssum <- Biomasssum+plotyeardata$IniBiomass[i]
                    if(k <= length(dominanceCut)){
                      if(prevBiomasssum <= plotyeardata$totalBiomass[1]*dominanceCut[k] &
                         Biomasssum > plotyeardata$totalBiomass[1]*dominanceCut[k]){  
                        cutpoint[k] <- i
                      }
                    }
                    
                  }
                  rm(k)
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
                  plotyeardata[,DominanceClass := as.character(cut(Biomassrank, breaks = c(1, cutpoint),
                                                                   labels = c(1:length(cutpoint)),
                                                                   include.lowest = TRUE))]
                  dominanceTable <- rbind(dominanceTable,plotyeardata[,.(uniTreeID, IniYear, DominanceClass)])
                }
              }
              inputData[,totalBiomass := NULL]
            }
            dominanceTable[, DominanceClass:=as.numeric(as.character(DominanceClass))]
            setnames(dominanceTable, "DominanceClass", paste("Dominance_", dominanceMethod, sep = ""))
            analysesData <- setkey(inputData, uniTreeID, IniYear)[setkey(dominanceTable, uniTreeID, IniYear),
                                                                  nomatch = 0]
            rm(inputData)
            return(analysesData)
          })

