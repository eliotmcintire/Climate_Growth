#' the function to calculate hegyi competition index
#'
#' @param data data.table which must have spatial comfiguration rows
#' 
#' @maxRadius numeric, the competition index will been calculated within this radius
#' 
#'
#' @return a data table that has three columns, i.e., active, mapcode and ecoregion
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
#' @rdname HeghyiCICalculation
#'
#' @author Yong Luo
#'
setGeneric("HeghyiCICalculation",
           function(data,
                    maxRadius) {
             standardGeneric("HeghyiCICalculation")
           })

#' @export
#' @rdname HeghyiCICalculation
setMethod(
  "HeghyiCICalculation",
  signature = c(data = "data.table",
                maxRadius = "numeric"),
  definition = function(data,
                        maxRadius){
    # calcuate coordination of each tree
    data[, ':='(coordX = sin(Angle*pi/180)*Distance,
                coordY = cos(Angle*pi/180)*Distance)]
    set(data, , c("Hegyi", "IntraHegyiRatio"), 0)
    for(i in 1:nrow(data)){
      # surrounding trees
      surroundingTrees <- data[PlotID == data$PlotID[i] & Year == data$Year[i] &
                                 TreeNumber != data$TreeNumber[i],]
      surroundingTrees[,XYDistance:=((coordX-data$coordX[i])^2+(coordY-data$coordY[i])^2)^0.5]
      surroundingTrees[XYDistance<0.1, XYDistance:=0.1]
      surroundingTrees <- surroundingTrees[XYDistance <= maxRadius,]
      surroundingTrees_SameSpecies <- surroundingTrees[Species == data$Species[i],]
      # tempHeigyi did not account for area correction
      tempHegyi <- sum(surroundingTrees$DBH/(data$DBH[i]*surroundingTrees$XYDistance))
      if(nrow(surroundingTrees_SameSpecies) > 0){
        tempIntraHegyi <- sum(surroundingTrees_SameSpecies$DBH/(data$DBH[i]*surroundingTrees_SameSpecies$XYDistance))
      } else {
        tempIntraHegyi <- 0
      }
      d <- data$Distance[i]
      A <- 2*(maxRadius^2)*acos(0.5*d/maxRadius)-0.5*d*((4*(maxRadius^2)-d^2)^0.5)
      prop <- A/500.3439
      data$Hegyi[i]<-tempHegyi/prop
      data$IntraHegyiRatio[i]<-tempIntraHegyi/tempHegyi
      if(i/1000==as.integer(i/1000))
      { cat("\n  The ",i,"th Row has been processed. \n\n") }
    }
    return(data)  
  })
