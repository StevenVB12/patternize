#' Performs PCA
#'
#' @param rList List of raster objects.
#' @param popList List of vectors including sampleIDs for eacht population.
#' @param colList List of colors for each population.
#'
#' @examples
#' data(rasterList_lanRGB)
#'
#' pop1 <- c('BC0077','BC0071')
#' pop2 <- c('BC0050','BC0049','BC0004')
#' popList <- list(pop1, pop2)
#' colList <- c("red", "blue")
#'
#' pcaOut <- patPCA(rasterList_lanRGB, popList, colList)
#' comp <- prcomp(pcaOut[[1]])
#' plot(comp$rotation, col=pcaOut[[2]]$col, pch=19)
#'
#' @export


patPCA <- function(rList, popList, colList){

  for(r in 1:length(rList)){

    rList[[r]][is.na(rList[[r]][])] <- 0
    ras <- as.data.frame(rList[[r]])
    colnames(ras) <- names(rList)[[r]]

    if(r == 1){
      rasDF <- ras
    }
    else{
      rasDF <- cbind(rasDF, ras)
    }
  }

  for(p in 1:length(popList)){

    for(ind in length(popList[[p]])){

      if(p == 1){

        groupCol <- c(popList[[p]][ind], colList[p])
      }

      else{

        groupCol <- rbind(groupCol, c(popList[[p]][ind], colList[p]))
      }
    }
  }

  groupCol<-as.data.frame(groupCol)
  colnames(groupCol) <- c('sampleID', 'col')

  return(list(rasDF, groupCol))
}




