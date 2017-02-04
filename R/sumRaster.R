#' sum rasterlist for set of samples
#'
#' @param rList List of raster objects.
#' @param IDlist List of sample IDs.
#' @param type Type of rasterlist; simple or nested (result from RGB or k-means analysis, respectively)
#'
#' @examples
#' data(rasterList_lanRGB)
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' summedRaster <- sumRaster(rasterList_lanRGB, IDlist, type = 'RGB')
#'
#' data(rasterList_lanK)
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' summedRasterList <- sumRaster(rasterList_lanK, IDlist, type = 'k')
#'
#' @export

sumRaster <- function(rList, IDlist, type){

  subRasterList <- list()

  for(e in 1:length(IDlist)){

    subRasterList[[IDlist[e]]] <- rList[[IDlist[e]]]

  }

  if(type == 'RGB'){

    names(subRasterList) <- NULL
    subRasterList$fun <- sum
    subRasterList$na.rm <- TRUE
    summedRaster <- do.call(mosaic,subRasterList)

    return(summedRaster)

  }

  if(type == 'k'){ ## this doesn't work correctly yet!

    rasterListList <- list()

    for(n in 1:length(subRasterList)){

      sample <- subRasterList[[n]]

      for(e in 1:length(sample)){

        if(n == 1){

          rasterListList[[e]] <- c(sample[[e]])

        }

        else{

        rasterListList[[e]] <- c(rasterListList[[e]], sample[[e]])

        }
      }
    }

    summedRasterList <- list()

    for(k in 1:length(rasterListList)){

      names(rasterListList[[k]]) <- NULL
      rasterListList[[k]]$fun <- sum
      rasterListList[[k]]$na.rm <- TRUE
      summedRaster <- do.call(mosaic,rasterListList[[k]])

      summedRasterList[[k]] <- summedRaster

    }

    return(summedRasterList)

  }


}
