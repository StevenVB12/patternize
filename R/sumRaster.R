#' This function sums the individual color pattern RasterLayes as obtained by the main patternize
#' functions.
#'
#' @param rList List of RasterLayers or list of RasterLayers for each k-means cluster.
#' @param IDlist List of sample IDs.
#' @param type Type of rasterlist; 'RGB' or 'k' (result from RGB or k-means analysis, respectively).
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

sumRaster <- function(rList,
                      IDlist,
                      type){

  subRasterList <- list()

  for(e in 1:length(IDlist)){

    subRasterList[[IDlist[e]]] <- rList[[IDlist[e]]]

  }

  if(type == 'RGB'){

    for(e in 1:length(IDlist)){

      if(!identical(raster::extent(subRasterList[[IDlist[e]]]), raster::extent(subRasterList[[1]]))){

        print(paste('raster extent set to extent of rasterLayer 1 in IDlist for sample', IDlist[e], sep=' '))
        raster::extent(subRasterList[[IDlist[e]]]) <- raster::extent(subRasterList[[1]])
      }
    }

    names(subRasterList) <- NULL
    subRasterList$fun <- sum
    subRasterList$na.rm <- TRUE
    summedRaster <- do.call(mosaic,subRasterList)

    return(summedRaster)

  }

  if(type == 'k'){

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

      for(e in 1:length(IDlist)){

        if(!identical(raster::extent(rasterListList[[k]][[e]]), raster::extent(rasterListList[[k]][[1]]))){

          print(paste('raster extent set to extent of rasterLayer 1 in IDlist for cluster', k , 'in sample', IDlist[e], sep=' '))
          raster::extent(rasterListList[[k]][[e]]) <- raster::extent(rasterListList[[k]][[1]])
        }
      }

      names(rasterListList[[k]]) <- NULL
      rasterListList[[k]]$fun <- sum
      rasterListList[[k]]$na.rm <- TRUE
      summedRaster <- do.call(mosaic,rasterListList[[k]])

      summedRasterList[[k]] <- summedRaster

    }

    return(summedRasterList)

  }


}
