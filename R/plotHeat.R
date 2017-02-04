#' Plot heatmap from summed rasterList
#'
#' @param summedRaster Summed raster or summedRasterList.
#' @param IDlist List of sample IDs.
#' @param colpalette Vector of colors for color palette (default = c("white","lightblue","blue","green", "yellow","red"))
#'
#'
#' @examples
#' data(rasterList_lanRGB)
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' summedRaster <- sumRaster(rasterList_lanRGB, IDlist, type = 'RGB')
#'
#' plotHeat(summedRaster, IDlist)
#'
#' data(rasterList_lanK)
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' summedRasterList <- sumRaster(rasterList_lanK, IDlist, type = 'k')
#'
#' plotHeat(summedRasterList, IDlist)
#'
#' plotHeat(summedRasterList[[2]], IDlist)
#'
#' @export
#' @import raster

plotHeat <- function(summedRaster, IDlist, colpalette = NULL){

  if(is.null(colpalette)){

    colfunc <- colorRampPalette(c("white","lightblue","blue","green", "yellow","red"))
  }

  else{

    if(!is.vector(colpalette)){

      stop('Specified color palette is not a vector')

    }

    colfunc <- colorRampPalette(colpalette)

  }

  if(!is.list(summedRaster)){

    par(mfrow=c(1,1), mai=c(0.05,0.8,0.05,0.8), oma=c(1,1,1,1)+1)

    plot(summedRaster/length(IDlist), col=colfunc(20), xaxt='n', yaxt='n', box=F, axes=F, zlim=c(0,1))

  }

  else{

    par(mfrow=c(2,round((length(summedRaster)+1)/2)), mai=c(0.05,0.8,0.05,0.8), oma=c(1,1,1,1)+1)

    for(k in 1:length(summedRaster)){

      plot(summedRaster[[k]]/length(IDlist), col=colfunc(20), xaxt='n', yaxt='n', box=F, axes=F, zlim=c(0,1))

    }
  }
}
