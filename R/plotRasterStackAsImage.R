#' Plot rasterStack as image.
#'
#' @param rasterStack A single rasterStack.
#' @param flipY Whether to flip the raster along the Y-axis.
#'
#' @export
#' @import raster

plotRasterstackAsImage <- function(rasterStack,
                                   flipY = FALSE){

  if(flipY){
    rasterStack <- raster::flip(rasterStack,'y')
  }

  x <- as.array(rasterStack)/255
  cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
  uniqueCols <- unique(cols)
  x2 <- match(cols, uniqueCols)
  dim(x2) <- dim(x)[1:2]
  raster::image(t(apply(x2, 2, rev)), col=uniqueCols, yaxt='n', xaxt='n')

}
