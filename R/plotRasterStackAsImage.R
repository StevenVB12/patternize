#' Plot rasterStack as image.
#'
#' @param rasterStack A single rasterStack.
#'
#' @export
#' @import raster

plotRasterstackAsImage <- function(rasterStack){

  x <- as.array(rasterStack[[1]])/255
  cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
  uniqueCols <- unique(cols)
  x2 <- match(cols, uniqueCols)
  dim(x2) <- dim(x)[1:2]
  raster::image(t(apply(x2, 2, rev)), col=uniqueCols, yaxt='n', xaxt='n')

}
