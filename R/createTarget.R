#' Create a target image (raster stack) from a polygon.
#'
#' @param outline xy coordinates that define outline.
#' @param image Image used in the analysis. This us used to extract the extant and dimensions for the raster layers.
#' @param color Color for the fill of the polygon (default = 'black'.
#' @param plot Zhether to plot the created target image.
#'
#' @examples
#' outline_BC0077 <- read.table(paste(system.file("extdata",  package = 'patternize'), '/BC0077_outline.txt', sep=''), h= F)
#' data(imageList)
#' target <- createTarget(outline_BC0077, imageList[[1]], plot =  TRUE)

createTarget <- function(outline, image, color = 'black', plot = FALSE){

  rasterEx <- raster::extent(image)

  poly <- sp::Polygons(list(Polygon(outline)),paste("r"))


  polyList  <- c(poly)
  polyNames <- c(paste("r"))
  sr=sp::SpatialPolygons(polyList)
  srdf=sp::SpatialPolygonsDataFrame(sr, data.frame(1:length(polyNames), row.names=polyNames))

  r <- raster::raster(rasterEx, nrow=dim(image)[1], ncol=dim(image)[2])

  rr1 <-raster::rasterize(srdf, r)
  rr2 <-raster::rasterize(srdf, r)
  rr3 <-raster::rasterize(srdf, r)

  if(is.character(color)){
    color <-  col2rgb(color)
  }

  rr1[rr1 == 1] <- color[1]
  rr2[rr2 == 1] <- color[2]
  rr3[rr3 == 1] <- color[3]

  rr <- raster::stack(rr1, rr2, rr3)
  rr[is.na(rr)] <- 255

  if(plot){
    x <- as.array(rr)/255
    cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
    uniqueCols <- unique(cols)
    x2 <- match(cols, uniqueCols)
    dim(x2) <- dim(x)[1:2]
    raster::image(t(apply(x2,2,rev)), col=uniqueCols,yaxt='n', xaxt='n')
  }

  return(rr)
}


