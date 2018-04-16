#' Create a target image (RasterStack) from a polygon.
#'
#' @param outline xy coordinates that define outline.
#' @param image Image imported as RasterStack used in the analysis. This is used to extract
#'    the extent and dimensions for the raster layers.
#' @param res Resolution for RasterStack (default = 300).
#' @param colorFill Color for the fill of the polygon (default = 'black').
#' @param colorBG Color for the background (default = 'white').
#' @param sigma Size of sigma for Gaussian blurring (default = 10).
#' @param plot Whether to plot the created target image (default = FALSE).
#'
#' @return RasterStack
#'
#' @examples
#'
#' \dontrun{
#' outline_BC0077 <- read.table(paste(system.file("extdata",  package = 'patternize'),
#' '/BC0077_outline.txt', sep=''), header = FALSE)
#'
#' data(imageList)
#'
#' target <- createTarget(outline_BC0077, imageList[[1]], plot =  TRUE)
#' }
#'
#' @export
#' @import raster
#' @importFrom grDevices col2rgb rgb

createTarget <- function(outline,
                         image,
                         res = 300,
                         colorFill = 'black',
                         colorBG = 'white',
                         sigma = 10,
                         plot = FALSE){

  if(is.character(colorFill)){
    colorFill <- col2rgb(colorFill)
  }

  if(is.character(colorBG)){
    colorBG <- col2rgb(colorBG)
  }

  rasterEx <- raster::extent(image)

  outline[,2] <- rasterEx[4] - outline[,2]

  poly <- sp::Polygons(list(sp::Polygon(outline)),paste("r"))

  polyList  <- c(poly)
  polyNames <- c(paste("r"))
  sr=sp::SpatialPolygons(polyList)
  srdf=sp::SpatialPolygonsDataFrame(sr, data.frame(1:length(polyNames), row.names=polyNames))

  r <- raster::raster(rasterEx, nrow = res, ncol = res)

  print('making raster layers')
  rr1 <-raster::rasterize(srdf, r, colorFill[1], background = colorBG[1])
  print('rasterized layer 1/3')
  rr2 <-raster::rasterize(srdf, r, colorFill[2], background = colorBG[2])
  print('rasterized layer 2/3')
  rr3 <-raster::rasterize(srdf, r, colorFill[3], background = colorBG[3])
  print('rasterized layer 3/3')

  gf <- focalWeight(r, sigma, "Gauss")

  rrr1 <- raster::focal(rr1, gf)
  rrr2 <- raster::focal(rr2, gf)
  rrr3 <- raster::focal(rr3, gf)

  rr <- raster::stack(rrr1, rrr2, rrr3)
  rr[is.na(rr)] <- 255

  if(plot){
    print('making plot...')
    x <- as.array(rr)/255
    cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
    uniqueCols <- unique(cols)
    x2 <- match(cols, uniqueCols)
    dim(x2) <- dim(x)[1:2]
    raster::image(t(apply(x2,2,rev)), col=uniqueCols,yaxt='n', xaxt='n')
  }

  print('resampling raster')
  r2 <- raster::raster(rasterEx, nrow=dim(image)[1], ncol=dim(image)[2])
  rrr <- raster::resample(rr,r2, datatype="INT1U", method='ngb')

  print('done')

  return(rrr)
}
