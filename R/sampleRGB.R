#' Interactive function to sample RGB value from pixel in an image.
#'
#' @param image Image imported as a RasterStack.
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the
#'    original image.
#'
#' @return RGB vector
#'
#' @examples
#' image <- raster::stack(system.file("extdata", "BC0077.jpg", package = "patternize"))
#' RGB <- sampleRGB(image, resampleFactor = 1)
#'
#' @export
#' @import raster
#' @importFrom imager as.cimg imsub

sampleRGB <- function(image,
                      resampleFactor = NULL,
                      crop = c(0,0,0,0)){

  # Reduce resolution
  if(!is.null(resampleFactor)){

    image <- redRes(image, resampleFactor)
  }

  # Crop image
  if(!identical(crop, c(0,0,0,0))){

    rasterExt <- crop
    image <- raster::crop(image, rasterExt)
  }

  # Transform to imager format
  imA <- raster::as.array(image)
  imR <- as.raster(imA, nrow = dim(image)[1], ncol = dim(image)[2], max = 255)
  im <- imager::as.cimg(imR)

  plot(im)

  # Pick pixel and return RGB
  print("Choose a point for which you want RGB values.")

  xy <- locator(n=1)

  print(paste('x: ', as.character(xy)[1], 'y: ', as.character(xy)[2]))

  RGB <- as.vector(imager::imsub(im,x = xy$x,y = xy$y))

  print(paste(c('RGB: ', RGB), collapse = ' '))

  return(RGB)
}
