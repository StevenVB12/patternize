#' Interactive function to sample RGB value from pixel or square area in an image.
#'
#' @param image Image imported as a RasterStack.
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the
#'    original image.
#' @param type Set 'point' to extract RGB from a single point or 'area' to extract from a square
#'    area defined by setting two points (default = 'point').
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
#' @importFrom graphics locator

sampleRGB <- function(image,
                      resampleFactor = NULL,
                      crop = c(0,0,0,0),
                      type = 'point'){

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
  imA[is.na(imA[])] <- 0
  imR <- as.raster(imA, nrow = dim(image)[1], ncol = dim(image)[2], max = 255)
  im <- imager::as.cimg(imR)

  plot(im)

  if(type == 'point'){

    # Pick pixel and return RGB
    print("Choose a point for which you want RGB values.")

    xy <- locator(n=1)

    print(paste('x: ', as.character(xy)[1], 'y: ', as.character(xy)[2]))

    RGB <- as.vector(imager::imsub(im,x = xy$x,y = xy$y))

    print(paste(c('RGB: ', RGB), collapse = ' '))
  }

  if(type == 'area'){

    print("Choose two points to define square area for which you want RGB values.")

    xy1 <- locator(n=1)

    print(paste('x: ', as.character(xy1)[1], 'y: ', as.character(xy1)[2]))

    xy2 <- locator(n=1)

    print(paste('x: ', as.character(xy2)[1], 'y: ', as.character(xy2)[2]))

    xy <- as.matrix(rbind(as.numeric(as.character(xy1)),as.numeric(as.character(xy2))))

    minX <- min(xy[,1])
    maxX <- max(xy[,1])
    minY <- min(xy[,2])
    maxY <- max(xy[,2])

    x <- 0 # Otherwis R CMD check returns note for imsub function
    y <- 0

    imS <- as.array(imager::imsub(im, x %in% c(round(minX,0):round(maxX,0)), y %in% c(round(minY,0):round(maxY,0))))

    R <- mean(imS[,,1,1])
    G <- mean(imS[,,1,2])
    B <- mean(imS[,,1,3])

    RGB <- c(R,G,B)

    print(paste(c('RGB: ', RGB), collapse = ' '))

  }

  return(RGB)
}
