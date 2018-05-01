#' Sample landmarks in an image.
#'
#' @param sampleList RasterStack or list of RasterStack objects as obtained
#'    by \code{\link{makeList}}.
#' @param resampleFactor Integer for downsampling the image(s) used by \code{\link{redRes}}.
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the
#'    original image.
#'
#' @return landmark matrix or landmark list
#'
#' @examples
#'
#' \dontrun{
#' IDlist <- c('BC0077','BC0071')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '.jpg'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#'
#' landmarkList <- sampleLandmarks(imageList)
#' }
#'
#' @export
#' @import raster
#' @importFrom imager as.cimg imsub
#' @importFrom graphics locator

sampleLandmarks <- function(sampleList,
                            resampleFactor = NULL,
                            crop = c(0,0,0,0)){

  objectList <- list()

  if(length(sampleList) > 1){

    for(n in 1:length(sampleList)){

      image <- sampleList[[n]]

      rasterExt <- raster::extent(image)

      # Reduce resolution
      if(!is.null(resampleFactor)){

        image <- redRes(image, resampleFactor)
        raster::extent(image) <- c(rasterExt[1]/resampleFactor,
                                   rasterExt[2]/resampleFactor,
                                   rasterExt[3]/resampleFactor,
                                   rasterExt[4]/resampleFactor)
      }

      if(is.null(resampleFactor)){
        resampleFactor <- 1
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

      print(paste('Choose landmarks for image', names(sampleList)[n], sep=' '))
      print('Click outside image area to continue.')

      xyCoords <- c()

      lN <- 0

      while(1){

        xy <- locator(n=1)

        if(as.numeric(xy)[1] > dim(im)[1] || as.numeric(xy)[1] < 0 || as.numeric(xy)[2] > dim(im)[1] || as.numeric(xy)[2] < 0){

          if(n != length(sampleList)){
            print('Loading next image...')
          }

          break
        }

        lN <- lN + 1

        xy$x <- xy$x * resampleFactor + crop[1]*resampleFactor
        xy$y <- xy$y * resampleFactor + crop[3]*resampleFactor

        print(paste('x', lN, ': ', as.character(xy)[1], '  y', lN, ': ', as.character(xy)[2], sep=''))

        xyCoords <- c(xyCoords, c(xy$x, xy$y))
      }

      landmarks <- matrix(xyCoords, ncol=2, byrow=T)

      objectList[[names(sampleList)[n]]] <- landmarks

    }
  }

  if(length(sampleList) == 1){

    image <- sampleList[[1]]

    rasterExt <- raster::extent(image)

    if(!is.null(resampleFactor)){

      image <- redRes(image, resampleFactor)
      raster::extent(image) <- c(rasterExt[1]/resampleFactor,
                                 rasterExt[2]/resampleFactor,
                                 rasterExt[3]/resampleFactor,
                                 rasterExt[4]/resampleFactor)
    }

    if(is.null(resampleFactor)){
      resampleFactor <- 1
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

    print(paste('Choose landmarks for image', names(sampleList)[1], sep=' '))
    print('Click outside image area to continue.')

    xyCoords <- c()

    lN <- 0

    while(1){

      xy <- locator(n=1)

      if(as.numeric(xy)[1] > dim(im)[1] || as.numeric(xy)[1] < 0 || as.numeric(xy)[2] > dim(im)[1] || as.numeric(xy)[2] < 0){
        break
      }

      lN <- lN + 1
      xy$x <- xy$x * resampleFactor + crop[1]*resampleFactor
      xy$y <- xy$y * resampleFactor + crop[3]*resampleFactor

      print(paste('x', lN, ': ', as.character(xy)[1], '  y', lN, ': ', as.character(xy)[2], sep=''))

      xyCoords <- c(xyCoords, c(xy$x, xy$y))
    }

    objectList <- matrix(xyCoords, ncol=2, byrow=T)
  }

  return(objectList)
}
