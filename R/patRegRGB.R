#' Image registration with RGB color extraction.
#'
#' @param sampleList List of RasterStack objects.
#' @param target Image used as target for registration.
#' @param RGB RGB values for color pattern extraction specified as vector.
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param useBlockPercentage Block percentage as used in NiftyReg.
#' @param colOffset Color offset for color pattern extraction.
#' @param removebg Whether to remove white background for image registration (default = FALSE).
#' @param plot Whether to plot transformed color patterns while processing (default = FALSE).
#'
#' @return List of raster objects.
#'
#' @examples
#' data(pictures)
#' target <- pictures[[1]]
#' RGB <- c(114,17,0)
#' rasterList <- patReg(pictures, target, RGB, colOffset= 0.15, removebg = TRUE, plot = TRUE)
#'
#' @export

patRegRGB <- function(sampleList, target, RGB, resampleFactor = 1, useBlockPercentage = 75, colOffset=0, removebg = FALSE, plot = FALSE){

  rasterList <- list()

  target <- redRes(target, resampleFactor)
  target <- raster::as.array(target)
  target <- apply(target,1:2,mean)

  if(removebg){
    target <- apply(target, 1:2, function(x) ifelse(x>100,0, x))
  }

  for(n in 1:length(sampleList)){

    sStack <- sampleList[[n]]

    extRaster <- raster::extent(sStack)
    # sourcePicture <- crop(sourcePicture, extRaster)
    source <- redRes(sStack, resampleFactor)
    source <- raster::as.array(source/255)
    sourceR <- apply(source,1:2,mean)

    if(removebg){
      sourceR <- apply(sourceR, 1:2, function(x) ifelse(x>100,0, x))
    }

    result <- RNiftyReg::niftyreg(sourceR, target, useBlockPercentage=useBlockPercentage)

    map <- apply(source, 1:2, function(x) all(abs(x-RGB) < colOffset*255))

    transformedMap <- RNiftyReg::applyTransform(RNiftyReg::forward(result), map, interpolation=0)
    transformedMapMatrix <- transformedMap[1:nrow(transformedMap),ncol(transformedMap):1]
    transformedMapMatrix[transformedMapMatrix == 0] <- NA

    r <- raster::raster(transformedMapMatrix)
    # extRaster <- extent(sourcePicture)etPicture[,,1]
    raster::extent(r) <- extRaster

    if(plot){

      if(n == 1){
        plot(1, type="n", axes=F, xlab='', ylab='')
      }

      par(new=T)
      plot(r, col=rgb(1,0,0,alpha=1/length(sampleList)),legend=F)
    }

    rasterList <- c(rasterList, r)
  }

  #   # Sum raster list
  #   rasterList$fun <- sum
  #   rasterList$na.rm <- TRUE
  #   rrr <- do.call(mosaic,rasterList)

  return(rasterList)
}
