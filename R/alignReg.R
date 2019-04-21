#' Aligns images using \code{\link[RNiftyReg]{niftyreg}} utilities for automated image registration..
#'
#' @param sampleList List of RasterStack objects.
#' @param target Image imported as RasterStack used as target for registration.
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}} (default = NULL).
#' @param useBlockPercentage Block percentage as used in \code{\link[RNiftyReg]{niftyreg}}
#'    (default = 75).
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the
#'    original image.
#' @param removebgR Integer indicating the range RGB treshold to remove from image (e.g. 100 removes
#'    pixels with average RGB > 100; default = NULL) for registration analysis. This works only to
#'    remove a white background.
#' @param maskOutline When outline is specified, everything outside of the outline will be masked for
#'    the color extraction (default = NULL).
#' @param plotTransformed Whether to plot transformed images while processing (default = FALSE).
#'
#' @return List of raster objects.
#'
#'
#' @export
#' @import raster

alignReg <- function(sampleList,
                     target,
                     resampleFactor = NULL,
                     useBlockPercentage = 75,
                     crop = c(0,0,0,0),
                     removebgR = NULL,
                     maskOutline = NULL,
                     plotTransformed = FALSE){

  rasterList <- list()

  if(!identical(crop, c(0,0,0,0))){

    targetExtRaster <- crop
    target <- raster::crop(target, targetExtRaster)
  }

  if(!is.null(resampleFactor)){
    target <- redRes(target, resampleFactor)
  }

  targetA <- apply(raster::as.array(target), 1:2, mean)

  if(is.numeric(removebgR)){

    targetA <- apply(targetA, 1:2, function(x) ifelse(x > removebgR, 0, x))
  }

  for(n in 1:length(sampleList)){

    sStack <- sampleList[[n]]
    extRaster <- raster::extent(sStack)

    if(!identical(crop, c(0,0,0,0))){

      extRaster <- crop
      sStack <- crop(sStack, extRaster)
    }

    sourceRaster <- redRes(sStack, 1)

    if(!is.null(resampleFactor)){
      sourceRaster <- redRes(sStack, resampleFactor)
    }


    sourceRasterK <- sourceRaster

    sourceRaster <- apply(raster::as.array(sourceRaster), 1:2, mean)

    if(is.numeric(removebgR)){

      sourceRaster <- apply(sourceRaster, 1:2, function(x) ifelse(x > removebgR, 0, x))
    }

    result <- RNiftyReg::niftyreg(sourceRaster, targetA, useBlockPercentage=useBlockPercentage)

    transformedMap <- RNiftyReg::applyTransform(RNiftyReg::forward(result), raster::as.array(sourceRasterK), interpolation=0)
    r1 <- raster::raster(transformedMap[1:nrow(transformedMap),ncol(transformedMap):1,1])
    r2 <- raster::raster(transformedMap[1:nrow(transformedMap),ncol(transformedMap):1,2])
    r3 <- raster::raster(transformedMap[1:nrow(transformedMap),ncol(transformedMap):1,3])


    transRaster <-raster::stack(r1,r2,r3)
    transRaster <- raster::flip(transRaster,'x')

    raster::extent(transRaster) <- raster::extent(sourceRasterK)


    if(!is.null(maskOutline)){

      transRaster <- maskOutline(transRaster, maskOutline, refShape = 'target', flipOutline = 'y', crop = crop,
                                 imageList = sampleList)
    }
    # transRaster[transRaster == 0] <- NA




    if(!identical(raster::extent(transRaster), raster::extent(target))){
      raster::extent(transRaster) <- raster::extent(target)
    }

    # transRaster <- raster::flip(transRaster,'y')

    if(plotTransformed){

      x <- as.array(transRaster)/255
      cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
      uniqueCols <- unique(cols)
      x2 <- match(cols, uniqueCols)
      dim(x2) <- dim(x)[1:2]
      raster::image(t(apply(x2, 2, rev)), col=uniqueCols, yaxt='n', xaxt='n')

    }

    rasterList[[names(sampleList)[n]]] <- transRaster

    print(paste('sample', names(sampleList)[n], 'done and added to imageList', sep=' '))
  }

  return(rasterList)
}
