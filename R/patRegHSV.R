#' Aligns images using \code{\link[RNiftyReg]{niftyreg}} utilities for automated image registration
#' and extracts colors using a predefined HSV values and cutoff value.
#'
#' @param sampleList List of RasterStack objects.
#' @param target Image imported as RasterStack used as target for registration.
#' @param HSV Values for color pattern extraction specified as HSV vector.
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}} (default = NULL).
#' @param useBlockPercentage Block percentage as used in \code{\link[RNiftyReg]{niftyreg}}
#'    (default = 75).
#' @param colOffset Color offset for color pattern extraction (default = 0.10).
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the
#'    original image.
#' @param removebgR Integer indicating the range RGB treshold to remove from image (e.g. 100 removes
#'    pixels with average RGB > 100; default = NULL) for registration analysis. This works only to
#'    remove a white background.
#' @param maskOutline When outline is specified, everything outside of the outline will be masked for
#'    the color extraction (default = NULL).
#' @param plot Whether to plot transformed color patterns while processing (default = FALSE).
#'    Transformed color patterns can be plot on top of each other ('stack') or next to the
#'    original image for each sample ('compare').
#' @param focal Whether to perform Gaussian blurring (default = FALSE).
#' @param sigma Size of sigma for Gaussian blurring (default = 3).
#' @param iterations Number of iterations for recalculating average color (default = 0). If set, the
#'    RGB value for pattern extraction will be iteratively recalculated to be the average of the
#'    extracted area. This may improve extraction of distinct color pattern, but fail for more
#'    gradually distributed (in color space) patterns.
#' @param ignoreHSVvalue Whether to ignore the HSV value (~darkness).
#' @param patternsToFile Name of directory to which the color pattern of each individual will be
#'    outputted (default = NULL).
#'
#' @return List of raster objects.
#'
#' @examples
#' \dontrun{
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '.jpg'
#'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#'
#' target <- imageList[[1]]
#'
#' HSV <- c(0.025,1,0.45)
#'
#' # Note that this example only aligns one image with the target,
#' # remove [2] to run a full examples.
#' rasterList_regHSV <- patRegRGB(imageList[2], target, HSV,
#' colOffset= 0.15, crop = c(100,400,40,250), removebgR = 100, plot = 'stack')
#' }
#'
#' @export
#' @import raster
#' @importFrom grDevices hsv rgb2hsv dev.off png

patRegHSV <- function(sampleList,
                      target,
                      HSV,
                      resampleFactor = NULL,
                      useBlockPercentage = 75,
                      colOffset=0.10,
                      crop = c(0,0,0,0),
                      removebgR = NULL,
                      maskOutline = NULL,
                      plot = FALSE,
                      focal =  FALSE,
                      sigma = 3,
                      iterations = 0,
                      ignoreHSVvalue = FALSE,
                      patternsToFile = NULL){

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

    if(focal){

      gf <- focalWeight(sourceRaster, sigma, "Gauss")

      rrr1 <- raster::focal(sourceRaster[[1]], gf)
      rrr2 <- raster::focal(sourceRaster[[2]], gf)
      rrr3 <- raster::focal(sourceRaster[[3]], gf)

      sourceRaster <- raster::stack(rrr1, rrr2, rrr3)
    }

    sourceRasterK <- sourceRaster

    sourceRaster <- apply(raster::as.array(sourceRaster), 1:2, mean)

    if(is.numeric(removebgR)){

      sourceRaster <- apply(sourceRaster, 1:2, function(x) ifelse(x > removebgR, 0, x))
    }

    result <- RNiftyReg::niftyreg(sourceRaster, targetA, useBlockPercentage=useBlockPercentage)

    # convert RGB values to HSV values
    sourceRasterK <- raster::overlay(sourceRasterK, fun = rgb2hsv)

    if(ignoreHSVvalue == TRUE){
      map <- apply(raster::as.array(sourceRasterK), 1:2, function(x) all(abs(x[1:2]-HSV[1:2]) < colOffset))
    }
    else{
      map <- apply(raster::as.array(sourceRasterK), 1:2, function(x) all(abs(x-HSV) < colOffset))
    }

    if(all(map == FALSE)){
      warning("The RGB range does not seem to overlap with any of the RGB values in the image")
    }

    if(iterations > 0){
      if(all(map == FALSE)){
        warning("Iterations can't be performed")
      }
    }

    if(!all(map == FALSE)){

      x <- 1
      while(x <= iterations){
        x <- x + 1

        mapRaster <- raster::raster(as.matrix(map))
        extent(mapRaster) <- extRaster
        mapRaster[mapRaster == 0] <- NA

        mapMASK<-raster::mask(sourceRasterK, mapRaster)

        RGBnew <- c(mean(na.omit(as.data.frame(mapMASK[[1]]))[,1]),
                    mean(na.omit(as.data.frame(mapMASK[[2]]))[,1]),
                    mean(na.omit(as.data.frame(mapMASK[[3]]))[,1]))

        if(ignoreHSVvalue == TRUE){
          map <- apply(raster::as.array(sourceRasterK), 1:2, function(x) all(abs(x[1:2]-HSV[1:2]) < colOffset))
        }
        else{
          map <- apply(raster::as.array(sourceRasterK), 1:2, function(x) all(abs(x-HSV) < colOffset))
        }
      }

      transformedMap <- RNiftyReg::applyTransform(RNiftyReg::forward(result), map, interpolation=0)
      transformedMapMatrix <- transformedMap[1:nrow(transformedMap),ncol(transformedMap):1]

      transRaster <- raster::raster(transformedMapMatrix)
      raster::extent(transRaster) <- extRaster

      if(!is.null(maskOutline)){

        transRaster <- maskOutline(transRaster, maskOutline, refShape = 'target', flipOutline = 'y', crop = crop,
                                   imageList = sampleList)
      }
      transRaster[transRaster == 0] <- NA

    }

    else{
      transRaster <- raster::raster(extRaster, nrow=dim(sStack)[1], ncol=dim(sStack)[2], vals = rep(NA, dim(sStack)[1]*dim(sStack)[2]))
    }

    if(!identical(raster::extent(transRaster), raster::extent(target))){
      raster::extent(transRaster) <- raster::extent(target)
    }

    if(plot == 'stack'){

      par(mfrow=c(1,1))
      if(n == 1){
        plot(1, type="n", axes = FALSE, xlab='', ylab='')
      }

      par(new = TRUE)
      raster::plot(transRaster, col=rgb(1,0,0,alpha=1/length(sampleList)), legend = FALSE)
    }

    if(plot == 'compare'){

      par(mfrow=c(1,2))
      plot(1, type="n", xlab='', ylab='', xaxt='n', yaxt='n', axes= FALSE, bty='n')
      par(new = TRUE)
      plot(raster::flip(transRaster,'x'), col='black', legend = FALSE, xaxt='n', yaxt='n', axes= FALSE, bty='n')

      x <- as.array(sStack)
      cols <- hsv(x[,,1], x[,,2], x[,,3])
      uniqueCols <- unique(cols)
      x2 <- match(cols, uniqueCols)
      dim(x2) <- dim(x)[1:2]
      raster::image(apply(x2, 1, rev), col=uniqueCols, yaxt='n', xaxt='n')

    }

    if(!is.null(patternsToFile)){

      dir.create(file.path(patternsToFile), showWarnings = FALSE)

      png(paste(patternsToFile, '/', names(sampleList)[n], '.png', sep=''))

      plot(1, type="n", axes = FALSE, xlab='', ylab='', main = names(sampleList)[n])
      par(new = TRUE)
      raster::plot(transRaster, col=rgb(1,0,0,alpha=1), legend = FALSE)

      dev.off()
    }

    rasterList[[names(sampleList)[n]]] <- transRaster

    print(paste('sample', names(sampleList)[n], 'done and added to rasterList', sep=' '))
  }

  return(rasterList)
}
