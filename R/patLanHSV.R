#' Aligns images usings transformations obtained from fixed landmarks and extracts colors
#' using a predefined RGB values and cutoff value.
#'
#' @param sampleList List of RasterStack objects.
#' @param landList Landmark list as returned by \code{\link[patternize]{makeList}}.
#' @param HSV HSV values for color pattern extraction specified as vector.
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param colOffset Color offset for color pattern extraction (default = 0.10).
#' @param crop Whether to use the landmarks range to crop the image. This can speed up the
#'    analysis (default = FALSE).
#' @param cropOffset Vector c(xmin, xmax, ymin, ymax) that specifies the number of pixels you
#'    want the cropping to be offset from the landmarks (in case the landmarks do not surround
#'    the entire color pattern). The values specified should present the percentage of the maximum
#'    landmark value along the x and y axis.
#' @param res Resolution for color pattern raster (default = 300). This should be reduced if
#'    the number of pixels in the image is lower than th raster.
#' @param transformRef ID of reference sample for shape to which color patterns will be transformed
#'    to. Can be 'meanshape' for transforming to mean shape of Procrustes analysis.
#' @param transformType Transformation type as used by \code{\link[Morpho]{computeTransform}}
#'    (default ='tps').
#' @param adjustCoords Adjust landmark coordinates in case they are reversed compared to pixel
#'    coordinates (default = FALSE).
#' @param plot Whether to plot transformed color patterns while processing (default = NULL).
#'    Transformed color patterns can be plot on top of each other ('stack') or next to the
#'    original image for each sample ('compare').
#' @param focal Whether to perform Gaussian blurring (default = FALSE).
#' @param sigma Size of sigma for Gaussian blurring (default = 3).
#' @param iterations Number of iterations for recalculating average color.
#' @param ignoreHSVvalue Whether to ignore the HSV value (~darkness).
#' @param patternsToFile Name of directory to which the color pattern of each individual will be
#'    outputted (default = NULL).
#'
#' @return  List of raster objects.
#'
#' @examples
#'
#' \dontrun{
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#'
#' landmarkList <- makeList(IDlist, 'landmark', prepath, extension)
#'
#' extension <- '.jpg'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#'
#' HSV <- c(0.025,1,0.45)
#' rasterList_lanHSV <- patLanRGB(imageList, landmarkList, HSV,
#' colOffset = 0.15, crop = TRUE, res = 100, adjustCoords = TRUE, plot = 'stack')
#' }
#'
#' @export
#' @import raster
#' @importFrom utils capture.output
#' @importFrom stats na.omit
#' @importFrom Morpho procSym computeTransform applyTransform
#' @importFrom grDevices hsv rgb2hsv dev.off png

patLanHSV <- function(sampleList,
                      landList,
                      HSV,
                      resampleFactor = NULL,
                      colOffset = 0.10,
                      crop = FALSE,
                      cropOffset = c(0,0,0,0),
                      res = 300,
                      transformRef = 'meanshape',
                      transformType = 'tps',
                      adjustCoords = FALSE,
                      plot = NULL,
                      focal =  FALSE,
                      sigma = 3,
                      iterations = 0,
                      ignoreHSVvalue = FALSE,
                      patternsToFile = NULL){

  rasterList <- list()

  if(length(sampleList) != length(landList)){
    stop("sampleList is not of the same length as lanArray")
  }

  for(n in 1:length(sampleList)){
    if(names(sampleList)[n] != names(landList)[n]){
      stop("samples are not in the same order in sampleList and lanArray")
    }
  }

  lanArray <- lanArray(landList, adjustCoords, sampleList)

  if(is.matrix(transformRef)){

    refShape <- transformRef

  }

  if(!is.matrix(transformRef)){

    if(transformRef == 'meanshape'){

      invisible(capture.output(transformed <- Morpho::procSym(lanArray)))
      refShape <- transformed$mshape

    }

    if(transformRef %in% names(landList)){

      e <- which(names(landList) == transformRef)
      refShape <- lanArray[,,e]
    }
  }


  for(n in 1:length(sampleList)){

    image <- sampleList[[n]]
    extRasterOr <- raster::extent(image)

    if(!is.null(resampleFactor)){
      image <- redRes(image, resampleFactor)
    }

    # convert RGB values to HSV values
    image <- raster::overlay(image, fun = rgb2hsv)

    if(crop){

      landm <- lanArray[,,n]

      extRaster <- raster::extent(min(landm[,1])-min(landm[,1])*cropOffset[1]/100,
                                  max(landm[,1])+max(landm[,1])*cropOffset[2]/100,
                                  min(landm[,2])-min(landm[,2])*cropOffset[3]/100,
                                  max(landm[,2])+max(landm[,2])*cropOffset[4]/100)

      imageC <- raster::crop(image, extRaster)

      y <- raster::raster(ncol = dim(image)[2], nrow = dim(image)[1])
      extent(y) <- extRasterOr
      image <- resample(imageC, y)
    }

    if(focal){
      gf <- focalWeight(image, sigma, "Gauss")

      rrr1 <- raster::focal(image[[1]], gf)
      rrr2 <- raster::focal(image[[2]], gf)
      rrr3 <- raster::focal(image[[3]], gf)

      image <- raster::stack(rrr1, rrr2, rrr3)
    }

    if(ignoreHSVvalue == TRUE){
      map <- apply(raster::as.array(image), 1:2, function(x) all(abs(x[1:2]-HSV[1:2]) < colOffset))
    }
    else{
      map <- apply(raster::as.array(image), 1:2, function(x) all(abs(x-HSV) < colOffset))
    }

    if(all(map == FALSE)){
      warning("The HSV range does not seem to overlap with any of the RGB values in the image")
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
        raster::extent(mapRaster) <- extRasterOr
        mapRaster[mapRaster == 0] <- NA

        mapMASK<-raster::mask(image, mapRaster)

        HSVnew <- c(mean(na.omit(as.data.frame(mapMASK[[1]]))[,1]),
                    mean(na.omit(as.data.frame(mapMASK[[2]]))[,1]),
                    mean(na.omit(as.data.frame(mapMASK[[3]]))[,1]))

        if(ignoreHSVvalue == TRUE){
          map <- apply(raster::as.array(image), 1:2, function(x) all(abs(x[1:2]-HSVnew[1:2]) < colOffset))
        }
        else{
          map <- apply(raster::as.array(image), 1:2, function(x) all(abs(x-HSVnew) < colOffset))
        }
      }

      mapR <- raster::raster(map)
      raster::extent(mapR) <- extRasterOr

      mapDF <- raster::as.data.frame(mapR, xy = TRUE)

      mapDFs <- subset(mapDF, mapDF$layer == TRUE)

      invisible(capture.output(transMatrix <- Morpho::computeTransform(refShape, as.matrix(lanArray[,,n]), type = transformType)))

      invisible(capture.output(mapTransformed <- Morpho::applyTransform(as.matrix(mapDFs[1:2]), transMatrix)))

      r <- raster::raster(ncol = res, nrow = res)

      if(any(c(transformRef == 'meanshape', is.matrix(transformRef)))){

        raster::extent(r) <- raster::extent(min(refShape[,2])-3*max(refShape[,1])*cropOffset[3]/100,
                                            max(refShape[,2])+3*max(refShape[,1])*cropOffset[4]/100,
                                            min(refShape[,1])-3*max(refShape[,2])*cropOffset[1]/100,
                                            max(refShape[,1])+3*max(refShape[,2])*cropOffset[2]/100)
      }

      else{

        raster::extent(r) <- raster::extent(min(refShape[,1])-max(refShape[,1])*cropOffset[1]/100,
                                            max(refShape[,1])+max(refShape[,1])*cropOffset[2]/100,
                                            min(refShape[,2])-max(refShape[,2])*cropOffset[3]/100,
                                            max(refShape[,2])+max(refShape[,2])*cropOffset[4]/100)
      }

      patternRaster <- raster::rasterize(mapTransformed, field = 1, r)

    }

    else{

      if(any(c(transformRef == 'meanshape', is.matrix(transformRef)))){

        patternRaster <- raster::extent(min(refShape[,1])-3*max(refShape[,1])*cropOffset[3]/100,
                                        max(refShape[,1])+3*max(refShape[,1])*cropOffset[4]/100,
                                        min(refShape[,2])-3*max(refShape[,2])*cropOffset[1]/100,
                                        max(refShape[,2])+3*max(refShape[,2])*cropOffset[2]/100)
      }

      else{

        patternRaster <- raster::raster(extent(min(refShape[,1])-max(refShape[,1])*cropOffset[1]/100,
                                               max(refShape[,1])+max(refShape[,1])*cropOffset[2]/100,
                                               min(refShape[,2])-max(refShape[,2])*cropOffset[3]/100,
                                               max(refShape[,2])+max(refShape[,2])*cropOffset[4]/100),
                                        ncol = res, nrow = res, vals = rep(NA, res*res))
      }
    }

    if(plot == 'stack'){

      par(mfrow=c(1,1))
      if(n == 1){
        plot(1, type="n", xlab='', ylab='', xaxt='n', yaxt='n', axes= FALSE, bty='n')
      }

      par(new = TRUE)

      if(any(c(transformRef == 'meanshape', is.matrix(transformRef)))){

        patternRasterP <- raster::flip(t(patternRaster), 'y')
      }

      else{
        patternRasterP <- patternRaster
      }

      plot(patternRasterP, col=rgb(1,0,0,alpha=1/length(sampleList)), legend = FALSE, xaxt='n', yaxt='n', axes= FALSE, bty='n')
    }

    if(plot == 'compare'){

      landm <- lanArray[,,n]

      par(mfrow=c(1,2))
      plot(1, type="n", xlab='', ylab='', xaxt='n', yaxt='n', axes= FALSE, bty='n')
      par(new = TRUE)
      plot(patternRaster, col='black', legend = FALSE, xaxt='n', yaxt='n', axes= FALSE, bty='n')

      rasterExt <- raster::extent(min(landm[,1])-min(landm[,1])*cropOffset[1]/100,
                                  max(landm[,1])+max(landm[,1])*cropOffset[2]/100,
                                  min(landm[,2])-min(landm[,2])*cropOffset[3]/100,
                                  max(landm[,2])+max(landm[,2])*cropOffset[4]/100)

      image <- raster::crop(image,rasterExt)

      image[is.na(image)] <- 1
      x <- as.array(image)
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
      raster::plot(patternRaster, col=rgb(1,0,0,alpha=1), legend = FALSE)

      dev.off()
    }

    rasterList[[names(landList)[n]]] <- patternRaster

    print(paste('sample', names(landList)[n], 'done and added to rasterList', sep=' '))
  }

  return(rasterList)
}
