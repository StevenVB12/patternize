#' Aligns images usings transformations obtained from fixed landmarks and extracts colors using
#' k-means clustering.
#'
#' @param sampleList List of RasterStack objects.
#' @param landList Landmark list as returned by \code{\link[patternize]{makeList}}.
#' @param k Integere for defining number of k-means clusters (default = 3).
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param crop Whether to use the landmarks range to crop the image. This can significantly speed
#'    up the analysis (default = FALSE).
#' @param cropOffset Vector c(xmin, xmax, ymin, ymax) that specifies the number of pixels you
#'    want the cropping to be offset from the landmarks (in case the landmarks do not surround
#'    the entire color pattern). The values specified should present the percentage of the maximum
#'    landmark value along the x and y axis.
#' @param res Resolution for color pattern raster (default = 300). This should be reduced if the
#'    number of pixels in the image is lower than th raster.
#' @param transformRef ID of reference sample for shape to which color patterns will be transformed
#'    to. Can be 'meanshape' for transforming to mean shape of Procrustes analysis.
#' @param transformType Transformation type as used by \code{\link[Morpho]{computeTransform}}
#'    (default ='tps').
#' @param removebgK Integer indicating the range RGB treshold to remove from image (e.g. 100
#'    removes pixels with average RGB > 100; default = NULL) for k-means analysis. This works only
#'    to remove a white background.
#' @param adjustCoords Adjust landmark coordinates in case they are reversed compared to pixel
#'    coordinates (default = FALSE).
#' @param plot Whether to plot transformed color patterns while processing (default = FALSE).
#' @param focal Whether to perform Gaussian blurring (default = FALSE).
#' @param sigma Size of sigma for Gaussian blurring (default = 3).
#'
#' @return  List of summed raster for each k-means cluster objects.
#'
#' @examples
#' \dontrun{
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#' landmarkList <- makeList(IDlist, 'landmark', prepath, extension)
#'
#' extension <- '.jpg'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#' # Note that this example only aligns two images with the target,
#' # remove [1:2] to run a full examples.
#' rasterList_lanK <- patLanK(imageList[1:2], landmarkList[1:2], k = 4, crop = TRUE,
#' res = 100, removebgK = 100, adjustCoords = TRUE, plot = TRUE)
#' }
#'
#' @export
#' @import raster
#' @importFrom utils capture.output
#' @importFrom Morpho procSym computeTransform applyTransform

patLanK <- function(sampleList,
                    landList,
                    k = 3,
                    resampleFactor = NULL,
                    crop = FALSE,
                    cropOffset = c(0,0,0,0),
                    res = 300,
                    transformRef = 'meanshape',
                    transformType='tps',
                    removebgK = NULL,
                    adjustCoords = FALSE,
                    plot = FALSE,
                    focal =  FALSE,
                    sigma = 3){

  rasterList <- list()

  if(length(sampleList) != length(landList)){
    stop("imageList is not of the same length as lanArray")
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

  if(transformRef == 'meanshape'){

    invisible(capture.output(transformed <- Morpho::procSym(lanArray)))
    refShape <- transformed$mshape

  }

  if(transformRef %in% names(landList)){

    e <- which(names(landList) == transformRef)
    refShape <- lanArray[e]
  }


  for(n in 1:length(sampleList)){

    image <- sampleList[[n]]
    extRasterOr <- raster::extent(image)

    if(!is.null(resampleFactor)){
      image <- redRes(image, resampleFactor)
    }

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

    if(is.vector(removebgK)){

      toMask <- apply(raster::as.array(image), 1:2, function(x) all(x > removebgK))

      toMaskR <- raster::raster(as.matrix(toMask))
      raster::extent(toMaskR) <- raster::extent(image)
      toMaskR[toMaskR == 0] <- NA

      image<-raster::mask(image, toMaskR, inverse = T)
      image[is.na(image)] <- 0
    }

    # k-means clustering of image

    if(n==1){
      startCenter = NULL
    }

    else{
      startCenter <- K$centers
    }

    image[is.na(image)] <- 255
    imageKmeans <- kImage(raster::as.array(image), k, startCenter)

    image.segmented <- imageKmeans[[1]]
    K <- imageKmeans[[2]]

    if(plot){
      x <- image.segmented/255
      cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
      uniqueCols <- unique(cols)
      x2 <- match(cols, uniqueCols)
      dim(x2) <- dim(x)[1:2]
      raster::image(t(apply(x2, 2, rev)), col=uniqueCols, yaxt='n', xaxt='n')
    }


    # Transform images and add to rasterList

    e=0

    rasterListInd <- list()

    for(i in 1:nrow(K$centers)){

      e=e+1

      rgb <- K$centers[i,]

      map <- apply(image.segmented, 1:2, function(x) all(x-rgb == 0))
      mapR <- raster::raster(map)
      raster::extent(mapR) <- extRasterOr

      mapDF <- raster::as.data.frame(mapR, xy = TRUE)

      mapDFs <- subset(mapDF, mapDF$layer == TRUE)

      invisible(capture.output(transMatrix <- Morpho::computeTransform(refShape, as.matrix(lanArray[,,n]), type = transformType)))

      invisible(capture.output(mapTransformed <- Morpho::applyTransform(as.matrix(mapDFs[1:2]), transMatrix)))

      r <- raster::raster(ncol = res, nrow = res)

      raster::extent(r) <- raster::extent(min(refShape[,1])-3*max(refShape[,1])*cropOffset[1]/100,
                                          max(refShape[,1])+3*max(refShape[,1])*cropOffset[2]/100,
                                          min(refShape[,2])-3*max(refShape[,2])*cropOffset[3]/100,
                                          max(refShape[,2])+3*max(refShape[,2])*cropOffset[4]/100)

      patternRaster <- raster::rasterize(mapTransformed, field = 1, r)

      rasterListInd[[e]] <- patternRaster


    rasterList[[names(landList)[n]]] <- rasterListInd
    }

    print(paste('sample', names(landList)[n], 'done and added to rasterList', sep=' '))
  }

  return(rasterList)

}

