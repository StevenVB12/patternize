#' Color pattern quantification using landmarks and k-means clustering for color extraction.
#'
#' @param imageList List of RasterStack objects.
#' @param lanArray Landmark array.
#' @param k Integere for defining number of k-means clusters (default = 3).
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param crop Whether to use the landmarks range to crop the image. This can significantly speed up the analysis (default = FALSE).
#' @param cropOffset Vector c(xmin, xmax, ymin, ymax) that specifies the number of pixels you want the cropping to be offset from the landmarks (in case the landmarks do not surround the entire color pattern).
#' @param res Resolution for color pattern raster (default = 300). This should be reduced if the number of pixels in the image is lower than th raster.
#' @param transformRef ID of reference sample for shape to which color patterns will be transformed to. Can be 'meanshape' for transforming to mean shape of Procrustes analysis.
#' @param transformType (default ='tps')
#' @param removebg Whether to remove white background rasterList (default = FALSE) for k-means analysis.
#' @param adjustCoords Adjust coordinates.
#' @param plot Whether to plot transformed color patterns while processing (default = FALSE).
#'
#' @return  List of summed raster for each k-means cluster objects.
#'
#' @examples
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#' landmarkList <- makeList(IDlist, 'landmark', prepath, extension)
#'
#' extension <- '.jpg'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#'
#' rasterList_lanK <- patLanK(imageList, landmarkList, k = 5, resampleFactor = 3, crop = TRUE, res = 150, removebg = TRUE, adjustCoords = TRUE, plot = TRUE)
#'
#' @export
#' @import raster
#'
patLanK <- function(imageList, landmarkList, k = 3, resampleFactor = 1, crop = FALSE, cropOffset = NULL, res = 300, transformRef = 'meanshape', transformType='tps', removebg = FALSE, adjustCoords = FALSE, plot = FALSE){

  rasterList <- list()

  if(length(imageList) != length(landmarkList)){
    stop("imageList is not of the same length as lanArray")
  }

  for(n in 1:length(imageList)){
    if(names(imageList)[n] != names(landmarkList)[n]){
      stop("samples are not in the same order in imageList and lanArray")
    }
  }

  lanArray <- lanArray(landmarkList, adjustCoords, imageList)

  if(transformRef == 'meanshape'){

    invisible(capture.output(transformed <- Morpho::procSym(lanArray)))
    refShape <- transformed$mshape

  }

  else{

    if(exists(landmarkList[[transformRef]])){

      e <- which(names(landmarks) == transformRef)
      refShape <- lanArray[e]
    }

    else{
      stop("specified ID for reference shape does not exist")
    }
  }

  for(n in 1:length(imageList)){

    image <- imageList[[n]]
    extRaster <- raster::extent(image)

    if(crop){

      landm <- lanArray[,,n]
      extRaster <- raster::extent(min(landm[,1]), max(landm[,1]), min(landm[,2]), max(landm[,2]))

      if(!is.null(cropOffset)){

        extRaster <- raster::extent(min(landm[,1])-cropOffset[1], max(landm[,1])+cropOffset[2], min(landm[,2])-cropOffset[3], max(landm[,2])+cropOffset[4])

      }

      image <- raster::crop(image, extRaster)
    }

    imageRaster <- redRes(image, resampleFactor)
    image <- raster::as.array(imageRaster)

    if(removebg){

      imagena <- apply(image, 1:2, function(x) ifelse(all(x>100),NA, x))
      imagenar <- raster::raster(as.matrix(imagena))
      extent(imagenar)<-extent(imageRaster)

      imageMASK<-raster::mask(imageRaster, imagenar)
      imageMASK[is.na(imageMASK)] <- 0
      image <- raster::as.array(imageMASK)
    }

    # k-means clustering of image

    if(n==1){
      startCenter = NULL
    }

    else{
      startCenter <- K$centers
    }

    imageKmeans <- kImage(image, k, startCenter)

    image.segmented <- imageKmeans[[1]]
    K <- imageKmeans[[2]]

    if(plot){
      x <- image.segmented/255
      cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
      uniqueCols <- unique(cols)
      x2 <- match(cols, uniqueCols)
      dim(x2) <- dim(x)[1:2]
      raster::image(t(apply(x2,2,rev)), col=uniqueCols,yaxt='n', xaxt='n')
    }


    # Transform images and add to rasterList

    e=0

    rasterListInd <- list()

    for(i in 1:nrow(K$centers)){

      e=e+1

      rgb <- K$centers[i,]

      map <- apply(image.segmented, 1:2, function(x) all(x-rgb == 0))
      mapR <- raster::raster(map)
      raster::extent(mapR) <- extRaster

      mapDF <- raster::as.data.frame(mapR, xy = TRUE)

      mapDFs <- subset(mapDF, layer == TRUE)

      invisible(capture.output(transMatrix <- Morpho::computeTransform(refShape, lanArray[,,n], type = transformType)))

      invisible(capture.output(mapTransformed <- Morpho::applyTransform(as.matrix(mapDFs[1:2]), transMatrix)))

      r <- raster::raster(ncol = res, nrow = res)

      raster::extent(r) <- extent(min(refShape[,1]),max(refShape[,1]),min(refShape[,2]),max(refShape[,2]))

      patternRaster <- raster::rasterize(mapTransformed, field = 1, r)

      # if(n==1){
      #   rasterList <- c(rasterList, patternRaster)
      # }
      #
      # else{
      rasterListInd[[e]] <- patternRaster
    # }

    rasterList[[names(landmarkList)[n]]] <- rasterListInd
    }
  }

  return(rasterList)

}

