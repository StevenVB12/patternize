#' Color pattern quantification using landmarks and RGB color extraction.
#'
#' @param imageList List of RasterStack objects.
#' @param landmarkList Landmark landmarkList.
#' @param RGB RGB values for color pattern extraction specified as vector.
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param colOffset Color offset for color pattern extraction (default = 0).
#' @param crop Whether to use the landmarks range to crop the image. This can significantly speed up the analysis (default = FALSE).
#' @param cropOffset Vector c(xmin, xmax, ymin, ymax) that specifies the number of pixels you want the cropping to be offset from the landmarks (in case the landmarks do not surround the entire color pattern).
#' @param res Resolution for color pattern raster (default = 300). This should be reduced if the number of pixels in the image is lower than th raster.
#' @param transformRef ID of reference sample for shape to which color patterns will be transformed to. Can be 'meanshape' for transforming to mean shape of Procrustes analysis.
#' @param transformType (default ='tps')
#' @param adjustCoords Adjust landmark coordinates.
#' @param plot Whether to plot transformed color patterns while processing (default = FALSE).
#'
#' @return  List of raster objects.
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
#' RGB <- c(114,17,0)
#' rasterList_lanRGB <- patLanRGB(imageList, landmarkList, RGB, resampleFactor = 3, colOffset = 0.15, crop = TRUE, res = 150, adjustCoords = TRUE, plot = TRUE)
#'
#' @export
#' @import raster


patLanRGB <- function(imageList, landmarkList, RGB, resampleFactor = 1, colOffset = 0, crop = FALSE, cropOffset = NULL, res = 300, transformRef = 'meanshape', transformType='tps', adjustCoords = FALSE, plot = FALSE){

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

    image <- redRes(image, resampleFactor)

    map <- apply(raster::as.array(image), 1:2, function(x) all(abs(x-RGB) < colOffset*255))


    mapR <- raster::raster(map)
    raster::extent(mapR) <- extRaster

    mapDF <- raster::as.data.frame(mapR, xy = TRUE)

    mapDFs <- subset(mapDF, layer == TRUE)

    invisible(capture.output(transMatrix <- Morpho::computeTransform(refShape, lanArray[,,n], type = transformType)))

    invisible(capture.output(mapTransformed <- Morpho::applyTransform(as.matrix(mapDFs[1:2]), transMatrix)))

    r <- raster::raster(ncol = res, nrow = res)

    raster::extent(r) <- extent(min(refShape[,1]),max(refShape[,1]),min(refShape[,2]),max(refShape[,2]))

    patternRaster <- raster::rasterize(mapTransformed, field = 1, r)
    if(plot){

      if(n == 1){
        plot(1, type="n", xlab='', ylab='', xaxt='n', yaxt='n', axes= FALSE, bty='n')
      }

      par(new=T)
      plot(patternRaster, col=rgb(1,0,0,alpha=1/length(imageList)), legend = FALSE, xaxt='n', yaxt='n', axes= FALSE, bty='n')
    }


    rasterList[[names(landmarkList)[n]]] <- patternRaster
  }


  #   # Sum raster list
  #   rasterList$fun <- sum
  #   rasterList$na.rm <- TRUE
  #   rrr <- do.call(mosaic,rasterList)

  return(rasterList)
}
