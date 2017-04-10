#' Aligns images usings transformations obtained from fixed landmarks and extracts colors
#' using a predefined RGB values and cutoff value.
#'
#' @param sampleList List of RasterStack objects.
#' @param landList Landmark list as returned by \code{\link[patternize]{makeList}}.
#' @param RGB RGB values for color pattern extraction specified as vector.
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param colOffset Color offset for color pattern extraction (default = 0.10).
#' @param crop Whether to use the landmarks range to crop the image. This can significantly
#'    speed up the analysis (default = FALSE).
#' @param cropOffset Vector c(xmin, xmax, ymin, ymax) that specifies the number of pixels you
#'    want the cropping to be offset from the landmarks (in case the landmarks do not surround
#'    the entire color pattern).
#' @param res Resolution for color pattern raster (default = 300). This should be reduced if
#'    the number of pixels in the image is lower than th raster.
#' @param transformRef ID of reference sample for shape to which color patterns will be transformed
#'    to. Can be 'meanshape' for transforming to mean shape of Procrustes analysis.
#' @param transformType Transformation type as used by \code{\link[Morpho]{computeTransform}}
#'    (default ='tps').
#' @param adjustCoords Adjust landmark coordinates in case they are reversed compared to pixel
#'    coordinates (default = FALSE).
#' @param plot Whether to plot transformed color patterns while processing (default = FALSE).
#' @param focal Whether to perform Gaussian blurring (default = FALSE).
#' @param sigma Size of sigma for Gaussian blurring (default = 3).
#' @param iterations Number of iterations for recalculating average color.
#'
#' @return  List of raster objects.
#'
#' @examples
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#'
#' landmarkList <- makeList(IDlist, 'landmark', prepath, extension)
#'
#' extension <- '.jpg'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#'
#' RGB <- c(114,17,0)
#' rasterList_lanRGB <- patLanRGB(imageList, landmarkList, RGB,
#' colOffset = 0.15, crop = TRUE, res = 100, adjustCoords = TRUE, plot = TRUE)
#'
#' @export
#' @import raster Morpho
#' @importFrom utils capture.output
#' @importFrom stats na.omit


patLanRGB <- function(sampleList,
                      landList,
                      RGB,
                      resampleFactor = NULL,
                      colOffset = 0.10,
                      crop = FALSE,
                      cropOffset = NULL,
                      res = 300,
                      transformRef = 'meanshape',
                      transformType='tps',
                      adjustCoords = FALSE,
                      plot = FALSE,
                      focal =  FALSE,
                      sigma = 3,
                      iterations = 0){

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

  else{
    if(transformRef == 'meanshape'){

      invisible(capture.output(transformed <- Morpho::procSym(lanArray)))
      refShape <- transformed$mshape

    }

    else{

      if(transformRef %in% names(landList)){

        e <- which(names(landList) == transformRef)
        refShape <- lanArray[,,e]
      }

      else{
        stop("specified ID for reference shape does not exist")
      }
    }
  }

  for(n in 1:length(sampleList)){

    image <- sampleList[[n]]
    extRaster <- raster::extent(image)

    if(crop){

      landm <- lanArray[,,n]
      extRaster <- raster::extent(min(landm[,1]),
                                  max(landm[,1]),
                                  min(landm[,2]),
                                  max(landm[,2]))

      if(!is.null(cropOffset)){

        extRaster <- raster::extent(min(landm[,1])-cropOffset[1],
                                    max(landm[,1])+cropOffset[2],
                                    min(landm[,2])-cropOffset[3],
                                    max(landm[,2])+cropOffset[4])

      }

      image <- raster::crop(image, extRaster)
    }

    if(!is.null(resampleFactor)){
      image <- redRes(image, resampleFactor)
    }

    if(focal){
      gf <- focalWeight(image, sigma, "Gauss")

      rrr1 <- raster::focal(image[[1]], gf)
      rrr2 <- raster::focal(image[[2]], gf)
      rrr3 <- raster::focal(image[[3]], gf)

      image <- raster::stack(rrr1, rrr2, rrr3)
    }

    map <- apply(raster::as.array(image), 1:2, function(x) all(abs(x-RGB) < colOffset*255))

    if(all(map == FALSE)){
      stop("The RGB range does not seem to overlap with any of the RGB values in the image")
    }

    x <- 1
    while(x <= iterations){
      x <- x + 1

      mapRaster <- raster::raster(as.matrix(map))
      raster::extent(mapRaster) <- extRaster
      mapRaster[mapRaster == 0] <- NA

      raster::extent(image) <- extRaster
      mapMASK<-raster::mask(image, mapRaster)

      RGB <- c(mean(na.omit(as.data.frame(mapMASK[[1]]))[,1]),
               mean(na.omit(as.data.frame(mapMASK[[2]]))[,1]),
               mean(na.omit(as.data.frame(mapMASK[[3]]))[,1]))

      map <- apply(raster::as.array(image), 1:2, function(x) all(abs(x-RGB) < colOffset*255))

    }

    mapR <- raster::raster(map)
    raster::extent(mapR) <- extRaster

    mapDF <- raster::as.data.frame(mapR, xy = TRUE)

    mapDFs <- subset(mapDF, mapDF$layer == TRUE)

    invisible(capture.output(transMatrix <- Morpho::computeTransform(refShape, as.matrix(lanArray[,,n]), type = transformType)))

    invisible(capture.output(mapTransformed <- Morpho::applyTransform(as.matrix(mapDFs[1:2]), transMatrix)))

    r <- raster::raster(ncol = res, nrow = res)

    raster::extent(r) <- extent(min(refShape[,1])*1.4,max(refShape[,1])*1.4,min(refShape[,2])*1.4,max(refShape[,2])*1.4)

    if(!is.null(cropOffset)){
      raster::extent(r) <- extent(min(refShape[,1]),max(refShape[,1]),min(refShape[,2]),max(refShape[,2]))
    }

    patternRaster <- raster::rasterize(mapTransformed, field = 1, r)
    if(plot){

      if(n == 1){
        plot(1, type="n", xlab='', ylab='', xaxt='n', yaxt='n', axes= FALSE, bty='n')
      }

      par(new = TRUE)
      plot(patternRaster, col=rgb(1,0,0,alpha=1/length(sampleList)), legend = FALSE, xaxt='n', yaxt='n', axes= FALSE, bty='n')
    }


    rasterList[[names(landList)[n]]] <- patternRaster

    print(paste('sample', names(landList)[n], 'done and added to rasterList', sep=' '))
  }

  return(rasterList)
}
