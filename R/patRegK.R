#' Image registration with k-means clustering for color extraction.
#'
#' @param sampleList List of RasterStack objects.
#' @param target Image used as target for registration.
#' @param k Integere for defining number of k-means clusters (default = 3).
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param useBlockPercentage Block percentage as used in NiftyReg.
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the original image.
#' @param removebg Whether to remove white background rasterList (default = FALSE) for registration and k-means analysis.
#' @param plot Whether to plot k-means clustered image while processing (default = FALSE).
#' @param focal Whether to perform Gaussian blurring (default = FALSE).
#' @param sigma Size of sigma for Gaussian blurring (default = 3).
#'
#' @return List of rasters for each k-means cluster objects.
#'
#' @examples
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '.jpg'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#' target <- imageList[[1]]
#' rasterList_regK <- patRegK(imageList, target, k = 5, resampleFactor = 10, crop = c(1000,4000,400,2500), removebg = TRUE, plot = TRUE)
#'
#' @export

patRegK <- function(sampleList, target, k = 3, resampleFactor = 1, useBlockPercentage = 75, crop = NULL, removebg = FALSE, plot = FALSE, focal =  FALSE, sigma = 3){

  imageList <- sampleList

  rasterList <- list()

  if(!is.null(crop)){

    targetExtRaster <- crop
    target <- raster::crop(target, targetExtRaster)

  }

  target <- redRes(target, resampleFactor)
  target <- raster::as.array(target)
  target <- apply(target,1:2,mean)

  if(removebg){
    target <- apply(target, 1:2, function(x) ifelse(x>100,0, x))
  }

  for(n in 1:length(sampleList)){

    sStack <- sampleList[[n]]
    extRaster <- raster::extent(sStack)

    if(!is.null(crop)){

      extRaster <- crop
      sStack <- crop(sStack, extRaster)

    }

    sourceRaster <- redRes(sStack, resampleFactor)

    if(focal){
      gf <- focalWeight(sourceRaster, sigma, "Gauss")

      rrr1 <- raster::focal(sourceRaster[[1]], gf)
      rrr2 <- raster::focal(sourceRaster[[2]], gf)
      rrr3 <- raster::focal(sourceRaster[[3]], gf)

      sourceRaster <- raster::stack(rrr1, rrr2, rrr3)
    }

    source <- raster::as.array(sourceRaster)
    sourceR <- apply(source,1:2,mean)

    if(removebg){

      sourceR <- apply(sourceR, 1:2, function(x) ifelse(x > 100, 0, x))
      sourceRna <- apply(source, 1:2, function(x) ifelse(all(x > 100), NA, x))
      sourceRnar <- raster::raster(as.matrix(sourceRna))
      extent(sourceRnar)<-extent(sourceRaster)

      sourceMASK<-raster::mask(sourceRaster, sourceRnar)
      sourceMASK[is.na(sourceMASK)] <- 0
      source <- raster::as.array(sourceMASK)
    }

    result <- RNiftyReg::niftyreg(sourceR, target, useBlockPercentage=useBlockPercentage)


    # k-means clustering of image

    if(n==1){
      startCenter = NULL
    }

    else{
      startCenter <- K$centers
    }

    imageKmeans <- kImage(source, k, startCenter)

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

      transformedMap <- RNiftyReg::applyTransform(RNiftyReg::forward(result), map, interpolation=0)
      transformedMapMatrix <- transformedMap[1:nrow(transformedMap),ncol(transformedMap):1]
      transformedMapMatrix[transformedMapMatrix == 0] <- NA

      r <- raster::raster(transformedMapMatrix)
      raster::extent(r) <- extRaster

      rasterListInd[[e]] <- r

    }
    rasterList[[names(imageList)[n]]] <- rasterListInd
  }

  return(rasterList)
}
