#' Image registration with k-means clustering.
#'
#' @param sampleList List of RasterStack objects.
#' @param target Image used as target for registration.
#' @param k Integere for defining number of k-means clusters (default = 3).
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param useBlockPercentage Block percentage as used in NiftyReg.
#' @param removebg Whether to remove white background rasterList (default = FALSE).
#' @param plot Whether to plot k-means clustered image while processing (default = FALSE).
#'
#' @return List of summed raster for each k-means cluster objects.
#'
#' @examples
#' data(pictures)
#' target <- pictures[[1]]
#' rasterList <- patRegK(pictures, target, k = 6, removebg = TRUE, plot = TRUE)
#'
#' @export

patRegK <- function(sampleList, target, k = 3, resampleFactor = 1, useBlockPercentage = 75, removebg = FALSE, plot = FALSE){

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
    source <- raster::as.array(source)
    sourceR <- apply(source,1:2,mean)

    if(removebg){
      sourceR <- apply(sourceR, 1:2, function(x) ifelse(x>100,0, x))
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
    if(n==1){

      for(i in 1:nrow(K$centers)){

        rgb <- K$centers[i,]

        map <- apply(image.segmented, 1:2, function(x) all(x-rgb == 0))

        transformedMap <- RNiftyReg::applyTransform(RNiftyReg::forward(result), map, interpolation=0)
        transformedMapMatrix <- transformedMap[1:nrow(transformedMap),ncol(transformedMap):1]

        r <- raster::raster(transformedMapMatrix)
        raster::extent(r)<-extRaster
        rasterList <- c(rasterList, r)
      }
    }

    else{

      e=0

      for(i in 1:nrow(K$centers)){

        e=e+1

        rgb <- K$centers[i,]

        map <- apply(image.segmented, 1:2, function(x) all(x-rgb == 0))

        transformedMap <- RNiftyReg::applyTransform(RNiftyReg::forward(result), map, interpolation=0)
        transtransformedMapMatrix <- transformedMap[1:nrow(transformedMap),ncol(transformedMap):1]
        # transtransformedMapMatrix[transtransformedMapMatrix == 0] <- NA

        r <- raster::raster(transtransformedMapMatrix)
        raster::extent(r)<-extRaster
        rasterList[[e]] <- rasterList[[e]] + r
      }
    }
  }

  #   # Sum raster list
  #   rasterList$fun <- sum
  #   rasterList$na.rm <- TRUE
  #   rrr <- do.call(mosaic,rasterList)

  return(rasterList)
}
