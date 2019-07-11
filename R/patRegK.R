#' Aligns images using \code{\link[RNiftyReg]{niftyreg}} utilities for automated image registration
#' and extracts colors using k-means clustering.
#'
#' @param sampleList List of RasterStack objects.
#' @param target Image imported as RasterStack used as target for registration.
#' @param k Integere for defining number of k-means clusters (default = 3).
#' @param fixedStartCenter Specify a dataframe with start centers for k-means clustering.
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}} (default = NULL).
#' @param useBlockPercentage Block percentage as used in \code{\link[RNiftyReg]{niftyreg}}
#'    (default = 75).
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the
#'    original image.
#' @param removebgR Integer indicating the range RGB treshold to remove from image (e.g. 100 removes
#'    pixels with average RGB > 100; default = NULL) for registration analysis. This works only to
#'    remove a white background.
#' @param removebgK Integer indicating the range RGB treshold to remove from image (e.g. 100 removes
#'    pixels with average RGB > 100; default = NULL) for k-means analysis. This works only to remove
#'    a white background.
#' @param maskOutline When outline is specified, everything outside of the outline will be masked for
#'    the color extraction (default = NULL).
#' @param maskColor Color the masked area gets. Set to 0 for black (default) or 255 for white.
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
#'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#'
#' target <- imageList[[1]]
#'
#' \dontrun{
#' rasterList_regK <- patRegK(imageList[3], target, k = 5,
#' crop = c(100,400,40,250), removebgR = 100, plot = TRUE)
#' }
#'
#' @export

patRegK <- function(sampleList,
                    target,
                    k = 3,
                    fixedStartCenter = NULL,
                    resampleFactor = NULL,
                    useBlockPercentage = 75,
                    crop = c(0,0,0,0),
                    removebgR = NULL,
                    removebgK = NULL,
                    maskOutline = NULL,
                    maskColor = 0,
                    plot = FALSE,
                    focal =  FALSE,
                    sigma = 3){

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

    transformedMap <- RNiftyReg::applyTransform(RNiftyReg::forward(result), raster::as.array(sourceRasterK), interpolation=0)

    r1 <- raster::raster(transformedMap[1:nrow(transformedMap),ncol(transformedMap):1,1])
    r2 <- raster::raster(transformedMap[1:nrow(transformedMap),ncol(transformedMap):1,2])
    r3 <- raster::raster(transformedMap[1:nrow(transformedMap),ncol(transformedMap):1,3])


    transRaster <-raster::stack(r1,r2,r3)
    transRaster <- raster::flip(transRaster,'x')

    raster::extent(transRaster) <- raster::extent(sourceRasterK)

    if(!is.null(maskOutline)){
      transRaster <- maskOutline(transRaster, maskOutline, refShape = 'target', flipOutline = 'y', crop = crop,
                                 maskColor = maskColor, imageList = sampleList)
    }

    # k-means clustering of image

    if(n==1 & is.null(fixedStartCenter)){
      startCenter = NULL
    }

    if(!is.null(fixedStartCenter)){
      startCenter <- fixedStartCenter
      print(paste('Fixed centers:', startCenter, sep = ' '))
    }

    # else{
    #   startCenter <- K$centers
    # }

    if(is.null(removebgK)){
      imageKmeans <- tryCatch(kImage(raster::as.array(transRaster), k, startCenter),
                              error = function(err) {
                                print(paste('sample', names(sampleList)[n], 'k-clustering failed and skipped', sep = ' '))
                                return(NULL)
                                })
      if(is.null(imageKmeans)){next}
    }

    if(is.numeric(removebgK)){

      transRasterMean <- apply(raster::as.array(transRaster),1:2,mean)
      toMask <- apply(raster::as.array(transRasterMean), 1:2, function(x) ifelse(x > removebgK, NA, x))

      toMaskR <- raster::raster(as.matrix(toMask))
      raster::extent(toMaskR) <- raster::extent(transRaster)

      transRaster<-raster::mask(transRaster, toMaskR)
      transRaster[is.na(transRaster)] <- NA

      imageKmeans <- tryCatch(kImage(raster::as.array(transRaster), k, startCenter),
                              error = function(err) {
                                print(paste('sample', names(sampleList)[n], 'k-clustering failed and skipped', sep = ' '))
                                return(NULL)
                              })
      if(is.null(imageKmeans)){next}
    }

    image.segmented <- imageKmeans[[1]]
    K <- imageKmeans[[2]]

    if(n==1 & is.null(fixedStartCenter)){
      startCenter <- K$centers
      print('start centers of first image:')
      print(startCenter)
    }

    if(plot){
      image.segmented[is.na(image.segmented)] <- 0
      x <- image.segmented/255
      cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
      uniqueCols <- unique(cols)
      x2 <- match(cols, uniqueCols)
      dim(x2) <- dim(x)[1:2]
      raster::image(t(apply(x2,2,rev)), col=uniqueCols,yaxt='n', xaxt='n')
    }

    # print(names(sampleList)[n])

    # Transform images and add to rasterList

    e=0

    rasterListInd <- list()

    for(i in 1:nrow(K$centers)){

      e=e+1

      rgb <- K$centers[i,]

      map <- apply(image.segmented, 1:2, function(x) all(x-rgb == 0))

      map[map == 0] <- NA

      r <- raster::raster(map)
      raster::extent(r) <- extRaster

      rasterListInd[[e]] <- r

      if(!identical(raster::extent(rasterListInd[[e]]), raster::extent(target))){
        raster::extent(rasterListInd[[e]]) <- raster::extent(target)
      }

    }

    rasterList[[names(sampleList)[n]]] <- rasterListInd

    print(paste('sample', names(sampleList)[n], 'done and added to rasterList', sep=' '))
  }

  return(rasterList)
}
