#' Extract colors using k-means clustering (for pre-aligned images).
#'
#' @param sampleList List of RasterStack objects.
#' @param k Integere for defining number of k-means clusters (default = 3).
#' @param fixedStartCenter Specify a dataframe with start centers for k-means clustering.
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param maskOutline When outline is specified, everything outside of the outline will be masked for
#'    the color extraction (default = NULL).
#' @param plot Whether to plot transformed color patterns while processing (default = FALSE).
#' @param focal Whether to perform Gaussian blurring (default = FALSE).
#' @param sigma Size of sigma for Gaussian blurring (default = 3).
#' @param maskToNA Replace the color value used for masking (i.e. 0 or 255) with NA.
#' @param kmeansOnAll Whether to perform the kmeans clusters on the combined set of pixels of all images
#'    first (default = FALSE).
#' @param ignoreHSVvalue Whether to ignore the HSV value (~darkness).
#'
#' @return  List of summed raster for each k-means cluster objects.
#'
#'
#' @export
#' @import raster
#' @importFrom utils capture.output
#' @importFrom grDevices hsv rgb2hsv

patK_HSV <- function(sampleList,
                 k = 3,
                 fixedStartCenter = NULL,
                 resampleFactor = NULL,
                 maskOutline = NULL,
                 plot = FALSE,
                 focal =  FALSE,
                 sigma = 3,
                 maskToNA = NULL,
                 kmeansOnAll = FALSE,
                 ignoreHSVvalue = FALSE){

  rasterList <- list()

  if(is.null(fixedStartCenter)){
    startCenter = NULL
  }

  if(!is.null(fixedStartCenter)){
    startCenter <- fixedStartCenter
    print('Fixed start centers:')
    print(startCenter)
  }

  for(n in 1:length(sampleList)){

    image <- sampleList[[n]]



    extRasterOr <- raster::extent(image)

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

    if(!is.null(maskOutline))(
      image <- maskOutline(image, maskOutline, refShape = 'target', flipOutline = 'y', imageList = sampleList)
    )



    # k-means clustering of image

    if(kmeansOnAll == FALSE){

      # convert RGB values to HSV values

      image <- raster::overlay(image, fun = rgb2hsv)

      if(!is.null(maskToNA)){
        image[image == maskToNA] <- NA

      }

      imageKmeans <- kImageHSV(raster::as.array(image), k, startCenter, ignoreHSVvalue = ignoreHSVvalue)

      imageKmeans <- tryCatch(kImageHSV(raster::as.array(image), k, startCenter, ignoreHSVvalue = ignoreHSVvalue),
                              error = function(err) {
                                print(paste('sample', names(sampleList)[n], 'k-clustering failed and skipped', sep = ' '))
                                return(NULL)
                              })
      # imageKmeans <- kImage(raster::as.array(image), k, startCenter)
      if(is.null(imageKmeans)){next}

      image.segmented <- imageKmeans[[1]]
      K <- imageKmeans[[2]]

      if(all(c(n==1, is.null(fixedStartCenter)))){
        startCenter <- K$centers
        print('start centers of first image:')
        print(startCenter)
      }

      if(plot){
        image.segmented[is.na(image.segmented)] <- 0
        x <- image.segmented
        cols <- hsv(x[,,1], x[,,2], x[,,3])
        uniqueCols <- unique(cols)
        x2 <- match(cols, uniqueCols)
        dim(x2) <- dim(x)[1:2]
        raster::image(t(apply(x2, 2, rev)), col=uniqueCols, yaxt='n', xaxt='n')
      }


      e=0

      rasterListInd <- list()

      for(i in 1:nrow(K$centers)){

        e=e+1

        if(ignoreHSVvalue == FALSE){
          rgb <- K$centers[i,]
        }
        else{
          rgb <- K$centers[i,]
          rgb <- c(rgb, 1)
        }

        map <- apply(image.segmented, 1:2, function(x) all(x-rgb == 0))
        mapR <- raster::raster(map)
        raster::extent(mapR) <- extRasterOr

        rasterListInd[[e]] <- mapR


        rasterList[[names(sampleList)[n]]] <- rasterListInd
      }

      print(paste('sample', names(sampleList)[n], 'done and added to rasterList', sep=' '))
    }
  }

  if(kmeansOnAll == TRUE){

    imageKmeans <- kImageHSV(sampleList, k, startCenter, maskToNA, kmeansOnAll, ignoreHSVvalue)

    images.segmented <- imageKmeans[[1]]
    K <- imageKmeans[[2]]

    # if(!is.null(fixedStartCenter)){
    #   print('start centers of all images:')
    #   print(startCenter)
    # }
    startCenter <- K$centers
    print('final k-means centers of all images:')
    print(startCenter)


    for(n in 1:length(images.segmented)){

      image.segmented <- images.segmented[[n]]

      if(plot){
        image.segmented[is.na(image.segmented)] <- 0
        x <- image.segmented
        cols <- hsv(x[,,1], x[,,2], x[,,3])
        uniqueCols <- unique(cols)
        x2 <- match(cols, uniqueCols)
        dim(x2) <- dim(x)[1:2]
        raster::image(t(apply(x2, 2, rev)), col=uniqueCols, yaxt='n', xaxt='n')
      }


      e=0

      rasterListInd <- list()

      for(i in 1:nrow(K$centers)){

        e=e+1

        if(ignoreHSVvalue == FALSE){
          rgb <- K$centers[i,]
        }
        else{
          rgb <- K$centers[i,]
          rgb <- c(rgb, 0.5)
        }

        map <- apply(image.segmented, 1:2, function(x) all(x-rgb == 0))
        mapR <- raster::raster(map)
        raster::extent(mapR) <- extRasterOr

        rasterListInd[[e]] <- mapR


        rasterList[[names(sampleList)[n]]] <- rasterListInd
      }
      print(paste('sample', names(sampleList)[n], 'done and added to rasterList', sep=' '))
    }

  }
  return(rasterList)

}

