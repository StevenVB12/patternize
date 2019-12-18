#' Extract colors using GMM clustering (for pre-aligned images).
#'
#' @param sampleList List of RasterStack objects.
#' @param k Integere for defining number of clusters (default = 3).
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param maskOutline When outline is specified, everything outside of the outline will be masked for
#'    the color extraction (default = NULL).
#' @param plot Whether to plot transformed color patterns while processing (default = FALSE).
#' @param focal Whether to perform Gaussian blurring (default = FALSE).
#' @param sigma Size of sigma for Gaussian blurring (default = 3).
#' @param maskToNA Replace the color value used for masking (i.e. 0 or 255) with NA.
#' @param kmeansOnAll Whether to perform the kmeans clusters on the combined set of pixels of all images
#'    first (default = FALSE).
#'
#' @return  List of summed raster for each k-means cluster objects.
#'
#'
#' @export
#' @import raster
#' @importFrom utils capture.output

patGMM <- function(sampleList,
                 k = 3,
                 resampleFactor = NULL,
                 maskOutline = NULL,
                 plot = FALSE,
                 focal =  FALSE,
                 sigma = 3,
                 maskToNA = NULL,
                 kmeansOnAll = FALSE){

  rasterList <- list()


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

    if(!is.null(maskToNA)){
      image[image == maskToNA] <- NA

    }

    # k-means clustering of image

    if(kmeansOnAll == FALSE){

      # imageKmeans <- tryCatch(GMMImage(raster::as.array(image), k),
      #                         error = function(err) {
      #                           print(paste('sample', names(sampleList)[n], 'k-clustering failed and skipped', sep = ' '))
      #                           return(NULL)
      #                         })
      imageKmeans <- GMMImage(raster::as.array(image), k)
      if(is.null(imageKmeans)){next}

      image.segmented <- imageKmeans[[1]]
      gmm <- imageKmeans[[2]]



      if(plot){
        image.segmented[is.na(image.segmented)] <- 0
        x <- image.segmented/255
        cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
        uniqueCols <- unique(cols)
        x2 <- match(cols, uniqueCols)
        dim(x2) <- dim(x)[1:2]
        raster::image(t(apply(x2, 2, rev)), col=uniqueCols, yaxt='n', xaxt='n')
      }


      e=0

      rasterListInd <- list()

      for(i in 1:nrow(gmm$centroids)){

        e=e+1

        rgb <- gmm$centroids[i,]

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

    imageKmeans <- GMMImage(sampleList, k, maskToNA, kmeansOnAll)

    images.segmented <- imageKmeans[[1]]
    gmm <- imageKmeans[[2]]



    for(n in 1:length(images.segmented)){

      image.segmented <- images.segmented[[n]]

      if(plot){
        image.segmented[is.na(image.segmented)] <- 0
        x <- image.segmented/255
        cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
        uniqueCols <- unique(cols)
        x2 <- match(cols, uniqueCols)
        dim(x2) <- dim(x)[1:2]
        raster::image(t(apply(x2, 2, rev)), col=uniqueCols, yaxt='n', xaxt='n')
      }

      e=0

      rasterListInd <- list()

      for(i in 1:nrow(gmm$centroids)){

        e=e+1

        rgb <- gmm$centroids[i,]

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

