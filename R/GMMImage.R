#' \code{\link[ClusterR]{GMM}} clustering of image imported as a RasterStack.
#'
#' @param image Image imported as a RasterStack for clustering.
#' @param k Integer for number of k clusters (default = 3).
#' @param maskToNA Replace the color value used for masking (i.e. 0 or 255) with NA.
#' @param kmeansOnAll Whether to perform the kmeans clusters on the combined set of pixels of all images
#'    first (default = FALSE).
#'
#' @return List including the clustered \code{RasterSatck} returned as an array and object
#'    of class "\code{GMM}".
#'
#' @export
#' @import sf
#' @importFrom stats kmeans
#' @importFrom ClusterR GMM predict_GMM
#' @importFrom methods is

GMMImage <- function(image,
                   k = 5,
                   maskToNA = NULL,
                   kmeansOnAll = FALSE){

  if(kmeansOnAll == FALSE){

    if(is(image)[1] == "RasterStack"){
      image <- raster::as.array(image)
    }

    df = data.frame(
      red = matrix(image[,,1], ncol=1),
      green = matrix(image[,,2], ncol=1),
      blue = matrix(image[,,3], ncol=1)
    )

    # if(is.null(startCenter)){
    #   K = kmeans(na.omit(df),k, nstart = 3)
    # }
    # else{
    #   K = kmeans(na.omit(df),startCenter)
    # }

    gmm = GMM(na.omit(df), k, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 0, em_iter = 10, verbose = F)

    pr = predict_GMM(na.omit(df), gmm$centroids, gmm$covariance_matrices, gmm$weights)

    df$label <- NA
    suppressWarnings(df$label[which(!is.na(df$red))] <- pr$cluster_labels)

    # df[is.na(df)] <- 0
    # df$label = K$cluster

    # Replace color of each pixel with mean RGB value of cluster

    # get the coloring

    colors = data.frame(label = 0:(nrow(gmm$centroids)-1),
                        R = gmm$centroids[,1],
                        G = gmm$centroids[,2],
                        B = gmm$centroids[,3])

    # merge color codes on df

    df$order = 1:nrow(df)
    df = merge(df, colors, all = TRUE)
    df = df[order(df$order),]
    df$order = NULL

    # Reshape data frame back into an image

    R = matrix(df$R, nrow=dim(image)[1])
    G = matrix(df$G, nrow=dim(image)[1])
    B = matrix(df$B, nrow=dim(image)[1])

    image.segmented = array(dim=dim(image))
    image.segmented[,,1] = R
    image.segmented[,,2] = G
    image.segmented[,,3] = B

    out <- list(image.segmented, gmm)
  }

  if(kmeansOnAll == TRUE){

    for(n in 1:length(image)){

      imageX <- image[[n]]
      imageX <- raster::as.array(imageX)

      if(!is.null(maskToNA)){
        imageX[imageX == maskToNA] <- NA
      }

      if(n==1){
        dfTot = data.frame(
          red = matrix(imageX[,,1], ncol=1),
          green = matrix(imageX[,,2], ncol=1),
          blue = matrix(imageX[,,3], ncol=1)
        )
        dfNrowTot <- c(nrow(dfTot))
      }
      else{
        df = data.frame(
          red = matrix(imageX[,,1], ncol=1),
          green = matrix(imageX[,,2], ncol=1),
          blue = matrix(imageX[,,3], ncol=1)
        )
        dfNrow <- c(nrow(df))
        dfNrowTot <- c(dfNrowTot, dfNrow)

        dfTot <- rbind(dfTot, df)

      }
    }

    # if(is.null(startCenter)){
    #   K = kmeans(na.omit(dfTot),k, nstart = 3)
    # }
    # else{
    #   K = kmeans(na.omit(dfTot),startCenter)
    # }

    gmm = GMM(na.omit(dfTot), k, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 0, em_iter = 10, verbose = F)

    pr = predict_GMM(na.omit(dfTot), gmm$centroids, gmm$covariance_matrices, gmm$weights)

    dfTot$label <- NA
    suppressWarnings(dfTot$label[which(!is.na(dfTot$red))] <- pr$cluster_labels)

    # get the coloring

    colors = data.frame(label = 0:(nrow(gmm$centroids)-1),
                        R = gmm$centroids[,1],
                        G = gmm$centroids[,2],
                        B = gmm$centroids[,3])

    # merge color codes on df

    dfTot$order = 1:nrow(dfTot)
    dfTot = merge(dfTot, colors, all = TRUE)
    dfTot = dfTot[order(dfTot$order),]
    dfTot$order = NULL



    # Reshape data frame back into an image
    s <- 1
    e <- 0
    image.segmented.list <- list()

    for(n in 1:length(image)){

      imageX <- image[[n]]
      imageX <- raster::as.array(imageX)

      e <- e + dfNrowTot[n]

      df <- dfTot[c(s:e),]

      R = matrix(df$R, nrow=dim(imageX)[1])
      G = matrix(df$G, nrow=dim(imageX)[1])
      B = matrix(df$B, nrow=dim(imageX)[1])

      image.segmented = array(dim=dim(imageX))
      image.segmented[,,1] = R
      image.segmented[,,2] = G
      image.segmented[,,3] = B

      image.segmented.list[[names(image)[n]]] <- image.segmented

      s <- s + dfNrowTot[n]
    }
    out <- list(image.segmented.list, gmm)
  }
  return(out)
}
