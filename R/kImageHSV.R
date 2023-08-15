#' \code{\link[stats]{kmeans}} clustering of image imported as a RasterStack. This function is
#' used by \code{patLanK} and \code{patRegK}.
#'
#' @param image HSV image imported as a RasterStack for k-means clustering.
#' @param k Integer for number of k-means clusters (default = 3).
#' @param startCenter A matrix of cluster centres to start k-means clustering from (default = NULL).
#' @param maskToNA Replace the color value used for masking (i.e. 0 or 255) with NA.
#' @param kmeansOnAll Whether to perform the kmeans clusters on the combined set of pixels of all images
#'    first (default = FALSE).
#' @param ignoreHSVvalue Whether to ignore the HSV value (~darkness).
#'
#' @return List including the k-means clustered \code{RasterSatck} returned as an array and object
#'    of class "\code{kmeans}".
#'
#' @examples
#' image <- raster::stack(system.file("extdata", "BC0077.jpg", package = "patternize"))
#' out <- kImage(image, 6)
#'
#' @export
#' @import sf
#' @importFrom stats kmeans
#' @importFrom grDevices hsv rgb2hsv
#' @importFrom methods is

kImageHSV <- function(image,
                      k = 5,
                      startCenter = NULL,
                      maskToNA = NULL,
                      kmeansOnAll = FALSE,
                      ignoreHSVvalue = FALSE){

  if(kmeansOnAll == FALSE){

    if(is(image)[1] == "RasterStack"){
      image <- raster::as.array(image)
    }

    if(ignoreHSVvalue == TRUE){
      df = data.frame(
        hue        = matrix(image[,,1], ncol=1),
        saturation = matrix(image[,,2], ncol=1)
      )
    }
    else{
      df = data.frame(
        hue        = matrix(image[,,1], ncol=1),
        saturation = matrix(image[,,2], ncol=1),
        value      = matrix(image[,,3], ncol=1)
      )
    }

    if(is.null(startCenter)){
      K = kmeans(na.omit(df),k, nstart = 3)
    }
    else{
      K = kmeans(na.omit(df),startCenter)
    }
    df$label <- NA
    suppressWarnings(df$label[which(!is.na(df$hue))] <- K$cluster)

    # Replace color of each pixel with mean RGB value of cluster

    # get the coloring

    if(ignoreHSVvalue == TRUE){
      colors = data.frame(label = 1:nrow(K$centers),
                          H = K$centers[,"hue"],
                          S = K$centers[,"saturation"])
    }
    else{
      colors = data.frame(label = 1:nrow(K$centers),
                          H = K$centers[,"hue"],
                          S = K$centers[,"saturation"],
                          V = K$centers[,"value"])
    }

    # merge color codes on df

    df$order = 1:nrow(df)
    df = merge(df, colors, all = TRUE)
    df = df[order(df$order),]
    df$order = NULL

    # Reshape data frame back into an image

    if(ignoreHSVvalue == TRUE){
      H = matrix(df$H, nrow=dim(image)[1])
      S = matrix(df$S, nrow=dim(image)[1])
      V = matrix(1, nrow=nrow(image), ncol=ncol(image))

      image.segmented = array(dim=dim(image))
      image.segmented[,,1] = H
      image.segmented[,,2] = S
      image.segmented[,,3] = V
    }
    else{
      H = matrix(df$H, nrow=dim(image)[1])
      S = matrix(df$S, nrow=dim(image)[1])
      V = matrix(df$V, nrow=dim(image)[1])

      image.segmented = array(dim=dim(image))
      image.segmented[,,1] = H
      image.segmented[,,2] = S
      image.segmented[,,3] = V
    }

    out <- list(image.segmented, K)
  }


  if(kmeansOnAll == TRUE){

    for(n in 1:length(image)){

      imageX <- image[[n]]

      # convert RGB values to HSV values
      imageX <- raster::overlay(imageX, fun = rgb2hsv)

      if(!is.null(maskToNA)){
        imageX[imageX == maskToNA] <- NA
      }

      imageX <- raster::as.array(imageX)


      if(ignoreHSVvalue == TRUE){
        if(n==1){
          dfTot = data.frame(
            hue = matrix(imageX[,,1], ncol=1),
            saturation = matrix(imageX[,,2], ncol=1)
          )
          dfNrowTot <- c(nrow(dfTot))
        }
        else{
          df = data.frame(
            hue = matrix(imageX[,,1], ncol=1),
            saturation = matrix(imageX[,,2], ncol=1)
          )
          dfNrow <- c(nrow(df))
          dfNrowTot <- c(dfNrowTot, dfNrow)

          dfTot <- rbind(dfTot, df)
        }
      }
      else{
        if(n==1){
          dfTot = data.frame(
            hue = matrix(imageX[,,1], ncol=1),
            saturation = matrix(imageX[,,2], ncol=1),
            value = matrix(imageX[,,3], ncol=1)
          )
          dfNrowTot <- c(nrow(dfTot))
        }
        else{
          df = data.frame(
            hue = matrix(imageX[,,1], ncol=1),
            saturation = matrix(imageX[,,2], ncol=1),
            value = matrix(imageX[,,3], ncol=1)
          )
          dfNrow <- c(nrow(df))
          dfNrowTot <- c(dfNrowTot, dfNrow)

          dfTot <- rbind(dfTot, df)
        }
      }
    }

    if(is.null(startCenter)){
      K = kmeans(na.omit(dfTot),k, nstart = 3)
    }
    else{
      K = kmeans(na.omit(dfTot),startCenter)
    }

    dfTot$label <- NA
    suppressWarnings(dfTot$label[which(!is.na(dfTot$hue))] <- K$cluster)

    # get the coloring

    if(ignoreHSVvalue == TRUE){
      colors = data.frame(label = 1:nrow(K$centers),
                          H = K$centers[,"hue"],
                          S = K$centers[,"saturation"])
    }
    else{
      colors = data.frame(label = 1:nrow(K$centers),
                          H = K$centers[,"hue"],
                          S = K$centers[,"saturation"],
                          V = K$centers[,"value"])
    }

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

      if(ignoreHSVvalue == TRUE){
        H = matrix(df$H, nrow=dim(imageX)[1])
        S = matrix(df$S, nrow=dim(imageX)[1])
        V = matrix(1, nrow=nrow(imageX), ncol=ncol(imageX))

        image.segmented = array(dim=dim(imageX))
        image.segmented[,,1] = H
        image.segmented[,,2] = S
        image.segmented[,,3] = V
      }
      else{
        H = matrix(df$H, nrow=dim(imageX)[1])
        S = matrix(df$S, nrow=dim(imageX)[1])
        V = matrix(df$V, nrow=dim(imageX)[1])

        image.segmented = array(dim=dim(imageX))
        image.segmented[,,1] = H
        image.segmented[,,2] = S
        image.segmented[,,3] = V
      }

      image.segmented.list[[names(image)[n]]] <- image.segmented

      s <- s + dfNrowTot[n]
    }
    out <- list(image.segmented.list, K)
  }

  return(out)
}
