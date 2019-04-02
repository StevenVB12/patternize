#' \code{\link[stats]{kmeans}} clustering of image imported as a RasterStack. This function is
#' used by \code{patLanK} and \code{patRegK}.
#'
#' @param image HSV image imported as a RasterStack for k-means clustering.
#' @param k Integer for number of k-means clusters (default = 3).
#' @param startCenter A matrix of cluster centres to start k-means clustering from (default = NULL).
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
#' @import rgdal
#' @importFrom stats kmeans

kImageHSV <- function(image,
                      k = 5,
                      startCenter = NULL,
                      ignoreHSVvalue = FALSE){

  if(class(image) == "RasterStack"){
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
    K = kmeans(df,k, nstart = 3)
  }
  else{
    K = kmeans(df,startCenter)
  }
  df$label = K$cluster

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
  df = merge(df, colors)
  df = df[order(df$order),]
  df$order = NULL

  # Reshape data frame back into an image

  if(ignoreHSVvalue == TRUE){
    H = matrix(df$H, nrow=dim(image)[1])
    S = matrix(df$S, nrow=dim(image)[1])
    V = matrix(0.5, nrow=nrow(image), ncol=ncol(image))

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
  return(out)
}
