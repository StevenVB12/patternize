#' \code{\link[stats]{kmeans}} clustering of image imported as a RasterStack. This function is
#' used by \code{patLanK} and \code{patRegK}.
#'
#' @param image Image imported as a RasterStack for k-means clustering.
#' @param k Integer for number of k-means clusters (default = 3).
#' @param startCenter A matrix of cluster centres to start k-means clustering from (default = NULL).
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

kImage <- function(image, k = 5, startCenter = NULL){

  if(class(image) == "RasterStack"){
    image <- raster::as.array(image)
  }

  df = data.frame(
    red = matrix(image[,,1], ncol=1),
    green = matrix(image[,,2], ncol=1),
    blue = matrix(image[,,3], ncol=1)
  )

  if(is.null(startCenter)){
    K = kmeans(df,k, nstart = 3)
  }
  else{
    K = kmeans(df,startCenter)
  }
  df$label = K$cluster

  # Replace color of each pixel with mean RGB value of cluster

  # get the coloring

  colors = data.frame(label = 1:nrow(K$centers),
                      R = K$centers[,"red"],
                      G = K$centers[,"green"],
                      B = K$centers[,"blue"])

  # merge color codes on df

  df$order = 1:nrow(df)
  df = merge(df, colors)
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

  out <- list(image.segmented, K)
  return(out)
}
