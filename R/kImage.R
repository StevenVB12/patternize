#' kmeans clustering of image.
#'
#' @param image Image for kmeans clustering.
#' @param k Integer for number of kmeans clusters (default = 3).
#' @param startCenter A matrix of cluster centres to start kmeans clustering from (default = NULL).
#'
#' @return List with kmeans clustered \code{image} and object of class "kmeans".
#'
#' @examples
#' image <- raster::stack(system.file("extdata", "BC0077.jpg", package = "patternize"))
#' out <- kImage(image, 6)
#'
#' @export

kImage <- function(image, k = 3, startCenter = NULL){

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
