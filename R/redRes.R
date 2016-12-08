#' Reduce the resolution of an iamge by downsampling.
#'
#' @param image Image for downsampling.
#' @param resampleFactor Integer for downsampling.
#'
#' @return Downsampled \code{image}.
#'
#' @examples
#' image <- raster::stack(system.file("extdata", "BC0077.jpg", package = "patternize"))
#' image.reduced <- redRes(image, 5)
#'
#' @export
#' @import rgdal

redRes <- function(image, resampleFactor){

  inCols <- ncol(image)
  inRows <- nrow(image)

  resampledRaster <- raster::raster(ncol=(inCols/resampleFactor), nrow=(inRows/resampleFactor))

  raster::extent(resampledRaster) <- raster::extent(image)

  resampled <- raster::resample(image,resampledRaster)

  return(resampled)

}
