#' Reduce the resolution of an image imported as a RasterStack by downsampling.
#'
#' @param image RasterStack for downsampling.
#' @param resampleFactor Integer for downsampling.
#'
#' @return Downsampled RasterStack
#'
#' @examples
#' image <- raster::stack(system.file("extdata", "BC0077.jpg", package = "patternize"))
#' image_reduced <- redRes(image, 5)
#'
#' @export
#' @import raster

redRes <- function(image,
                   resampleFactor){

  inCols <- ncol(image)
  inRows <- nrow(image)

  resampledRaster <- raster::raster(ncol=(inCols/resampleFactor), nrow=(inRows/resampleFactor))

  raster::extent(resampledRaster) <- raster::extent(image)

  resampled <- raster::resample(image, resampledRaster)

  return(resampled)

}
