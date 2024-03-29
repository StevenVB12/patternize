% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kImageHSV.R
\name{kImageHSV}
\alias{kImageHSV}
\title{\code{\link[stats]{kmeans}} clustering of image imported as a RasterStack. This function is
used by \code{patLanK} and \code{patRegK}.}
\usage{
kImageHSV(
  image,
  k = 5,
  startCenter = NULL,
  maskToNA = NULL,
  kmeansOnAll = FALSE,
  ignoreHSVvalue = FALSE
)
}
\arguments{
\item{image}{HSV image imported as a RasterStack for k-means clustering.}

\item{k}{Integer for number of k-means clusters (default = 3).}

\item{startCenter}{A matrix of cluster centres to start k-means clustering from (default = NULL).}

\item{maskToNA}{Replace the color value used for masking (i.e. 0 or 255) with NA.}

\item{kmeansOnAll}{Whether to perform the kmeans clusters on the combined set of pixels of all images
first (default = FALSE).}

\item{ignoreHSVvalue}{Whether to ignore the HSV value (~darkness).}
}
\value{
List including the k-means clustered \code{RasterSatck} returned as an array and object
   of class "\code{kmeans}".
}
\description{
\code{\link[stats]{kmeans}} clustering of image imported as a RasterStack. This function is
used by \code{patLanK} and \code{patRegK}.
}
\examples{
image <- raster::stack(system.file("extdata", "BC0077.jpg", package = "patternize"))
out <- kImage(image, 6)

}
