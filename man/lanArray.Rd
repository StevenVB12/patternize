% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lanArray.R
\name{lanArray}
\alias{lanArray}
\title{Build landmark array for \code{\link[Morpho]{Morpho}}.}
\usage{
lanArray(sampleList, adjustCoords = FALSE, imageList = NULL, imageIDs = NULL)
}
\arguments{
\item{sampleList}{List of landmark matrices as returned by \code{\link{makeList}}.}

\item{adjustCoords}{Adjust landmark coordinates in case they are reversed compared to pixel
coordinates (default = FALSE).}

\item{imageList}{List of RasterStacks as returned by \code{\link{makeList}} should be given
when \code{adjustCoords = TRUE}.}

\item{imageIDs}{A list of IDs to match landmarks to images if landmarkList and imageList don't
have the same length.}
}
\value{
X x Y x n array, where X and Y define the coordinates of the landmark points and n
   is the sample size.
}
\description{
Build landmark array for \code{\link[Morpho]{Morpho}}.
}
\examples{

\dontrun{
IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')

prepath <- system.file("extdata",  package = 'patternize')
extension <- '_landmarks_LFW.txt'

landmarkList <- makeList(IDlist, 'landmark', prepath, extension)

landmarkArray <- lanArray(landmarkList)
}

}
