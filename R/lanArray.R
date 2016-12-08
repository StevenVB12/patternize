#' Build landmark array for Morpho.
#'
#' @param sampleList List of landmark matrices.
#' @param adjustLandmarks Adjust the coordinates of the lanmarks (default = FALSE).
#'
#' @return  X x Y x n array, where X and Y define the coordinates of the landmark points and n is the sample size.
#'
#' @examples
#' data(landmarks)
#' landmarkArray <- lanArray(landmarks)
#'
#' @export


lanArray <- function(sampleList, adjustLandmarks=F){

  # Make datastructure for Generailzed Procrustis Analysis
  for(n in 1:length(sampleList)){
    print(paste('sample', n,  names(sampleList)[n], 'added to array', sep=' '))
    # Read in landmark files and build array for Morpho
    landmarks <- as.data.frame(sampleList[[n]])
#
#     if(adjustLandmarks){
#       library(raster)
#
#       pictureFile <- stack(paste(prePathPictures,"/",sampleList[sampleID],postPathPictures, sep=""))
#       landmarks <- adjustCoordinates(landmarks, pictureFile)
#     }

    if(n == 1){
      landmarksArray <- array(dim=c(1,2,1), data=landmarks)
    }
    else{
      landmarksArray1 <- array(dim=c(1,2,1), data=landmarks)
      landmarksArray  <- abind::abind(landmarksArray,landmarksArray1,along=3)
    }
  }
  return(landmarksArray)
}
