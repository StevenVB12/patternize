#' Build list of landmarks or picture stacks using filepath and file extension.
#'
#' @param IDlist List of sample IDs.
#' @param type 'landmark' or 'image'
#' @param prepath Prepath (default = NULL).
#' @param extension Extension (default = NULL).
#'
#' @return Landmark or image list.
#'
#' @examples
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#' landmarkList <- makeList(IDlist, 'landmark', prepath, extension)
#'
#' extension <- '.jpg'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#'
#' @export


makeList <- function(IDlist, type, prepath = NULL, extension = NULL){

  objectList <- list()

  for(n in 1:length(IDlist)){

    if(type == 'landmark'){

      if(is.null(prepath)){
        landmarks <- read.table(paste(IDlist[n], extension, sep=''), h= F)
      }
      else{
        landmarks <- read.table(paste(prepath,'/',IDlist[n], extension, sep=''), h= F)
      }

      landmarks <- as.matrix(landmarks)
      colnames(landmarks) <- NULL

      objectList[[IDlist[n]]] <- landmarks
    }


    if(type == 'image'){

      if(is.null(prepath)){
        image <- raster::stack(paste(IDlist[n], extension, sep=''))
      }
      else{
        image <- raster::stack(paste(prepath,'/',IDlist[n], extension, sep=''))
      }

      objectList[[IDlist[n]]] <- image
    }
  }

  return(objectList)
}

