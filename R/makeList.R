#' Build list of landmarks or RasterStacks from images using filepath and file extension.
#'
#' @param IDlist List of sample IDs.
#' @param type 'landmark' or 'image' depending on what type of list to make.
#' @param prepath Prepath (default = NULL).
#' @param extension Extension (default = NULL).
#'
#' @return Landmark or RasterStack list.
#'
#' @examples
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#'
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#'
#' landmarkList <- makeList(IDlist, 'landmark', prepath, extension)
#'
#' extension <- '.jpg'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#'
#' @export
#' @importFrom utils read.table


makeList <- function(IDlist,
                     type,
                     prepath = NULL,
                     extension = NULL){

  objectList <- list()

  for(n in 1:length(IDlist)){

    print(paste('sample', n,  IDlist[n], 'added to list', sep=' '))

    if(type == 'landmark'){

      if(is.null(prepath)){
        landmarks <- read.table(paste(IDlist[n], extension, sep=''), header = FALSE,
                                stringsAsFactors = FALSE, colClasses = c('numeric', 'numeric'))
      }
      else{
        landmarks <- read.table(paste(prepath,'/',IDlist[n], extension, sep=''), header = FALSE,
                                stringsAsFactors = FALSE, colClasses = c('numeric', 'numeric'))
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


