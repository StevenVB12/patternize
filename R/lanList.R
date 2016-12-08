#' Build landmark list using filepath and file extension.
#'
#' @param IDlist List of sample IDs.
#' @param prepath Prepath (default = NULL).
#' @param extension Extension (default = NULL).
#'
#' @return Landmark list.
#'
#' @examples
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#' landmarkList <- lanList(IDlist, prepath, extension)
#'
#' @export


lanList <- function(IDlist, prepath = NULL, extension = NULL){

  landmarkList <- list()

  for(n in 1:length(IDlist)){

    if(is.null(prepath)){
      landmarks <- read.table(paste(IDlist[n], extension, sep=''), h= F)
    }
    else{
      landmarks <- read.table(paste(prepath,'/',IDlist[n], extension, sep=''), h= F)
    }

    colnames(landmarks)<-c('X','Y')

    landmarkList[[IDlist[n]]] <- landmarks
  }

  return(landmarkList)
}

