#' Build list of landmarks or RasterStacks from images using filepath and file extension.
#'
#' @param IDlist List of sample IDs.
#' @param type 'landmark' or 'image' depending on what type of list to make.
#' @param prepath Prepath (default = NULL).
#' @param extension Extension (default = NULL).
#' @param format ImageJ (Fiji) or tps format (default = 'imageJ').
#' @param tpsFile Provide filename of tps file ff format is 'tps'.
#' @param skipLandmark Vector of rownumbers of landmarks to skip.
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
#' @importFrom geomorph readland.tps

makeList <- function(IDlist,
                     type,
                     prepath = NULL,
                     extension = NULL,
                     format = 'imageJ',
                     tpsFile = NULL,
                     skipLandmark = NULL){

  objectList <- list()

  if(!is.null(skipLandmark)){
    skipLandmark <- -1*skipLandmark
  }

  for(n in 1:length(IDlist)){

    if(format == 'imageJ'){
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

        if(!is.null(skipLandmark)){
          landmarks <- landmarks[skipLandmark,]
        }

        objectList[[IDlist[n]]] <- landmarks
      }
    }


    if(type == 'image'){

      if(is.null(prepath)){
        suppressWarnings(image <- raster::stack(paste(IDlist[n], extension, sep='')))
        crs(image) <- sp::CRS('+init=EPSG:4326')
      }
      else{
        suppressWarnings(image <- raster::stack(paste(prepath,'/',IDlist[n], extension, sep='')))
        crs(image) <- sp::CRS('+init=EPSG:4326')
      }

      objectList[[IDlist[n]]] <- image
    }
  }

  if(all(c(type == 'landmark', format == 'tps'))){

    objectListX <- readland.tps(tpsFile, specID = 'imageID', warnmsg = FALSE)
    objectList <- lapply(1:dim(objectListX)[3],function(i) objectListX[,,i])
    names(objectList) <- IDlist
  }

  return(objectList)
}


