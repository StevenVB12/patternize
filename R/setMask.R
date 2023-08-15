#' Interactive function to to draw an outline for masking.
#'
#' @param summedRaster Summed raster of extracted patterns.
#' @param IDlist List of sample IDs.
#' @param filename Name of file to which mask will be written.
#' @param ... additional arguments for plotHeat function.
#'
#' @return file
#'
#' @export
#' @import raster
#' @importFrom graphics locator
#' @importFrom utils write.table

setMask <- function(summedRaster,
                    IDlist,
                    filename,
                    ...){

  plotHeat(summedRaster, IDlist, ...)

  print("Choose points to mask patterns. Click outside image area to stop.")

  n = 1

  outline <- c()

  while(1){

    xy <- locator(n=1)

    n <- n + 1

    if(any(c(as.numeric(xy)[1] < raster::extent(summedRaster)[1],
             as.numeric(xy)[1] > raster::extent(summedRaster)[2],
             as.numeric(xy)[2] < raster::extent(summedRaster)[3],
             as.numeric(xy)[2] > raster::extent(summedRaster)[4]))){
      print("done")
      break
    }

    outline <- rbind(outline, as.numeric(xy))
    colnames(outline) <- c("x", "y")

    print(paste('x: ', as.character(xy)[1], 'y: ', as.character(xy)[2]))

    if(n > 1){
      lines(outline[c(n-1:n),], col = 'green', lwd = 2)
    }
  }

  write.table(outline, file = filename, quote = FALSE)

}
