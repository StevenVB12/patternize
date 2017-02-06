#' Plot heatmap from summed rasterList
#'
#' @param summedRaster Summed raster or summedRasterList.
#' @param IDlist List of sample IDs.
#' @param colpalette Vector of colors for color palette (default = c("white","lightblue","blue","green", "yellow","red"))
#' @param plotCartoon Whether to plot a cartoon. This cartoon should be drawn on one of the samples used in the analysis.
#' @param type This can be 'target' in case the reference shape is a single sample (for registration analysis) or 'mean' if the images were transformed to a mean shape (only for meanshape when using landmark transformation)
#' @param outline xy coordinates that define outline.
#' @param lines list of files with xy coordinates of line objects to be added to cartoon.
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the original image used in landmark or registration analysis.
#' @param flipRaster Whether to flip raster along xy axis (in case there is an inconsistency between raster and outline coordinates).
#' @param flipOutline Whether to flip plot along x, y or xy axis.
#' @param imageList List of image should be given if one wants to flip the outline.
#'
#' @examples
#' data(rasterList_lanRGB)
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' outline_BC0077 <- read.table(paste(system.file("extdata",  package = 'patternize'), '/BC0077_outline.txt', sep=''), h= F)
#' lines_BC0077 <- list.files(path=paste(system.file("extdata",  package = 'patternize')), pattern='vein', full.names = T)
#'
#' summedRaster_lanRGB <- sumRaster(rasterList_lanRGB, IDlist, type = 'RGB')
#' plotHeat(summedRaster_lanRGB, IDlist)
#'
#' summedRaster_regRGB <- sumRaster(rasterList_regRGB, IDlist, type = 'RGB')
#' data(imageList)
#' plotHeat(summedRaster_regRGB, IDlist, plotCartoon = TRUE, type = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipRaster = 'xy', imageList = imageList)
#' plotHeat(summedRaster_regRGB, IDlist, plotCartoon = TRUE, type = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipRaster = 'x', flipOutline = 'y', imageList = imageList)
#' plotHeat(summedRaster_regRGB, IDlist, plotCartoon = TRUE, type = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipRaster = 'y', flipOutline = 'x', imageList = imageList)
#' plotHeat(summedRaster_regRGB, IDlist, plotCartoon = TRUE, type = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipOutline = 'xy', imageList = imageList)
#'
#' data(rasterList_lanK)
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' summedRasterList <- sumRaster(rasterList_lanK, IDlist, type = 'k')
#' plotHeat(summedRasterList, IDlist)
#'
#' summedRasterList_regK <- sumRaster(rasterList_regK, IDlist, type = 'k')
#' plotHeat(summedRasterList_regK, IDlist, plotCartoon = TRUE, type = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipRaster = 'xy', imageList = imageList)
#'
#' plotHeat(summedRasterList_regK[[1]], IDlist, plotCartoon = TRUE, type = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipRaster = 'xy', imageList = imageList)
#'
#' @export
#' @import raster

plotHeat <- function(summedRaster, IDlist, colpalette = NULL, plotCartoon = FALSE, type = NULL, outline = NULL, lines = NULL, crop = NULL, flipRaster = NULL, flipOutline = NULL, imageList = NULL){

  if(is.null(colpalette)){

    colfunc <- colorRampPalette(c("white","lightblue","blue","green", "yellow","red"))

  }

  else{

    if(!is.vector(colpalette)){

      stop('Specified color palette is not a vector')

    }

    colfunc <- colorRampPalette(colpalette)

  }

  if(!is.null(lines)){

    lineList <- list()

    for(e in 1:length(lines)){
      lineList[[e]] <- read.table(lines[e], h= FALSE)
    }

  }
  if(!is.null(flipOutline) || !is.null(flipRaster)){

    imageEx <- raster::extent(imageList[[1]])
    outline[,2] <- outline[,2] - crop[3]

    for(e in 1:length(lineList)){

      lineList[[e]][[2]] <- lineList[[e]][[2]] - crop[3]
    }

  }

  if(!is.null(flipOutline) && flipOutline == 'y' || !is.null(flipOutline) && flipOutline == 'xy'){

    outline[,2] <- outline[,2] + crop[3]

    for(e in 1:length(lineList)){

      lineList[[e]][[2]] <- lineList[[e]][[2]] + crop[3]
    }

  }

  if(is.null(flipOutline) && !is.null(flipRaster) || !is.null(flipOutline) && flipOutline == 'x'){

    outline[,2] <- outline[,2] + ((crop[3] - imageEx[3]) - (imageEx[4] - crop[4])) + crop[3]

    for(e in 1:length(lineList)){

      lineList[[e]][[2]] <- lineList[[e]][[2]] + ((crop[3] - imageEx[3]) - (imageEx[4] - crop[4])) + crop[3]
    }

  }

  if(!is.null(flipOutline)){

    if(flipOutline == 'x'){

        outline[,1] <- imageEx[2] - outline[,1] + (crop[1] - imageEx[1]) - (imageEx[2] - crop[2])

        for(e in 1:length(lineList)){

          lineList[[e]][[1]] <- imageEx[2] - lineList[[e]][[1]] + (crop[1] - imageEx[1]) - (imageEx[2] - crop[2])
        }

    }

    if(flipOutline == 'y'){

        outline[,2] <- imageEx[4] - outline[,2]


        for(e in 1:length(lineList)){

          lineList[[e]][[2]] <- imageEx[4] - lineList[[e]][[2]]

        }
    }

    if(flipOutline == 'xy'){

        outline[,1] <- imageEx[2] - outline[,1] + (crop[1] - imageEx[1]) - (imageEx[2] - crop[2])
        outline[,2] <- imageEx[4] - outline[,2]

        for(e in 1:length(lineList)){

          lineList[[e]][[1]] <- imageEx[2] - lineList[[e]][[1]] + (crop[1] - imageEx[1]) - (imageEx[2] - crop[2])
          lineList[[e]][[2]] <- imageEx[4] - lineList[[e]][[2]]

        }
    }
  }

  if(!is.list(summedRaster)){

    rasterEx <- raster::extent(summedRaster)

    if(is.null(outline)){
      XLIM <- c(rasterEx[1],rasterEx[2])
      YLIM <- c(rasterEx[3],rasterEx[4])
    }
    else{
      XLIM <- c(min(outline[,1]),max(outline[,1]))
      YLIM <- c(min(outline[,2]),max(outline[,2]))
    }

    if(!is.null(flipRaster)){
      if(flipRaster == 'x'){
          summedRaster <- raster::flip(summedRaster,'x')
      }
      if(flipRaster == 'y'){
          summedRaster <- raster::flip(summedRaster,'y')
      }
      if(flipRaster == 'xy'){
          summedRaster <- raster::flip(summedRaster,'x')
          summedRaster <- raster::flip(summedRaster,'y')
      }
    }

    par(mfrow=c(1,1), mai=c(0.05,0.8,0.05,0.8), oma=c(1,1,1,1)+1)

    plot(summedRaster/length(IDlist), col=colfunc(20), xaxt='n', yaxt='n', box=F, axes=F, xlim = XLIM, ylim= YLIM, zlim=c(0,1))

    if(plotCartoon){

      if(is.null(type) || is.null(outline)){

        stop('Not all paramters are set to plot the cartoon.')

      }

      if(type == 'target'){

        polygon(outline, col=NA, border='gray', xlim = XLIM, ylim= YLIM)

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lines(lineList[[e]], col='gray')

          }
        }
      }
    }
  }

  else{

    rasterEx <- raster::extent(summedRaster[[1]])

    if(is.null(outline)){
      XLIM <- c(rasterEx[1],rasterEx[2])
      YLIM <- c(rasterEx[3],rasterEx[4])
    }
    else{
      XLIM <- c(min(outline[,1]),max(outline[,1]))
      YLIM <- c(min(outline[,2]),max(outline[,2]))
    }

    par(mfrow=c(2,round((length(summedRaster)+1)/2)), mai=c(0.05,0.8,0.05,0.8), oma=c(1,1,1,1)+1)

    for(k in 1:length(summedRaster)){

      if(!is.null(flipRaster)){
        if(flipRaster == 'x'){
            summedRaster[[k]] <- raster::flip(summedRaster[[k]],'x')
        }
        if(flipRaster == 'y'){
            summedRaster[[k]] <- raster::flip(summedRaster[[k]],'y')
        }
        if(flipRaster == 'xy'){
            summedRaster[[k]] <- raster::flip(summedRaster[[k]],'x')
            summedRaster[[k]] <- raster::flip(summedRaster[[k]],'y')
        }
      }

      plot(summedRaster[[k]]/length(IDlist), col=colfunc(20), xaxt='n', yaxt='n', box=F, axes=F, xlim = XLIM, ylim= YLIM, zlim=c(0,1))

      if(plotCartoon){

      if(is.null(type) || is.null(outline)){

        stop('Not all paramters are set to plot the cartoon.')

      }

      if(type == 'target'){

        polygon(outline, col=NA, border='gray', xlim = XLIM, ylim= YLIM)

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lines(lineList[[e]], col='gray')

            }
          }
        }
      }
    }
  }
}
