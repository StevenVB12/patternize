#' Plot heatmap from summed rasterList
#'
#' @param summedRaster Summed raster or summedRasterList.
#' @param IDlist List of sample IDs.
#' @param colpalette Vector of colors for color palette (default = c("white","lightblue","blue","green", "yellow","red"))
#' @param plotCartoon Whether to plot a cartoon. This cartoon should be drawn on one of the samples used in the analysis.
#' @param refShape This can be 'target' in case the reference shape is a single sample (for registration analysis) or 'mean' if the images were transformed to a mean shape (only for meanshape when using landmark transformation)
#' @param outline xy coordinates that define outline.
#' @param lines list of files with xy coordinates of line objects to be added to cartoon.
#' @param landList Landmark landmarkList.
#' @param adjustCoords Adjust landmark coordinates.
#' @param cartoonID ID of the sample for which the cartoon was drawn.
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the original image used in landmark or registration analysis.
#' @param flipRaster Whether to flip raster along xy axis (in case there is an inconsistency between raster and outline coordinates).
#' @param flipOutline Whether to flip plot along x, y or xy axis.
#' @param imageList List of image should be given if one wants to flip the outline or adjust landmark coordinates.
#' @param cartoonOrder Whether to plot the cartoon outline 'above' or 'under' the pattern raster (default = 'above'). Set to 'under' for filled outlines.
#' @param lineOrder Whether to plot the cartoon lines 'above' or 'under' the pattern raster (default = 'above').
#' @param cartoonCol Outline and line color for cartoon (deafault = 'gray').
#' @param cartoonFill Fill color for outline of cartoon (default = NULL).
#' @param zlim z-axis limit (default = c(0,1))
#' @param legend.title Title of the raster legend (default = 'Proportion')
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param main Optional main title.
#'
#' @examples
#' data(rasterList_lanRGB)
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' outline_BC0077 <- read.table(paste(system.file("extdata",  package = 'patternize'), '/BC0077_outline.txt', sep=''), h= F)
#' lines_BC0077 <- list.files(path=paste(system.file("extdata",  package = 'patternize')), pattern='vein', full.names = T)
#'
#'
#' summedRaster_regRGB <- sumRaster(rasterList_regRGB, IDlist, type = 'RGB')
#' data(imageList)
#' plotHeat(summedRaster_regRGB, IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipRaster = 'xy', imageList = imageList, cartoonOrder = 'under', cartoonFill = 'black', main = 'registration_example')
#' #plotHeat(summedRaster_regRGB, IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipRaster = 'x', flipOutline = 'y', imageList = imageList, cartoonOrder = 'under', cartoonFill = 'black', main = 'registration_example')
#' #plotHeat(summedRaster_regRGB, IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipRaster = 'y', flipOutline = 'x', imageList = imageList, cartoonOrder = 'under', cartoonFill = 'black', main = 'registration_example')
#' #plotHeat(summedRaster_regRGB, IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipOutline = 'xy', imageList = imageList, cartoonOrder = 'under', cartoonFill = 'black', main = 'registration_example')
#'
#' #data(rasterList_lanK)
#' #IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' #summedRasterList <- sumRaster(rasterList_lanK, IDlist, type = 'k')
#' #plotHeat(summedRasterList, IDlist)
#'
#' #summedRasterList_regK <- sumRaster(rasterList_regK, IDlist, refShape = 'k')
#' #plotHeat(summedRasterList_regK, IDlist, plotCartoon = TRUE, refShape = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipRaster = 'xy', imageList = imageList, cartoonOrder = 'under', cartoonFill = 'black', main = 'kmeans_example')
#'
#' #plotHeat(summedRasterList_regK[[1]], IDlist, plotCartoon = TRUE, type = 'target', outline = outline_BC0077, lines = lines_BC0077, crop = c(1000,4000,400,2500), flipRaster = 'xy', imageList = imageList, cartoonOrder = 'under', cartoonFill = 'black', main = 'kmeans_example')
#'
#'
#' prepath <- system.file("extdata", package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#' landmarkList <- makeList(IDlist, 'landmark', prepath, extension)
#' summedRaster_lanRGB <- sumRaster(rasterList_lanRGB, IDlist, type = 'RGB')
#' plotHeat(summedRaster_lanRGB, IDlist, plotCartoon = TRUE, refShape = 'mean', outline = outline_BC0077, lines = lines_BC0077, landList = landmarkList, adjustCoords = TRUE, imageList = imageList, cartoonID = 'BC0077', cartoonOrder = 'under', cartoonFill= 'black', main = 'Landmark_example')
#'
#' summedRaster_lanK <- sumRaster(rasterList_lanK, IDlist, type = 'k')
#' plotHeat(summedRaster_lanK, IDlist, plotCartoon = TRUE, refShape = 'mean', outline = outline_BC0077, lines = lines_BC0077, landList = landmarkList, adjustCoords = TRUE, imageList = imageList, cartoonID = 'BC0077', cartoonOrder = 'under', cartoonFill= 'black', main = 'Landmark_example')
#' plotHeat(summedRaster_lanK[[2]], IDlist, plotCartoon = TRUE, refShape = 'mean', outline = outline_BC0077, lines = lines_BC0077, landList = landmarkList, adjustCoords = TRUE, imageList = imageList, cartoonID = 'BC0077', cartoonOrder = 'under', cartoonFill= 'black', main = 'Landmark_example')
#'
#' @export
#' @import raster

plotHeat <- function(summedRaster, IDlist, colpalette = NULL, plotCartoon = FALSE, refShape = NULL, outline = NULL, lines = NULL, landList = NULL,  adjustCoords = FALSE, cartoonID = NULL, crop = c(0,0,0,0), flipRaster = NULL, flipOutline = NULL, imageList = NULL, cartoonOrder = 'above', lineOrder = 'above', cartoonCol = 'gray', cartoonFill = NULL, zlim = c(0,1), legend.title = 'Proportion', xlab='', ylab='', main=''){

  if(!is.list(summedRaster)){

    rasterEx <- raster::extent(summedRaster)

  }

  else{

    rasterEx <- raster::extent(summedRaster[[1]])

  }

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

    if(refShape != 'mean'){

      outline[,2] <- outline[,2] - crop[3]

      if(!is.null(lines)){

        for(e in 1:length(lineList)){

          lineList[[e]][[2]] <- lineList[[e]][[2]] - crop[3]
        }
      }
    }
  }

  if(!identical(crop, c(0,0,0,0))){
    if(refShape != 'mean'){

      if(!is.null(flipOutline) && flipOutline == 'y' || !is.null(flipOutline) && flipOutline == 'xy'){

        outline[,2] <- outline[,2] + crop[3]

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lineList[[e]][[2]] <- lineList[[e]][[2]] + crop[3]
          }
        }

      }

      if(is.null(flipOutline) && !is.null(flipRaster) || !is.null(flipOutline) && flipOutline == 'x'){

        outline[,2] <- outline[,2] + ((crop[3] - imageEx[3]) - (imageEx[4] - crop[4])) + crop[3]

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lineList[[e]][[2]] <- lineList[[e]][[2]] + ((crop[3] - imageEx[3]) - (imageEx[4] - crop[4])) + crop[3]
          }
        }
      }

    }

    if(!is.null(flipOutline)){

      if(flipOutline == 'x'){

        outline[,1] <- imageEx[2] - outline[,1] + (crop[1] - imageEx[1]) - (imageEx[2] - crop[2])

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lineList[[e]][[1]] <- imageEx[2] - lineList[[e]][[1]] + (crop[1] - imageEx[1]) - (imageEx[2] - crop[2])
          }
        }

      }

      if(flipOutline == 'y'){

        outline[,2] <- imageEx[4] - outline[,2]

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lineList[[e]][[2]] <- imageEx[4] - lineList[[e]][[2]]

          }
        }
      }

      if(flipOutline == 'xy'){

        outline[,1] <- imageEx[2] - outline[,1] + (crop[1] - imageEx[1]) - (imageEx[2] - crop[2])
        outline[,2] <- imageEx[4] - outline[,2]

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lineList[[e]][[1]] <- imageEx[2] - lineList[[e]][[1]] + (crop[1] - imageEx[1]) - (imageEx[2] - crop[2])
            lineList[[e]][[2]] <- imageEx[4] - lineList[[e]][[2]]

          }
        }
      }
    }
  }

  if(plotCartoon && refShape == 'mean'){

    indx <- which(IDlist == cartoonID)
    invisible(capture.output(landArray <- lanArray(landList, adjustCoords, imageList)))

    if(adjustCoords){

      extPicture <- extent(imageList[[indx]])
      outline[,2] <- extPicture[4]-outline[,2]

      if(!is.null(lines)){

        for(e in 1:length(lineList)){

          extPicture <- extent(imageList[[indx]])
          lineList[[e]][,2] <- extPicture[4]-lineList[[e]][,2]

        }
      }
    }

    invisible(capture.output(transformed <- Morpho::procSym(landArray)))


    invisible(capture.output(cartoonLandTrans <- Morpho::computeTransform(transformed$mshape, as.matrix(landArray[,,indx]), type="tps")))
    outlineTrans <- Morpho::applyTransform(as.matrix(outline), cartoonLandTrans)


    if(!is.null(lines)){

      cartoonLinesTrans <- list()
      for(e in 1:length(lineList)){

        cartoonLinesTrans[[e]] <- Morpho::applyTransform(as.matrix(lineList[[e]]), cartoonLandTrans)
      }
    }

    if(!is.null(flipOutline)){

      if(flipOutline == 'x'){
        outlineTrans[,1] = rasterEx[2] - outlineTrans[,1] + rasterEx[1]

        if(!is.null(cartoonLinesTrans)){
          for(e in 1:length(cartoonLinesTrans)){
            cartoonLinesTrans[[e]][,1] <- rasterEx[2] - cartoonLinesTrans[[e]][,1] + rasterEx[1]
          }
        }
      }

      if(flipOutline == 'y'){
        outlineTrans[,2] = rasterEx[4] - outlineTrans[,2] + rasterEx[3]

        if(!is.null(cartoonLinesTrans)){
          for(e in 1:length(cartoonLinesTrans)){
            cartoonLinesTrans[[e]][,2] <- rasterEx[4] - cartoonLinesTrans[[e]][,2] + rasterEx[3]
          }
        }
      }
      if(flipOutline == 'xy'){
        outlineTrans[,1] = rasterEx[2] - outlineTrans[,1] + rasterEx[1]
        outlineTrans[,2] = rasterEx[4] - outlineTrans[,2] + rasterEx[3]

        if(!is.null(cartoonLinesTrans)){
          for(e in 1:length(cartoonLinesTrans)){
            cartoonLinesTrans[[e]][,1] <- rasterEx[2] - cartoonLinesTrans[[e]][,1] + rasterEx[1]
            cartoonLinesTrans[[e]][,2] <- rasterEx[4] - cartoonLinesTrans[[e]][,2] + rasterEx[3]
          }
        }

      }
    }

    XLIM <- c(min(outlineTrans[,1]),max(outlineTrans[,1]))
    YLIM <- c(min(outlineTrans[,2]),max(outlineTrans[,2]))
  }




  if(!is.list(summedRaster)){

    if(is.null(outline) && is.null(lines)){
      XLIM <- c(rasterEx[1],rasterEx[2])
      YLIM <- c(rasterEx[3],rasterEx[4])
    }
    else{
      if(refShape == 'target'){
        XLIM <- c(min(outline[,1]),max(outline[,1]))
        YLIM <- c(min(outline[,2]),max(outline[,2]))
      }
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

    par(mfrow=c(1,1), mai=c(0.05,0.8,0.15,0.8), oma=c(1,1,1,1)+1)

    if(is.null(refShape) || refShape == 'target'){
      plot(NULL, type="n", axes=F, xlim = XLIM, ylim= YLIM, main=main, xlab = '', ylab='')
    }

    if(plotCartoon){

      if(is.null(refShape) || is.null(outline)){

        stop('Not all paramters are set to plot the cartoon.')

      }
    }

    if(plotCartoon && cartoonOrder == 'under'){

      summedRaster[summedRaster == 0] <- NA

      if(refShape == 'target'){

        polygon(outline, col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM)

      }

      if(refShape == 'mean'){

        plot(NULL, type="n", axes=F, xlim = XLIM, ylim= YLIM, main=main, xlab = '', ylab='')

        polygon(outlineTrans,col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM)

      }
    }

    if(plotCartoon && lineOrder == 'under'){

      if(refShape == 'target'){

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lines(lineList[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM)

          }
        }
      }

      if(refShape =='mean'){

        if(!is.null(lines)){

          for(e in 1:length(cartoonLinesTrans)){

            lines(cartoonLinesTrans[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM)
          }
        }
      }
    }

    plot(summedRaster/length(IDlist), col=colfunc(21), xaxt='n', yaxt='n', box=F, axes=F, xlim = XLIM, ylim= YLIM, zlim=zlim, legend.args=list(text=legend.title, side=4, line=3), add= TRUE)
    mtext(side = 1, text = xlab, line = 0)
    mtext(side = 2, text = ylab, line = 0)

    if(plotCartoon && cartoonOrder == 'above'){

      if(refShape == 'target'){

        polygon(outline, col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM)

      }

      if(refShape == 'mean'){

        polygon(outlineTrans,col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM)

      }
    }

    if(plotCartoon && lineOrder == 'above'){

      if(refShape == 'target'){

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lines(lineList[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM)

          }
        }
      }
      if(refShape =='mean'){

        if(!is.null(lines)){

          for(e in 1:length(cartoonLinesTrans)){

            lines(cartoonLinesTrans[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM)
          }
        }
      }
    }
  }

  else{

    if(is.null(outline)){
      XLIM <- c(rasterEx[1],rasterEx[2])
      YLIM <- c(rasterEx[3],rasterEx[4])
    }
    else{
      if(refShape == 'target'){
        XLIM <- c(min(outline[,1]),max(outline[,1]))
        YLIM <- c(min(outline[,2]),max(outline[,2]))
      }
    }

    par(mfrow=c(2,trunc((length(summedRaster)+1)/2)), mai=c(0.05,0.8,0.15,0.8), oma=c(1,1,1,1)+1)

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

      if(is.null(refShape) || refShape == 'target'){
        plot(NULL, type="n", axes=F, xlab="", ylab="", xlim = XLIM, ylim= YLIM, main= main)
      }

      if(plotCartoon){

        if(is.null(refShape) || is.null(outline)){

          stop('Not all paramters are set to plot the cartoon.')

        }
      }

      if(plotCartoon && cartoonOrder == 'under'){

        summedRaster[[k]][summedRaster[[k]] == 0] <- NA

        if(refShape == 'target'){

          polygon(outline, col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM)

        }

        if(refShape == 'mean'){

          plot(NULL, type="n", axes=F, xlim = XLIM, ylim= YLIM, main=main, xlab = '', ylab='')

          polygon(outlineTrans,col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM)

        }
      }

      if(plotCartoon && lineOrder == 'under'){

        if(refShape == 'target'){

          if(!is.null(lines)){

            for(e in 1:length(lineList)){

              lines(lineList[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM)

            }
          }
        }

        if(refShape =='mean'){

          if(!is.null(lineList)){

            for(e in 1:length(cartoonLinesTrans)){

              lines(cartoonLinesTrans[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM)
            }
          }
        }
      }

      plot(summedRaster[[k]]/length(IDlist), col=colfunc(21), xaxt='n', yaxt='n', box=F, axes=F, xlim = XLIM, ylim= YLIM, zlim=zlim, legend.args=list(text=legend.title, side=4, line=3), add= TRUE)
      mtext(side = 1, text = xlab, line = 0)
      mtext(side = 2, text = ylab, line = 0)

      if(plotCartoon && cartoonOrder == 'above'){

        if(refShape == 'target'){

          polygon(outline, col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM)

        }

        if(refShape == 'mean'){

          polygon(outlineTrans,col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM)

        }
      }

      if(plotCartoon && lineOrder == 'above'){

        if(refShape == 'target'){

          if(!is.null(lines)){

            for(e in 1:length(lineList)){

              lines(lineList[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM)

            }
          }
        }

        if(refShape =='mean'){

          if(!is.null(lines)){

            for(e in 1:length(cartoonLinesTrans)){

              lines(cartoonLinesTrans[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM)
            }
          }
        }
      }
    }
  }
}



