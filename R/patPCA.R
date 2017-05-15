#' This function transforms the individual color pattern rasters as obtained by the main
#' patternize functions to a dataframe of 0 and 1 values that can be used for Principal
#' Component Analysis (\code{\link[stats]{prcomp}}). This function also allows to plot the
#' analysis including a visualization of the shape changes along the axis. Pixel values
#' are predicted by multiplying the rotation matrix (eigenvectors) with a vector that has
#' the same length as the number of rows in the rotation matrix and in which all values are
#' set to zero except for the PC value for which we want to predict the pixel values.
#'
#' @param rList List of raster objects.
#' @param popList List of vectors including sampleIDs for eacht population.
#' @param colList List of colors for each population.
#' @param plot Whether to plot the PCA analysis.
#' @param plotType Plot 'points' or sample 'labels' (default = 'points')
#' @param plotChanges Wether to include plots of the changes along the PC axis (default = FALSE).
#' @param PCx PC axis to be presented for x-axis (default PC1).
#' @param PCy PC axis to be presented for y-axis (default PC2).
#' @param colpalette Vector of colors for color palette
#'    (default = c("white","lightblue","blue","green", "yellow","red"))
#' @param plotCartoon Whether to plot a cartoon. This cartoon should be drawn on one of the
#'    samples used in the analysis.
#' @param refShape This can be 'target' in case the reference shape is a single sample (for
#'    registration analysis) or 'mean' if the images were transformed to a mean shape (only for
#'    meanshape when using landmark transformation)
#' @param outline xy coordinates that define outline.
#' @param lines list of files with xy coordinates of line objects to be added to cartoon.
#' @param landList Landmark landmarkList.
#' @param adjustCoords Adjust landmark coordinates.
#' @param cartoonID ID of the sample for which the cartoon was drawn.
#' @param normalized Set this to true in case the summed rasters are already devided by the
#'    sample number.
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop
#'    the original image used in landmark or registration analysis.
#' @param flipRaster Whether to flip raster along xy axis (in case there is an inconsistency
#'    between raster and outline coordinates).
#' @param flipOutline Whether to flip plot along x, y or xy axis.
#' @param imageList List of image should be given if one wants to flip the outline or adjust
#'    landmark coordinates.
#' @param cartoonOrder Whether to plot the cartoon outline 'above' or 'under' the pattern raster
#'    (default = 'above'). Set to 'under' for filled outlines.
#' @param lineOrder Whether to plot the cartoon lines 'above' or 'under' the pattern raster
#'    (default = 'above').
#' @param cartoonCol Outline and line color for cartoon (deafault = 'gray').
#' @param cartoonFill Fill color for outline of cartoon (default = NULL).
#' @param plotLandmarks Whether to plot the landmarks from the target image or mean shape
#'    landmarks (default = FALSE).
#' @param landCol Color for ploting landmarks (default = 'black').
#' @param zlim z-axis limit (default = c(0,1))
#' @param legendTitle Title of the raster legend (default = 'Proportion')
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param main Optional main title.
#'
#' @return  List including a [1] dataframe of the binary raster values that can be used for
#'    principle component analysis, [2] a dataframe of sample IDs and specified population
#'    colors and [3] prcomp results.
#'
#' @seealso \code{\link[stats]{prcomp}}
#'
#' @examples
#' data(rasterList_lanRGB)
#'
#' pop1 <- c('BC0077','BC0071')
#' pop2 <- c('BC0050','BC0049','BC0004')
#' popList <- list(pop1, pop2)
#' colList <- c("red", "blue")
#'
#' pcaOut <- patPCA(rasterList_lanRGB, popList, colList, plot = TRUE)
#'
#' @export
#' @import raster
#' @importFrom stats prcomp


patPCA <- function(rList,
                   popList,
                   colList,
                   plot = FALSE,
                   plotType = 'points',
                   plotChanges = FALSE,
                   PCx = 1,
                   PCy = 2,
                   plotCartoon = FALSE,
                   refShape = NULL,
                   outline = NULL,
                   lines = NULL,
                   landList = NULL,
                   adjustCoords = FALSE,
                   crop = c(0,0,0,0),
                   flipRaster = NULL,
                   flipOutline = NULL,
                   imageList = NULL,
                   cartoonID = NULL,
                   colpalette = NULL,
                   normalized = NULL,
                   cartoonOrder = 'above',
                   lineOrder = 'above',
                   cartoonCol = 'gray',
                   cartoonFill = NULL,
                   plotLandmarks = FALSE,
                   landCol = 'black',
                   zlim = c(-1,1),
                   legendTitle = 'Predicted',
                   xlab='',
                   ylab='',
                   main=''){

  for(r in 1:length(rList)){

    rList[[r]][is.na(rList[[r]])] <- 0
    ras <- raster::as.data.frame(rList[[r]])
    colnames(ras) <- names(rList)[[r]]

    if(r == 1){
      rasDF <- ras
    }
    else{
      rasDF <- cbind(rasDF, ras)
    }
  }
  groupCol <- c()

  for(p in 1:length(popList)){

    for(ind in 1:length(popList[[p]])){

        groupCol <- rbind(groupCol, c(popList[[p]][ind], colList[p]))

    }
  }

  groupCol<-as.data.frame(groupCol)
  colnames(groupCol) <- c('sampleID', 'col')


  if(plot == TRUE){

    comp <- prcomp(t(rasDF))

    pcdata <- comp$x
    rotation <- comp$rotation

    summ <- summary(comp)

    if(plotChanges){

      PCxmin <- min(pcdata[,PCx])
      PCxmax <- max(pcdata[,PCx])

      PCymin <- min(pcdata[,PCy])
      PCymax <- max(pcdata[,PCy])


      pc.vecMix <- rep(0, dim(pcdata)[1])
      pc.vecMix[PCx] <- PCxmin

      pc.vecMax <- rep(0, dim(pcdata)[1])
      pc.vecMax[PCx] <- PCxmax

      pc.vecMiy <- rep(0, dim(pcdata)[1])
      pc.vecMiy[PCy] <- PCymin

      pc.vecMay <- rep(0, dim(pcdata)[1])
      pc.vecMay[PCy] <- PCymax


      xMi <- pc.vecMix %*%  t(rotation)
      xMa <- pc.vecMax %*%  t(rotation)

      x2Mi <- t(matrix(xMi,ncol = dim(rList[[1]])[1], nrow = dim(rList[[1]])[2]))
      x2Ma <- t(matrix(xMa,ncol = dim(rList[[1]])[1], nrow = dim(rList[[1]])[2]))

      yMi <- pc.vecMiy %*%  t(rotation)
      yMa <- pc.vecMay %*%  t(rotation)

      y2Mi <- t(matrix(yMi,ncol = dim(rList[[1]])[1], nrow = dim(rList[[1]])[2]))
      y2Ma <- t(matrix(yMa,ncol = dim(rList[[1]])[1], nrow = dim(rList[[1]])[2]))


      mapMix <-raster::raster(x2Mi)
      mapMax <-raster::raster(x2Ma)

      mapMiy <-raster::raster(y2Mi)
      mapMay <-raster::raster(y2Ma)

      raster::extent(mapMix) <- raster::extent(rList[[1]])
      raster::extent(mapMax) <- raster::extent(rList[[1]])

      raster::extent(mapMiy) <- raster::extent(rList[[1]])
      raster::extent(mapMay) <- raster::extent(rList[[1]])
    }

    if(!plotChanges){

      if(plotType == 'points'){

        plot(comp$x[,PCx:PCy], col=as.vector(groupCol$col), pch=20, cex=2)
      }

      if(plotType == 'labels'){

        plot(comp$x[,PCx:PCy], col=NA, pch=19)
        text(comp$x[,PCx:PCy], col=as.vector(groupCol$col), as.character(groupCol$sampleID))
      }
    }

    if(plotChanges){

      mat <- matrix(c(4,1,1,5,1,1,6,2,3), 3, 3, byrow = TRUE)
      layout(mat, widths=c(1,1,1), heights=c(1,1,1))


      # par(mar=c(4,4,1,1))

      if(plotType == 'points'){

        plot(comp$x[,PCx], comp$x[,PCy], col=as.vector(groupCol$col), pch=20, cex=3,
             xlab=paste('PC',PCx,' (', round(summ$importance[2,PCx]*100, 1), ' %)'),
             ylab=paste('PC',PCy,' (', round(summ$importance[2,PCy]*100, 1), ' %)'))
      }

      if(plotType == 'labels'){

        plot(comp$x[,PCx], comp$x[,PCy], col=NA, pch=19,
             xlab=paste('PC',PCx,' (', round(summ$importance[2,PCx]*100, 1), ' %)'),
             ylab=paste('PC',PCy,' (', round(summ$importance[2,PCy]*100, 1), ' %)'))

        text(comp$x[,PCx], comp$x[,PCy], col=as.vector(groupCol$col), as.character(groupCol$sampleID))
      }

      if(is.null(colpalette)){

        colpalette <- c("blue","lightblue","black","pink","red")
      }

      else{

        if(!is.vector(colpalette)){

          stop('Specified color palette is not a vector')

        }
      }

      plotHeat(mapMix, rList, plotCartoon = plotCartoon, refShape = refShape, outline = outline, lines = lines,
               adjustCoords = adjustCoords, landList = landList, crop = crop, flipRaster = flipRaster,
               flipOutline = flipOutline, imageList = imageList, cartoonID = cartoonID, colpalette = colpalette,
               normalized = normalized, cartoonOrder = cartoonOrder, lineOrder = lineOrder, cartoonCol = cartoonCol,
               cartoonFill = cartoonFill, plotLandmarks = plotLandmarks, landCol = landCol, zlim = zlim, xlab=xlab,
               ylab=ylab, main=main, plotPCA = TRUE)

      mtext(paste('min PC', PCx, sep=' '), 1)

      plotHeat(mapMax, rList, plotCartoon = plotCartoon, refShape = refShape, outline = outline, lines = lines,
               adjustCoords = adjustCoords, landList = landList, crop = crop, flipRaster = flipRaster,
               flipOutline = flipOutline, imageList = imageList, cartoonID = cartoonID, colpalette = colpalette,
               normalized = normalized, cartoonOrder = cartoonOrder, lineOrder = lineOrder, cartoonCol = cartoonCol,
               cartoonFill = cartoonFill, plotLandmarks = plotLandmarks, landCol = landCol, zlim = zlim, xlab=xlab,
               ylab=ylab, main=main, plotPCA = TRUE)

      mtext(paste('max PC', PCx, sep=' '), 1)

      plotHeat(mapMay, rList, plotCartoon = plotCartoon, refShape = refShape, outline = outline, lines = lines,
               adjustCoords = adjustCoords, landList = landList, crop = crop, flipRaster = flipRaster,
               flipOutline = flipOutline, imageList = imageList, cartoonID = cartoonID, colpalette = colpalette,
               normalized = normalized, cartoonOrder = cartoonOrder, lineOrder = lineOrder, cartoonCol = cartoonCol,
               cartoonFill = cartoonFill, plotLandmarks = plotLandmarks, landCol = landCol, zlim = zlim, xlab=xlab,
               ylab=ylab, main=main, plotPCA = TRUE)

      mtext(paste('max PC', PCy, sep=' '), 2)

      plotHeat(mapMiy, rList, plotCartoon = plotCartoon, refShape = refShape, outline = outline, lines = lines,
               adjustCoords = adjustCoords, landList = landList, crop = crop, flipRaster = flipRaster,
               flipOutline = flipOutline, imageList = imageList, cartoonID = cartoonID, colpalette = colpalette,
               normalized = normalized, cartoonOrder = cartoonOrder, lineOrder = lineOrder, cartoonCol = cartoonCol,
               cartoonFill = cartoonFill, plotLandmarks = plotLandmarks, landCol = landCol, zlim = zlim, xlab=xlab,
               ylab=ylab, main=main, plotPCA = TRUE)

      mtext(paste('min PC', PCy, sep=' '), 2)

      colfunc <- colorRampPalette(colpalette)

      plot(NULL, xlim=c(0,1), ylim=c(0,1), type="n", axes = FALSE, xlab = '', ylab='')

      plot(mapMay, col=colfunc(21), zlim = c(-1,1), legend.only = TRUE, legend.width = 5, horizontal = TRUE,
           smallplot = c(0.3, 1, 0.5, 0.6), legend.args = list(text=legendTitle, side = 3, font = 2, line = 2.5, cex = 1))
    }
  }

  return(list(t(rasDF), groupCol, comp))
}




