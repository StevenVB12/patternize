#' This function transforms the individual color pattern rasters as obtained by the main
#' patternize functions to a dataframe of 0 and 1 values that can be used for constrained
#' Redundancy Analysis (RDA) (\code{\link[vegan]{rda}}). This function also allows to plot the
#' analysis including a visualization of the shape changes along the axis.
#'
#' @param rList List of raster objects.
#' @param popList List of vectors including sampleIDs for each population.
#' @param colList List of colors for each population.
#' @param symbolList List with graphical plotting symbols (default = NULL).
#' @param rListPredict List of raster objects to predict into DFA space (default = NULL).
#' @param popListPredict List of vectors including sampleIDs for each set of predict samples
#' (default = NULL). Note to that this also has to be a list if only one population is included.
#' @param colListPredict List of colors for each set of predict samples (default = NULL).
#' @param symbolListPredict List with graphical plotting symbols for predict sets (default = NULL).
#' @param plot Whether to plot the PCA analysis (default = FALSE).
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
#' @param landCol Color for plotting landmarks (default = 'black').
#' @param zlim z-axis limit (default = c(0,1))
#' @param legendTitle Title of the raster legend (default = 'Proportion')
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param main Optional main title.
#'
#' @return  If plot = TRUE: List including a [1] dataframe of the binary raster values that can be used for
#'    discriminant function analysis, [2] a dataframe of sample IDs and specified population
#'    colors and [3] lda results. if rListPredict not empty: [4] class prediction of samples. If plot = FALSE:
#'    lda result only.
#'
#' @seealso \code{\link[MASS]{lda}}
#'
#' @examples
#' data(rasterList_lanRGB)
#'
#' pop1 <- c('BC0077','BC0071')
#' pop2 <- c('BC0050','BC0049','BC0004')
#' popList <- list(pop1, pop2)
#' colList <- c("red", "blue")
#'
#' pcaOut <- patRDA(rasterList_lanRGB, popList, colList, plot = TRUE)
#'
#' @export
#' @import raster vegan
#' @importFrom grDevices adjustcolor
#' @importFrom stats var

patRDA <- function(rList,
                   popList,
                   colList,
                   symbolList = NULL,
                   rListPredict = NULL,
                   popListPredict = NULL,
                   colListPredict = NULL,
                   symbolListPredict = NULL,
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

  # data for DFA
  # make dataframe of rasters
  print("making dataframe from rasters")
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

  # remove non-variant rows in each population
  print("removing non-variant rows")
  popSub <- list()
  rwN <- list()

  # 1. remove invariant rows over all data

  popSub1 <- rasDF[apply(rasDF, 1, var, na.rm=TRUE) != 0,]
  cols <- rownames(popSub1)

  ## This is for lda
  # 2. add jitter to nonvariant rows within populations
  # for(e in 1:length(popList)){
  #   popSub[[e]] <- popSub1[,popList[[e]]]
  #   cols <- ncol(popSub[[e]])
  #   for(l in 1:nrow(popSub[[e]])){
  #     if(mean(as.numeric(popSub[[e]][l,])) == 0){
  #       list0 <- rep(0,cols)
  #       indx <- sample(1:cols, 1)
  #       list0[indx] <- 0.001
  #       popSub[[e]][l,] <- list0
  #     }
  #     if(mean(as.numeric(popSub[[e]][l,])) == 1){
  #       list0 <- rep(1,cols)
  #       indx <- sample(1:cols, 1)
  #       list0[indx] <- 0.999
  #       popSub[[e]][l,] <- list0
  #     }
  #   }
  #   # popSub[[e]] <- popSub1[,popList[[e]]][apply(popSub1[,popList[[e]]], 1, var, na.rm=TRUE) != 0,]
  #   # rwN[[e]] <- rownames(popSub[[e]])
  # }
  #
  # # cols <- Reduce(intersect, rwN)
  #
  # for(e in 1:length(popSub)){
  #
  #   if(e == 1){pcaInVar <- popSub[[e]]}
  #   else{pcaInVar <- cbind(pcaInVar, popSub[[e]])}
  # }

  # make group vector
  grp <- c()

  for(p in 1:length(popList)){

    grp <- c(grp, rep(p, length(popList[[p]])))
  }
  grpM <- as.data.frame(cbind(colnames(popSub1),grp))
  grpM <- as.factor(as.character(grpM[,2]))

  # make color and symbol object
  groupCol <- c()

  for(p in 1:length(popList)){

    for(ind in 1:length(popList[[p]])){

      if(!is.null(symbolList)){

        groupCol <- rbind(groupCol, c(popList[[p]][ind], colList[p], symbolList[p]))
      }

      if(is.null(symbolList)){

        groupCol <- rbind(groupCol, c(popList[[p]][ind], colList[p]))
      }

    }
  }

  groupCol <- as.data.frame(groupCol)

  if(!is.null(symbolList)){
    colnames(groupCol) <- c('sampleID', 'col', 'symbol')
  }
  if(is.null(symbolList)){
    colnames(groupCol) <- c('sampleID', 'col')
  }


  print("calculating rda")
  # suppressWarnings(ldaOut <- lda(t(pcaInVar), grp))
  # plda <- predict(object = ldaOut, newdata = t(pcaInVar))
  # pcdata <- plda$x

  ldaOut <- vegan::rda(t(popSub1), grpM)
  pcdata <- predict(object = ldaOut, x=t(popSub1), y=grpM,  type = "wa")

  sc <- scores(ldaOut, choices = c(1:ncol(pcdata)))
  # pcdata <- sc$sites


  if(length(popList) > 2){
    xmin <- min(pcdata[,PCx])
    xmax <- max(pcdata[,PCx])
    ymin <- min(pcdata[,PCy])
    ymax <- max(pcdata[,PCy])
  }
  if(length(popList) <= 2){
    xmin <- min(pcdata[,1])
    xmax <- max(pcdata[,1])
  }


  if(!is.null(rListPredict)){
    # data for DFA predict
    # make dataframe of rasters for predict
    print("making predict dataframe from rasters")
    for(r in 1:length(rListPredict)){

      rListPredict[[r]][is.na(rListPredict[[r]])] <- 0
      ras <- raster::as.data.frame(rListPredict[[r]])
      colnames(ras) <- names(rListPredict)[[r]]

      if(r == 1){
        rasDFPredict <- ras
      }
      else{
        rasDFPredict <- cbind(rasDFPredict, ras)
      }
    }

    # remove non-variant rows in each population
    print("removing non-variant rows")

    groupColPredict <- c()

    for(p in 1:length(popListPredict)){

      for(ind in 1:length(popListPredict[[p]])){

        if(!is.null(symbolListPredict)){

          groupColPredict <- rbind(groupColPredict, c(popListPredict[[p]][ind], colListPredict[p], symbolListPredict[p]))
        }

        if(is.null(symbolListPredict)){

          groupColPredict <- rbind(groupColPredict, c(popListPredict[[p]][ind], colListPredict[p]))
        }

      }
    }

    groupColPredict <- as.data.frame(groupColPredict)

    if(!is.null(symbolListPredict)){
      colnames(groupColPredict) <- c('sampleID', 'col', 'symbol')
    }
    if(is.null(symbolListPredict)){
      colnames(groupColPredict) <- c('sampleID', 'col')
    }

    # pldaPredict <- predict(object = ldaOut, newdata = t(rasDFPredict[cols,]))
    # pcdataPredict <- pldaPredict$x

    pcdataPredict <- predict(ldaOut, t(rasDFPredict[cols,]), type="wa")

    if(length(popList) > 2){
      xmin <- min(pcdata[,PCx], pcdataPredict[,PCx])
      xmax <- max(pcdata[,PCx], pcdataPredict[,PCx])
      ymin <- min(pcdata[,PCy], pcdataPredict[,PCy])
      ymax <- max(pcdata[,PCy], pcdataPredict[,PCy])
    }
    if(ncol(pcdata) <= 2){
      xmin <- min(pcdata[,1], pcdataPredict[,1])
      xmax <- max(pcdata[,1], pcdataPredict[,1])
    }
  }


  if(plot == TRUE){

    # amount of the between-group variance that is explained by each linear discriminant
    # importance = ldaOut$svd^2/sum(ldaOut$svd^2)

    importance <- summary(eigenvals(ldaOut))


    if(plotChanges){
      print("calculating changes")

      PCxmin <- min(pcdata[,PCx])
      PCxmax <- max(pcdata[,PCx])

      if(length(popList) > 2){
        PCymin <- min(pcdata[,PCy])
        PCymax <- max(pcdata[,PCy])
      }


      pc.vecMix <- rep(0, ncol(sc$species))
      pc.vecMix[PCx] <- PCxmin

      pc.vecMax <- rep(0, ncol(sc$species))
      pc.vecMax[PCx] <- PCxmax

      if(length(popList) > 2){
        pc.vecMiy <- rep(0, ncol(sc$species))
        pc.vecMiy[PCy] <- PCymin

        pc.vecMay <- rep(0, ncol(sc$species))
        pc.vecMay[PCy] <- PCymax
      }

      # fill in empty rows
      rowN <- rownames(sc$species)
      colNr <- ncol(sc$species)
      rowNr <- nrow(rasDF)
      ldaF <- matrix(0, nrow = rowNr, ncol =colNr)
      n <- 1
      for(e in 1:nrow(rasDF)){
        if(n <= length(rowN)){
          line <- sc$species[n,]
          if(rowN[n] == e){
            ldaF[e,] <- line
            n <- n + 1
          }
          # else{
          #   ldaF[e,] <- rep(0, colN)
          # }
        }
        # else{
        #   ldaF <- rep(0, colN)
        # }
      }
      # ldaF <- as.data.frame(ldaF)
      rownames(ldaF) <- c(1:nrow(rasDF))


      xMi <- pc.vecMix %*%  t(ldaF)
      xMa <- pc.vecMax %*%  t(ldaF)

      x2Mi <- t(matrix(xMi,ncol = dim(rList[[1]])[1], nrow = dim(rList[[1]])[2]))
      x2Ma <- t(matrix(xMa,ncol = dim(rList[[1]])[1], nrow = dim(rList[[1]])[2]))

      if(length(popList) > 2){
        yMi <- pc.vecMiy %*%  t(ldaF)
        yMa <- pc.vecMay %*%  t(ldaF)

        y2Mi <- t(matrix(yMi,ncol = dim(rList[[1]])[1], nrow = dim(rList[[1]])[2]))
        y2Ma <- t(matrix(yMa,ncol = dim(rList[[1]])[1], nrow = dim(rList[[1]])[2]))
      }

      mapMix <-raster::raster(x2Mi)
      mapMax <-raster::raster(x2Ma)

      if(length(popList) > 2){
        mapMiy <-raster::raster(y2Mi)
        mapMay <-raster::raster(y2Ma)
      }

      raster::extent(mapMix) <- raster::extent(rList[[1]])
      raster::extent(mapMax) <- raster::extent(rList[[1]])

      if(length(popList) > 2){
        raster::extent(mapMiy) <- raster::extent(rList[[1]])
        raster::extent(mapMay) <- raster::extent(rList[[1]])
      }
    }

    print("plotting")
    if(plotChanges){

      mat <- matrix(c(4,1,1,5,1,1,6,2,3), 3, 3, byrow = TRUE)
      layout(mat, widths=c(1,1,1), heights=c(1,1,1))

      if(length(popList) <= 2){
        mat <- matrix(c(1,1,1,1,1,1,4,2,3), 3, 3, byrow = TRUE)
        layout(mat, widths=c(1,1,1), heights=c(1,1,1))
      }
    }

    if(!plotChanges){
      mat <- matrix(c(1,1,1,1,1,1,1,1,1), 3, 3, byrow = TRUE)
      layout(mat, widths=c(1,1,1), heights=c(1,1,1))
    }

    if(length(popList) > 2){
      if(all(c(plotType == 'points', is.null(symbolList)))){

        plot(pcdata[,c(PCx,PCy)], col=groupCol$col, pch=20, cex=3,
             xlim = c(xmin, xmax), ylim = c(ymin, ymax),
             xlab=paste('rda',PCx,' (', round(importance[2,PCx]*100, 1), ' %)'),
             ylab=paste('rda',PCy,' (', round(importance[2,PCy]*100, 1), ' %)'))

        if(!is.null(rListPredict)){
          points(pcdataPredict[,c(PCx,PCy)], col = groupColPredict$col, pch=20, cex=3)
        }
      }

      if(all(c(plotType == 'points', !is.null(symbolList)))){

        plot(pcdata[,c(PCx,PCy)], col=groupCol$col, pch=as.numeric(groupCol$symbol), cex=3,
             xlim = c(xmin, xmax), ylim = c(ymin, ymax),
             xlab=paste('rda',PCx,' (', round(importance[2,PCx]*100, 1), ' %)'),
             ylab=paste('rda',PCy,' (', round(importance[2,PCy]*100, 1), ' %)'))

        if(!is.null(rListPredict)){
          points(pcdataPredict[,c(PCx,PCy)], col = groupColPredict$col, pch=as.numeric(groupColPredict$symbol), cex=3)
        }

      }

      if(plotType == 'labels'){

        plot(pcdata[,c(PCx,PCy)], col=NA, pch=19,
             xlim = c(xmin, xmax), ylim = c(ymin, ymax),
             xlab=paste('rda',PCx,' (', round(importance[2,PCx]*100, 1), ' %)'),
             ylab=paste('rda',PCy,' (', round(importance[2,PCy]*100, 1), ' %)'))
        text(pcdata[,c(PCx,PCy)], col=groupCol$col, as.character(groupCol$sampleID))

        if(!is.null(rListPredict)){
          text(pcdataPredict[,c(PCx,PCy)], col = groupColPredict$col, as.character(groupColPredict$sampleID))
        }
      }
    }

    # plotting if only two groups
    if(length(popList) <= 2){

      d1 <- density(pcdata[,1][grp == 1])
      d2 <- density(pcdata[,1][grp == 2])

      if(all(c(plotType == 'points', is.null(symbolList)))){

        plot(pcdata[,1], rep(0, nrow(pcdata)), col=groupCol$col, pch=20, cex=3,
             xlim = c(xmin, xmax), ylim = c(0, max(1, d1$y,d2$y)),
             xlab=paste('rda',PCx,' (', round(importance[2,PCx]*100, 1), ' %)'),
             ylab = "Density")
        polygon(d1, col=adjustcolor(colList[1], alpha.f = 0.2), border=colList[1])
        polygon(d2, col=adjustcolor(colList[2], alpha.f = 0.2), border=colList[2])

        if(!is.null(rListPredict)){
          points(pcdataPredict[,1], rep(0, nrow(pcdataPredict)), col = groupColPredict$col, pch=20, cex=3)
        }
      }

      if(all(c(plotType == 'points', !is.null(symbolList)))){

        plot(pcdata[,1], rep(0, nrow(pcdata)), col=groupCol$col, pch=as.numeric(groupCol$symbol), cex=3,
             xlim = c(xmin, xmax), ylim = c(0, max(1, d1$y,d2$y)),
             xlab=paste('rda',PCx,' (', round(importance[2,PCx]*100, 1), ' %)'),
             ylab = "Density")
        polygon(d1, col=adjustcolor(colList[1], alpha.f = 0.2), border=colList[1])
        polygon(d2, col=adjustcolor(colList[2], alpha.f = 0.2), border=colList[2])

        if(!is.null(rListPredict)){
          points(pcdataPredict[,1], rep(0, nrow(pcdataPredict)), col = groupColPredict$col, pch=as.numeric(groupColPredict$symbol), cex=3)
        }

      }

      if(plotType == 'labels'){

        plot(pcdata[,1], rep(0, nrow(pcdata)), col=NA, pch=19,
             xlim = c(xmin, xmax), ylim = c(0, max(1, d1$y,d2$y)),
             xlab=paste('rda',PCx,' (', round(importance[2,PCx]*100, 1), ' %)'),
             ylab = "Density")
        polygon(d1, col=adjustcolor(colList[1], alpha.f = 0.2), border=colList[1])
        polygon(d2, col=adjustcolor(colList[2], alpha.f = 0.2), border=colList[2])

        text(pcdata[,1], rep(0, nrow(pcdata)), col=groupCol$col, as.character(groupCol$sampleID), srt = 45, adj = c(0,0))

        if(!is.null(rListPredict)){
          text(pcdataPredict[,1], rep(0, nrow(pcdataPredict)), col = groupColPredict$col, as.character(groupColPredict$sampleID), srt = 45, adj = c(0,0))
        }
      }
    }

    if(plotChanges){
      if(is.null(colpalette)){

        colpalette <- c("blue","lightblue","black","pink","red")
      }

      else{

        if(!is.vector(colpalette)){

          stop('Specified color palette is not a vector')

        }
      }

      plotHeat(mapMix/max(abs(xMi)), rList, plotCartoon = plotCartoon, refShape = refShape, outline = outline, lines = lines,
               adjustCoords = adjustCoords, landList = landList, crop = crop, flipRaster = flipRaster,
               flipOutline = flipOutline, imageList = imageList, cartoonID = cartoonID, colpalette = colpalette,
               normalized = normalized, cartoonOrder = cartoonOrder, lineOrder = lineOrder, cartoonCol = cartoonCol,
               cartoonFill = cartoonFill, plotLandmarks = plotLandmarks, landCol = landCol, zlim = zlim, xlab=xlab,
               ylab=ylab, main=main, plotType = 'PCA')

      mtext(paste('min DF', PCx, sep=' '), 1)

      plotHeat(mapMax/max(abs(xMa)), rList, plotCartoon = plotCartoon, refShape = refShape, outline = outline, lines = lines,
               adjustCoords = adjustCoords, landList = landList, crop = crop, flipRaster = flipRaster,
               flipOutline = flipOutline, imageList = imageList, cartoonID = cartoonID, colpalette = colpalette,
               normalized = normalized, cartoonOrder = cartoonOrder, lineOrder = lineOrder, cartoonCol = cartoonCol,
               cartoonFill = cartoonFill, plotLandmarks = plotLandmarks, landCol = landCol, zlim = zlim, xlab=xlab,
               ylab=ylab, main=main, plotType = 'PCA')

      mtext(paste('max DF', PCx, sep=' '), 1)

      if(length(popList) > 2){
        plotHeat(mapMay/max(abs(yMa)), rList, plotCartoon = plotCartoon, refShape = refShape, outline = outline, lines = lines,
                 adjustCoords = adjustCoords, landList = landList, crop = crop, flipRaster = flipRaster,
                 flipOutline = flipOutline, imageList = imageList, cartoonID = cartoonID, colpalette = colpalette,
                 normalized = normalized, cartoonOrder = cartoonOrder, lineOrder = lineOrder, cartoonCol = cartoonCol,
                 cartoonFill = cartoonFill, plotLandmarks = plotLandmarks, landCol = landCol, zlim = zlim, xlab=xlab,
                 ylab=ylab, main=main, plotType = 'PCA')

        mtext(paste('max DF', PCy, sep=' '), 2)

        plotHeat(mapMiy/max(abs(yMi)), rList, plotCartoon = plotCartoon, refShape = refShape, outline = outline, lines = lines,
                 adjustCoords = adjustCoords, landList = landList, crop = crop, flipRaster = flipRaster,
                 flipOutline = flipOutline, imageList = imageList, cartoonID = cartoonID, colpalette = colpalette,
                 normalized = normalized, cartoonOrder = cartoonOrder, lineOrder = lineOrder, cartoonCol = cartoonCol,
                 cartoonFill = cartoonFill, plotLandmarks = plotLandmarks, landCol = landCol, zlim = zlim, xlab=xlab,
                 ylab=ylab, main=main, plotType = 'PCA')

        mtext(paste('min DF', PCy, sep=' '), 2)
      }

      colfunc <- colorRampPalette(colpalette)

      plot(NULL, xlim=c(0,1), ylim=c(0,1), type="n", axes = FALSE, xlab = '', ylab='')

      plot(mapMax, col=colfunc(21), zlim = c(-1,1), legend.only = TRUE, legend.width = 5, horizontal = TRUE,
           smallplot = c(0.3, 1, 0.5, 0.6), legend.args = list(text=legendTitle, side = 3, font = 2, line = 2.5, cex = 1))
    }

    if(is.null(rListPredict)){
      return(list(t(rasDF), groupCol, ldaOut))
    }
    if(!is.null(rListPredict)){
      return(list(t(rasDF), groupCol, ldaOut, pcdataPredict))
    }

  }

  else{
    return(ldaOut)
  }


}
