#' Plots the color pattern heatmaps from \code{sumRaster} output.
#'
#' @param summedRaster Summed raster or summedRasterList.
#' @param IDlist List of sample IDs.
#' @param colpalette Vector of colors for color palette
#'    (default = c("white","lightblue","blue","green", "yellow","red"))
#' @param plotCartoon Whether to plot a cartoon. This cartoon should be drawn on one of the samples
#'    used in the analysis.
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
#' @param imageList List of images should be given if one wants to flip the outline or adjust
#'    landmark coordinates.
#' @param refImage Image (RasterStack) used for target. Use raster::stack('filename').
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
#' @param legend Whether to plot legend with heatmaps.
#' @param legendTitle Title of the raster legend (default = 'Proportion')
#' @param legend.side Side to plot legend (default = 4)
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param main Optional main title.
#' @param plotType Set as 'PCA' when visualizing shape changes along PCA axis in \
#'    code{\link[patternize]{patPCA}}, as 'one' when visualizing single image or as 'multi' for multi
#'    plotting or when setting customized margins (default = 'multi').
#' @param imageIDs A list of IDs to match landmarks to images if landmarkList and imageList don't
#'    have the same length.
#' @param format ImageJ (Fiji) or tps format (default = 'imageJ').
#'
#' @examples
#' data(rasterList_lanRGB)
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' outline_BC0077 <- read.table(paste(system.file("extdata",  package = 'patternize'),
#' '/BC0077_outline.txt', sep=''), header = FALSE)
#' lines_BC0077 <- list.files(path=paste(system.file("extdata",  package = 'patternize')),
#' pattern='vein', full.names = TRUE)
#'
#' summedRaster_regRGB <- sumRaster(rasterList_regRGB, IDlist, type = 'RGB')
#' data(imageList)
#'
#' plotHeat(summedRaster_regRGB, IDlist, plotCartoon = TRUE, refShape = 'target',
#' outline = outline_BC0077, lines = lines_BC0077, crop = c(100,400,40,250),
#' flipRaster = 'xy', imageList = imageList, cartoonOrder = 'under', cartoonID = 'BC0077',
#' cartoonFill = 'black', main = 'registration_example')
#'
#' \dontrun{
#' data(rasterList_lanK)
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' summedRasterList <- sumRaster(rasterList_lanK, IDlist, type = 'k')
#' plotHeat(summedRasterList, IDlist)
#'
#' summedRasterList_regK <- sumRaster(rasterList_regK, IDlist, type = 'k')
#' plotHeat(summedRasterList_regK, IDlist, plotCartoon = TRUE, refShape = 'target',
#' outline = outline_BC0077, lines = lines_BC0077, crop = c(100,400,40,250),
#' flipRaster = 'y', imageList = imageList, cartoonOrder = 'under',
#' cartoonFill = 'black', main = 'kmeans_example')
#'
#' plotHeat(summedRasterList_regK[[1]], IDlist, plotCartoon = TRUE, refShape = 'target',
#' outline = outline_BC0077, lines = lines_BC0077, crop = c(100,400,40,250),
#' flipRaster = 'y', imageList = imageList, cartoonOrder = 'under',
#' cartoonFill = 'black', main = 'kmeans_example')
#'
#'
#' prepath <- system.file("extdata", package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#' landmarkList <- makeList(IDlist, 'landmark', prepath, extension)
#'
#' summedRaster_lanRGB <- sumRaster(rasterList_lanRGB, IDlist, type = 'RGB')
#'
#' plotHeat(summedRaster_lanRGB, IDlist, plotCartoon = TRUE, refShape = 'mean',
#' outline = outline_BC0077, lines = lines_BC0077, landList = landmarkList,
#' adjustCoords = TRUE, imageList = imageList, cartoonID = 'BC0077',
#' cartoonOrder = 'under', cartoonFill= 'black', main = 'Landmark_example')
#'
#' summedRaster_lanK <- sumRaster(rasterList_lanK, IDlist, type = 'k')
#'
#' plotHeat(summedRaster_lanK, IDlist, plotCartoon = TRUE, refShape = 'mean',
#' outline = outline_BC0077, lines = lines_BC0077, landList = landmarkList,
#' adjustCoords = TRUE, imageList = imageList, cartoonID = 'BC0077',
#' cartoonOrder = 'under', cartoonFill= 'black', main = 'Landmark_example')
#'
#' plotHeat(summedRaster_lanK[[2]], IDlist, plotCartoon = TRUE, refShape = 'mean',
#' outline = outline_BC0077, lines = lines_BC0077, landList = landmarkList,
#' adjustCoords = TRUE, imageList = imageList, cartoonID = 'BC0077',
#' cartoonOrder = 'under', cartoonFill= 'black', main = 'Landmark_example')
#' }
#'
#' @export
#' @import raster
#' @importFrom graphics layout mtext par points polygon
#' @importFrom grDevices col2rgb colorRampPalette rgb

plotHeat <- function(summedRaster,
                     IDlist,
                     colpalette = NULL,
                     plotCartoon = FALSE,
                     refShape = NULL,
                     outline = NULL,
                     lines = NULL,
                     landList = NULL,
                     adjustCoords = FALSE,
                     cartoonID = NULL,
                     normalized = FALSE,
                     crop = c(0,0,0,0),
                     flipRaster = NULL,
                     flipOutline = NULL,
                     imageList = NULL,
                     refImage = NULL,
                     cartoonOrder = 'above',
                     lineOrder = 'above',
                     cartoonCol = 'gray',
                     cartoonFill = NULL,
                     plotLandmarks = FALSE,
                     landCol = 'black',
                     zlim = c(0,1),
                     legend = TRUE,
                     legendTitle = 'Proportion',
                     legend.side = 4,
                     xlab='',
                     ylab='',
                     main='',
                     plotType = 'multi',
                     imageIDs = NULL,
                     format = 'imageJ'){



  if(!is.list(summedRaster)){
    rasterEx <- raster::extent(summedRaster)
  }
  else{
    rasterEx <- raster::extent(summedRaster[[1]])
  }

  if(is.null(cartoonID)){
    if(!is.list(summedRaster)){
      imageEx <- raster::extent(summedRaster)
    }
    else{
      imageEx <- raster::extent(summedRaster[[1]])
    }
  }
  else{
    if(cartoonID %in% names(imageList)){
      imageEx <- raster::extent(imageList[[cartoonID]])
    }
    else{
      imageEx <- raster::extent(imageList[[1]])
    }
    if(!is.null(refImage)){
      imageEx <- raster::extent(refImage)
    }
  }


  if(!is.null(lines)){

    lineList <- list()

    for(e in 1:length(lines)){

      lineList[[e]] <- read.table(lines[e], header = FALSE)

    }
  }

  if(is.null(landList)){
    # need to fix the outline for the shift that happens in the raster extent when flipping.
    # Not an issue if crop wasn't used.
    exDiff1 <- abs(rasterEx[2]-imageEx[2])-rasterEx[1]
    outline[,1] <- outline[,1] + exDiff1

    # exDiff2 <- abs(rasterEx[4]-imageEx[4])-rasterEx[3]
    # outline[,2] <- outline[,2] + exDiff2

    if(!is.null(lines)){
      for(e in 1:length(lineList)){

        lineList[[e]][[1]] <- lineList[[e]][[1]] + exDiff1
        # lineList[[e]][[2]] <- lineList[[e]][[2]] + exDiff2
      }
    }
  }

  if(!is.null(cartoonID)){
    if(cartoonID %in% names(imageList)){
      imageExRef <- raster::extent(imageList[[cartoonID]])
    }
  }
  else if(!is.null(refImage)){
    imageExRef <- raster::extent(refImage)
  }
  else{
    if(plotCartoon != FALSE){
      stop("if cartoonID not in image list, please provide image to refImage argument.")
    }
  }

  if(format == 'tps'){
    outline[,2] <- (imageExRef[4]-outline[,2])
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


  if(normalized){
    divide <- 1
  }
  else{
    divide <- length(IDlist)
  }

  if(any(c(!is.null(flipOutline), !is.null(flipRaster)))){


    if(all(c(refShape[1] != 'mean', !is.matrix(refShape)))){

      outline[,2] <- outline[,2] - crop[3]

      if(!is.null(lines)){

        for(e in 1:length(lineList)){

          lineList[[e]][[2]] <- lineList[[e]][[2]] - crop[3]
        }
      }
    }
    if(all(c(refShape[1] != 'mean', !is.matrix(refShape), identical(crop, c(0,0,0,0))))){
      if(all(c(!is.null(flipOutline), flipOutline == 'y'))){
        outline[,2] <- imageEx[4] - outline[,2]

        if(!is.null(lines)){
          for(e in 1:length(lineList)){

            lineList[[e]][[2]] <- imageEx[4] - lineList[[e]][[2]]
          }
        }
      }
    }

  }

  if(!identical(crop, c(0,0,0,0))){
    if(all(c(refShape[1] != 'mean', !is.matrix(refShape)))){

      if(all(c(!is.null(flipOutline), flipOutline == 'y'))){
        outline[,2] <- outline[,2] + crop[3]
        if(!is.null(lines)){
          for(e in 1:length(lineList)){
            lineList[[e]][[2]] <- lineList[[e]][[2]] + crop[3]
          }
        }
      }
      if(all(c(!is.null(flipOutline), flipOutline == 'xy'))){
        outline[,2] <- outline[,2] + crop[3]
        if(!is.null(lines)){
          for(e in 1:length(lineList)){
            lineList[[e]][[2]] <- lineList[[e]][[2]] + crop[3]
          }
        }
      }
      # if(is.null(flipOutline) && !is.null(flipRaster) || !is.null(flipOutline) && flipOutline == 'x'){
      if(all(c(!is.null(flipOutline), flipOutline == 'x'))){
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

  if(all(c(plotCartoon, refShape[1] != 'target'))){

    indx <- which(names(imageList) == cartoonID)
    invisible(capture.output(landArray <- lanArray(landList, adjustCoords, imageList, imageIDs = imageIDs)))

    if(is.character(refShape)){
      indx <- which(names(landList) == cartoonID)
      transRefLan <- as.matrix(lanArray[,,indx])
    }
    else{
      transRefLan <- refShape
    }

    if(is.matrix(refShape)){
      invisible(capture.output(cartoonLandTrans <- Morpho::computeTransform(refShape,
                                                                            as.matrix(transRefLan),
                                                                            type="tps")))
    }

    if(!is.matrix(refShape)){

      invisible(capture.output(transformed <- Morpho::procSym(landArray)))

      refShape <- transformed$mshape

      invisible(capture.output(cartoonLandTrans <- Morpho::computeTransform(refShape,
                                                                            as.matrix(transRefLan),
                                                                            type="tps")))
    }

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


    outlineTrans <- Morpho::applyTransform(as.matrix(outline), cartoonLandTrans)

    cartoonLinesTrans <- NULL

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
      if(any(c(refShape[1] == 'target', is.matrix(refShape)))){
        XLIM <- c(min(outline[,1]),max(outline[,1]))
        YLIM <- c(min(outline[,2]),max(outline[,2]))
      }
    }

    if(!is.null(flipRaster)){

      if(!is.null(refImage)){
        y <- raster::raster(ncol = dim(refImage)[2], nrow = dim(refImage)[1])
        raster::extent(y) <- raster::extent(refImage)
        summedRaster <- raster::resample(summedRaster, y)
      }

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

    if(plotType == 'one'){
      par(mfrow = c(1,1), mai=c(0.05,0.8,0.15,0.8), oma=c(1,1,1,1)+1)
    }
    if(plotType == 'PCA'){
      par(mar=c(4,4,2,2))
    }

    if(any(c(is.null(refShape), refShape[1] == 'target'))){
      plot(NULL, type="n", axes=F, xlim = XLIM, ylim= YLIM, main=main, xlab = '', ylab='', asp=1)
    }

    if(plotCartoon){

      if(any(c(is.null(refShape), is.null(outline)))){

        stop('Not all paramters are set to plot the cartoon.')

      }
    }

    if(plotCartoon && cartoonOrder == 'under'){

      summedRaster[summedRaster == 0] <- NA

      if(refShape[1] == 'target'){

        polygon(outline, col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

      }

      if(any(c(refShape[1] == 'mean', is.matrix(refShape)))){

        plot(NULL, type="n", axes=F, xlim = XLIM, ylim= YLIM, main=main, xlab = '', ylab='', asp=1)

        polygon(outlineTrans,col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

      }
    }

    if(all(c(plotCartoon, lineOrder == 'under'))){

      if(refShape[1] == 'target'){

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lines(lineList[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

          }
        }
      }

      if(any(c(refShape[1] =='mean', is.matrix(refShape)))){

        if(!is.null(lines)){

          for(e in 1:length(cartoonLinesTrans)){

            lines(cartoonLinesTrans[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)
          }
        }
      }
    }


    if(plotType == 'PCA'){

      image(summedRaster/divide, col=colfunc(99), xaxt='n', yaxt='n', box=F, axes=F,
            xlim = XLIM, ylim= YLIM, zlim=zlim, add= TRUE, useRaster= FALSE, legend = FALSE, asp=1)
    }
    else{

      plot(summedRaster/divide, col=colfunc(99), xaxt='n', yaxt='n', box=F, axes=F,
           xlim = XLIM, ylim= YLIM, zlim=zlim, add= TRUE, useRaster= FALSE, legend = FALSE, asp=1)

      if(legend == TRUE){
        plot(summedRaster/divide, legend.only=TRUE, zlim=zlim, col=colfunc(99),legend.width=1, legend.shrink=0.75,
             legend.args=list(text=legendTitle, side=legend.side, font=2, line=2.5, cex=1))
        }
      }
    mtext(side = 1, text = xlab, line = 0)
    mtext(side = 2, text = ylab, line = 0)

    if(all(c(plotCartoon, cartoonOrder == 'above'))){

      if(refShape[1] == 'target'){

        polygon(outline, col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

      }

      if(any(c(refShape[1] == 'mean', is.matrix(refShape)))){

        polygon(outlineTrans,col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

      }
    }

    if(all(c(plotCartoon, lineOrder == 'above'))){

      if(refShape[1] == 'target'){

        if(!is.null(lines)){

          for(e in 1:length(lineList)){

            lines(lineList[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

          }
        }
      }
      if(any(c(refShape[1] =='mean', is.matrix(refShape)))){

        if(!is.null(lines)){

          for(e in 1:length(cartoonLinesTrans)){

            lines(cartoonLinesTrans[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)
          }
        }
      }
    }

    if(plotLandmarks){
      if(refShape[1] == 'target'){
        points(as.matrix(landArray[,,indx]), pch = 19, col = landCol)
      }
      if(any(c(refShape[1] =='mean', is.matrix(refShape)))){
        points(transformed$mshape, pch = 19, col = landCol)
      }
    }
  }

  else{

    if(is.null(outline)){
      XLIM <- c(rasterEx[1],rasterEx[2])
      YLIM <- c(rasterEx[3],rasterEx[4])
    }
    else{
      if(any(c(refShape[1] == 'target', is.matrix(refShape)))){
        XLIM <- c(min(outline[,1]),max(outline[,1]))
        YLIM <- c(min(outline[,2]),max(outline[,2]))
      }
    }


    par(mfrow=c(2,trunc((length(summedRaster)+1)/2)), mai=c(0.05,0.8,0.15,0.8), oma=c(1,1,1,1)+1)

    for(k in 1:length(summedRaster)){

      if(!is.null(flipRaster)){

        if(!is.null(refImage)){
          y <- raster::raster(ncol = dim(refImage)[2], nrow = dim(refImage)[1])
          raster::extent(y) <- raster::extent(refImage)
          summedRaster[[k]] <- raster::resample(summedRaster[[k]], y)
        }

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

      if(any(c(is.null(refShape), refShape[1] == 'target'))){
        plot(NULL, type="n", axes=F, xlab="", ylab="", xlim = XLIM, ylim= YLIM, main= main, asp=1)

      }

      if(plotCartoon){

        if(any(c(is.null(refShape), is.null(outline)))){

          stop('Not all paramters are set to plot the cartoon.')

        }
      }

      if(plotCartoon && cartoonOrder == 'under'){

        summedRaster[[k]][summedRaster[[k]] == 0] <- NA

        if(refShape[1] == 'target'){

          polygon(outline, col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

        }

        if(any(c(refShape[1] == 'mean', is.matrix(refShape)))){

          plot(NULL, type="n", axes=F, xlim = XLIM, ylim= YLIM, main=main, xlab = '', ylab='', asp=1)

          polygon(outlineTrans,col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

        }
      }

      if(all(c(plotCartoon, lineOrder == 'under'))){

        if(refShape[1] == 'target'){

          if(!is.null(lines)){

            for(e in 1:length(lineList)){

              lines(lineList[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

            }
          }
        }

        if(any(c(refShape[1] =='mean', is.matrix(refShape)))){

          if(!is.null(lineList)){

            for(e in 1:length(cartoonLinesTrans)){

              lines(cartoonLinesTrans[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)
            }
          }
        }
      }

      plot(summedRaster[[k]]/divide, col=colfunc(99), xaxt='n', yaxt='n', box=F, axes=F,
           xlim = XLIM, ylim= YLIM, zlim=zlim, legend.args=list(text=legendTitle, side=legend.side, line=3),
           add= TRUE, asp=1)

      mtext(side = 1, text = xlab, line = 0)
      mtext(side = 2, text = ylab, line = 0)

      if(all(c(plotCartoon, cartoonOrder == 'above'))){

        if(refShape[1] == 'target'){

          polygon(outline, col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

        }

        if(any(c(refShape[1] == 'mean', is.matrix(refShape)))){

          polygon(outlineTrans,col=cartoonFill, border=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

        }
      }

      if(all(c(plotCartoon, lineOrder == 'above'))){

        if(refShape[1] == 'target'){

          if(!is.null(lines)){

            for(e in 1:length(lineList)){

              lines(lineList[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)

            }
          }
        }

        if(any(c(refShape[1] =='mean', is.matrix(refShape)))){

          if(!is.null(lines)){

            for(e in 1:length(cartoonLinesTrans)){

              lines(cartoonLinesTrans[[e]], col=cartoonCol, xlim = XLIM, ylim= YLIM, asp=1)
            }
          }

          if(plotLandmarks){

            points(transformed$mshape, pch = 19, col = landCol)
          }
        }
      }

      if(plotLandmarks){
        if(refShape[1] == 'target'){
          points(as.matrix(landArray[,,indx]), pch = 19, col = landCol)
        }
        if(any(c(refShape[1] =='mean', is.matrix(refShape)))){
          points(transformed$mshape, pch = 19, col = landCol)
        }
      }
    }
  }
}
