#' Plot color pattern prediction for specified PCA values
#'
#' @param PCAdata Output of PCA analysis. List item 3 of patPCA.
#' @param PCApredict A vector with the PCA values for which to predict the phenotype. This vector
#'    only needs to include the values upto the last PCA axis to predict along, other values are
#'    set to zero.
#' @param IDlist List of sample IDs.
#' @param rasterList rasterList used for PCA.
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
#' @param cartoonOrder Whether to plot the cartoon outline 'above' or 'under' the pattern raster
#'    (default = 'above'). Set to 'under' for filled outlines.
#' @param lineOrder Whether to plot the cartoon lines 'above' or 'under' the pattern raster
#'    (default = 'above').
#' @param cartoonCol Outline and line color for cartoon (deafault = 'gray').
#' @param cartoonFill Fill color for outline of cartoon (default = NULL).
#' @param legendTitle Title of the raster legend (default = 'Proportion').
#' @param zlim zlim values for predicted pattern.
#'
#' @export
#' @import raster

createPhenotype <- function(PCAdata,
                            PCApredict,
                            IDlist,
                            rasterList,
                            colpalette = NULL,
                            plotCartoon = FALSE,
                            refShape = NULL,
                            outline = NULL,
                            lines = NULL,
                            landList = NULL,
                            adjustCoords = FALSE,
                            cartoonID = NULL,
                            normalized = TRUE,
                            crop = c(0,0,0,0),
                            flipRaster = NULL,
                            flipOutline = NULL,
                            imageList = NULL,
                            cartoonOrder = 'above',
                            lineOrder = 'above',
                            cartoonCol = 'gray',
                            cartoonFill = NULL,
                            legendTitle = 'Proportion',
                            zlim = NULL){

  pc <- PCAdata$x
  rotation <- PCAdata$rotation

  pc.vec <- rep(0, dim(pc)[1])
  pc.vec[1:length(PCApredict)] <- PCApredict

  pc.pred <- pc.vec %*%  t(rotation)

  pc.pred.image <- t(matrix(pc.pred, ncol = dim(rasterList[[1]])[1], nrow = dim(rasterList[[1]])[2]))

  pc.pred.image.raster <-raster::raster(pc.pred.image)

  raster::extent(pc.pred.image.raster) <- raster::extent(rasterList[[1]])

  if(is.null(zlim)){
    plotHeat(pc.pred.image.raster, IDlist, plotCartoon = TRUE, refShape = refShape, outline = outline, lines = lines,
             adjustCoords = TRUE, landList = landList, imageList = imageList, cartoonID = cartoonID, colpalette = colpalette,
             cartoonFill = 'black', cartoonOrder = 'under',
             zlim = c(min(raster::values(pc.pred.image.raster)),max(raster::values(pc.pred.image.raster))),
             normalized = normalized)
  }
  else{
    plotHeat(pc.pred.image.raster, IDlist, plotCartoon = TRUE, refShape = refShape, outline = outline, lines = lines,
             adjustCoords = TRUE, landList = landList, imageList = imageList, cartoonID = cartoonID, colpalette = colpalette,
             cartoonFill = 'black', cartoonOrder = 'under', zlim = zlim, normalized = normalized)
  }
}





