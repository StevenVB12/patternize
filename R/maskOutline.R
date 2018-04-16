#' Intersects a RasterStack with an outline. Everything outside of the outline will be removed
#' from the raster.
#'
#' @param RasterStack RasterStack to be masked.
#' @param outline xy coordinates that define outline.
#' @param refShape This can be 'target' in case the reference shape is a single sample (for
#'    registration analysis) or 'mean' if the images were transformed to a mean shape (only
#'    for meanshape when using landmark transformation)
#' @param landList Landmark list to be given when type = 'mean'.
#' @param adjustCoords Adjust landmark coordinates in case they are reversed compared to
#'    pixel coordinates (default = FALSE).
#' @param cartoonID ID of the sample for which the cartoon was drawn. Only has to be given when
#'    refShape is 'mean'.
#' @param IDlist List of sample IDs should be specified when refShape is 'mean'.
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to
#'    crop the original image used in landmark or registration analysis.
#' @param flipRaster Whether to flip raster along xy axis (in case there is an inconsistency
#'    between raster and outline coordinates).
#' @param flipOutline Whether to flip plot along x, y or xy axis.
#' @param imageList List of image as obtained from \code{\link[patternize]{makeList}} should
#'    be given if one wants to flip the outline or adjust landmark coordinates.
#' @param maskColor Color the masked area gets. Set to 0 for black (default) or 255 for white.
#'
#' @examples
#'
#' \dontrun{
#' data(imageList)
#' outline_BC0077 <- read.table(paste(system.file("extdata",  package = 'patternize'),
#' '/BC0077_outline.txt', sep=''), header = FALSE)
#'
#' masked <- maskOutline(imageList[[1]], outline_BC0077, refShape = 'target', flipOutline = 'y')
#' }
#'
#' @export
#' @import raster
#' @importFrom utils capture.output

maskOutline <-function(RasterStack,
                       outline,
                       refShape,
                       landList = NULL,
                       adjustCoords = FALSE,
                       cartoonID = NULL,
                       IDlist = NULL,
                       crop = c(0,0,0,0),
                       flipRaster = NULL,
                       flipOutline = NULL,
                       imageList = NULL,
                       maskColor = 0){

  if(is.list(imageList)){

    imageEx <- raster::extent(imageList[[1]])
  }
  else{
    imageEx <- raster::extent(imageList)
  }

  if(!is.null(flipOutline) || !is.null(flipRaster)){

    if(refShape != 'mean'){

      outline[,2] <- outline[,2] - crop[3]

    }
  }

  if(refShape != 'mean'){


    if(!is.null(flipOutline) && flipOutline == 'y' || !is.null(flipOutline) && flipOutline == 'xy'){

      outline[,2] <- outline[,2] + crop[3]

    }

    if(is.null(flipOutline) && !is.null(flipRaster) || !is.null(flipOutline) && flipOutline == 'x'){

      outline[,2] <- outline[,2] + ((crop[3] - imageEx[3]) - (imageEx[4] - crop[4])) + crop[3]

    }

    if(!is.null(flipOutline)){

      if(flipOutline == 'x'){

        outline[,1] <- imageEx[2] - outline[,1] + (crop[1] - imageEx[1]) - (imageEx[2] - crop[2])

      }

      if(flipOutline == 'y'){

        outline[,2] <- imageEx[4] - outline[,2]

      }

      if(flipOutline == 'xy'){

        outline[,1] <- imageEx[2] - outline[,1] + (crop[1] - imageEx[1]) - (imageEx[2] - crop[2])
        outline[,2] <- imageEx[4] - outline[,2]

      }
    }
  }

  if(refShape == 'mean'){

    indx <- which(IDlist == cartoonID)
    invisible(capture.output(landArray <- lanArray(landList, adjustCoords, imageList)))

    if(adjustCoords){

      extPicture <- raster::extent(imageList[[indx]])
      outline[,2] <- extPicture[4]-outline[,2]
    }

    invisible(capture.output(transformed <- Morpho::procSym(landArray)))


    invisible(capture.output(cartoonLandTrans <- Morpho::computeTransform(transformed$mshape,
                                                                          as.matrix(landArray[,,indx]), type="tps")))

    if(!is.null(flipOutline)){

      if(flipOutline == 'x'){
        outline[,1] = imageEx[2] - outline[,1] + imageEx[1]


      }

      if(flipOutline == 'y'){
        outline[,2] = imageEx[4] - outline[,2] + imageEx[3]

      }

      if(flipOutline == 'xy'){
        outline[,1] = imageEx[2] - outline[,1] + imageEx[1]
        outline[,2] = imageEx[4] - outline[,2] + imageEx[3]

      }
    }

    outline <- Morpho::applyTransform(as.matrix(outline), cartoonLandTrans)


  }


  poly <- sp::Polygons(list(sp::Polygon(outline)),paste("r"))

  polyList  <- c(poly)
  polyNames <- c(paste("r"))
  sr <- sp::SpatialPolygons(polyList)
  srdf <- sp::SpatialPolygonsDataFrame(sr, data.frame(1:length(polyNames), row.names=polyNames))

  imageExr <- raster::extent(RasterStack)
  r <- raster::raster(imageExr, nrow=dim(RasterStack)[1], ncol=dim(RasterStack)[2])
  rr <- raster::rasterize(srdf, r)

  if(!is.null(flipRaster)){
    if(flipRaster == 'x'){
      RasterStack <- raster::flip(RasterStack,'x')
    }
    if(flipRaster == 'y'){
      RasterStack <- raster::flip(RasterStack,'y')
    }
    if(flipRaster == 'xy'){
      RasterStack <- raster::flip(RasterStack,'x')
      RasterStack <- raster::flip(RasterStack,'y')
    }
  }

  RasterStack <- raster::mask(RasterStack, rr)
  RasterStack[is.na(RasterStack)] <- maskColor

  return(RasterStack)
}






