#' Align images using landmarks
#'
#' @param imageList List of RasterStack objects.
#' @param landList Landmark list as returned by \code{\link[patternize]{makeList}}.
#' @param IDlist List of sample IDs should be specified when masking outline and transformRef
#'    is 'meanshape'.
#' @param adjustCoords Adjust landmark coordinates in case they are reversed compared to pixel
#'    coordinates (default = FALSE).
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#' @param res Resolution vector c(x,y) for output rasters (default = c(300,300)). This should be
#'    reduced if the number of pixels in the image is lower than th raster.
#' @param transformRef ID or landmark matrix of reference sample for shape to which color patterns
#'    will be transformed to. Can be 'meanshape' for transforming to mean shape of Procrustes
#'    analysis.
#' @param transformType Transformation type as used by \code{\link[Morpho]{computeTransform}}
#'    (default ='tps').
#' @param maskOutline When outline is specified, everything outside of the outline will be masked for
#'    the color extraction (default = NULL). This can be a list of multiple outlines.
#' @param removebg Integer or RGB vector indicating the range of RGB threshold to remove from
#'    image (e.g. 100 removes pixels with average RGB > 100; default = NULL).
#' @param removebgColOffset Color offset for color background extraction (default = 0.10).
#' @param inverse If TRUE, areas withing the outline will be masked. If maskOutline is a list, this
#'    should also be a list.
#' @param cartoonID ID of the sample for which the cartoon was drawn and will be used for masking
#'    (should be set when transformRef = 'meanShape').
#' @param refImage Image (RasterStack) used for target. Use raster::stack('filename').
#' @param plotTransformed Plot transformed image (default = FALSE).
#' @param format ImageJ (Fiji) or tps format (default = 'imageJ').
#'
#' @return List of aligned RasterStack objects.
#'
#' @export
#' @import raster
#' @importFrom Morpho procSym computeTransform applyTransform

alignLan <- function(imageList,
                     landList,
                     IDlist = NULL,
                     adjustCoords = FALSE,
                     resampleFactor = NULL,
                     res = c(300,300),
                     transformRef = 'meanshape',
                     transformType = 'tps',
                     maskOutline = NULL,
                     removebg = NULL,
                     removebgColOffset = 0.1,
                     inverse = FALSE,
                     cartoonID = NULL,
                     refImage = NULL,
                     plotTransformed = FALSE,
                     format = 'imageJ'){

  rasterList <- list()


  if(!is.null(cartoonID)){
    if(cartoonID %in% names(imageList)){
      imageEx <- raster::extent(imageList[[cartoonID]])
    }
  }
  else if(!is.null(refImage)){
    imageEx <- raster::extent(refImage)
  }
  else{
    stop("if cartoonID not in image list, please provide image to refImage argument.")
  }

  if(!is.data.frame(maskOutline)){
    maskOutlineList <- maskOutline
    maskOutline <- maskOutline[[1]]
    inverseList <- inverse
    inverse <- inverse[[1]]

    if(format == 'tps'){
      maskOutline[,2] <- (imageEx[4]-maskOutline[,2])

      for(e in 2:length(maskOutlineList)){
        maskOutlineList[[e]][,2] <- (imageEx[4]-maskOutlineList[[e]][,2])
      }
    }
  }

  else{
    maskOutline <- as.matrix(maskOutline)
    if(format == 'tps'){
      maskOutline[,2] <- (imageEx[4]-maskOutline[,2])
    }
  }

  # Check whether imageList and landList have the same length
  if(length(imageList) != length(landList)){
    stop("imageList is not of the same length as lanArray")
  }

  for(n in 1:length(imageList)){
    if(names(imageList)[n] != names(landList)[n]){
      stop("samples are not in the same order in imageList and lanArray")
    }
  }

  # Make landmark array
  lanArray <- lanArray(landList, adjustCoords, imageList)


  # Set the reference shape
  if(is.matrix(transformRef)){

    refShape <- transformRef
  }

  if(!is.matrix(transformRef)){

    if(transformRef == 'meanshape'){

      invisible(capture.output(transformed <- Morpho::procSym(lanArray)))
      refShape <- transformed$mshape
    }

    if(transformRef %in% names(landList)){

      e <- which(names(landList) == transformRef)
      refShape <- lanArray[,,e]
    }
  }

  # Transform the outline for masking if 'meanShape'
  if(!is.null(cartoonID)){
    if(any(c(!is.null(maskOutline), transformRef == 'meanshape'))){

    maskOutlineNew <- maskOutline
    maskOutlineNew[,2] <- imageEx[4]-maskOutlineNew[,2]
    }
  }

  if(is.null(cartoonID)){
    maskOutlineNew <- maskOutline
  }

  if(all(c(!is.null(maskOutline), transformRef == 'meanshape'))){

    if(is.character(transformRef)){
      indx <- which(names(landList) == cartoonID)
      transRefLan <- as.matrix(lanArray[,,indx])
    }
    else{
      transRefLan <- transformRef
    }
    invisible(capture.output(cartoonLandTrans <- Morpho::computeTransform(refShape,
                                                                          as.matrix(transRefLan),
                                                                          type='tps')))

    maskOutlineMean <- Morpho::applyTransform(as.matrix(maskOutlineNew), cartoonLandTrans)
  }

  if(all(c(!is.null(maskOutline), transformRef[1] != 'meanshape', !is.null(cartoonID)))){

    maskOutlineNew[,2] <- imageEx[4] - maskOutlineNew[,2]
  }

  # Function used in focal to calculate average for NA values in tranformed matrix
  fill.na <- function(x, i=5) {
    if( is.na(x)[i] ) {
      return( round(mean(x, na.rm=TRUE),0) )
    } else {
      return( round(x[i],0) )
    }
  }

  # Run the loop for each sample
  for(n in 1:length(imageList)){

    image <- imageList[[n]]
    extRaster <- raster::extent(image)

    # Reduce resolution
    if(!is.null(resampleFactor)){
      image <- redRes(image, resampleFactor)
    }

    if(is.vector(removebg)){

      if(length(removebg) > 1){
        toMask <- apply(raster::as.array(image), 1:2, function(x) all(abs(x-removebg) < removebgColOffset*255))
      }
      else{
        toMask <- apply(raster::as.array(image), 1:2, function(x) all(x > removebg))
      }

      toMaskR <- raster::raster(as.matrix(toMask))
      raster::extent(toMaskR) <- raster::extent(image)
      toMaskR[toMaskR == 0] <- NA

      image <- raster::mask(image, toMaskR, inverse = TRUE)
      # image[is.na(image)] <- 0
    }

    # Transform image using landmarks
    invisible(capture.output(transMatrix <- Morpho::computeTransform(refShape,
                                                                     as.matrix(lanArray[,,n]),
                                                                     type = 'tps')))

    imageDF1 <- raster::as.data.frame(image[[1]], xy = TRUE)
    imageDF2 <- raster::as.data.frame(image[[2]], xy = TRUE)
    imageDF3 <- raster::as.data.frame(image[[3]], xy = TRUE)

    invisible(capture.output(imageT <- Morpho::applyTransform(as.matrix(imageDF1)[,1:2], transMatrix)))

    r <- raster::raster(nrow = dim(image)[1], ncol = dim(image)[2])


    raster::extent(r) <- imageEx #c(min(imageT[,1]),max(imageT[,1]),min(imageT[,2]),max(imageT[,2]))


    # Rasterize the transformed image and fill in NA values using
    imageT1r <- raster::rasterize(imageT, field = imageDF1[,3], r, fun = mean)
    imageT1rf <- focal(imageT1r, w=matrix(1,nrow=3,ncol=3), fun=fill.na, pad = TRUE, na.rm = FALSE)

    imageT2r <- raster::rasterize(imageT, field = imageDF2[,3], r, fun = mean)
    imageT2rf <- focal(imageT2r, w=matrix(1,nrow=3,ncol=3), fun=fill.na, pad = TRUE, na.rm = FALSE)

    imageT3r <- raster::rasterize(imageT, field = imageDF3[,3], r, fun = mean)
    imageT3rf <- focal(imageT3r, w=matrix(1,nrow=3,ncol=3), fun=fill.na, pad = TRUE, na.rm = FALSE)

    imageTr <- raster::stack(imageT1rf, imageT2rf, imageT3rf)

    # imageTr[is.na(imageTr)] <- 255

    # MaskOutline
    if(!is.null(maskOutline)){

      if(transformRef[1] != 'meanshape'){

        imageTr <- maskOutline(imageTr, maskOutlineNew, refShape = 'target', crop = c(0,0,0,0),
                               maskColor = 0, imageList = imageList, cartoonID = cartoonID, inverse = inverse)

        cropEx <- c(min(maskOutlineNew[,1]), max(maskOutlineNew[,1]), min(maskOutlineNew[,2]), max(maskOutlineNew[,2]))

        if(exists('maskOutlineList')){
          for(e in 2:length(maskOutlineList)){
            imageTr <- maskOutline(imageTr, maskOutlineList[[e]], refShape = 'target', crop = c(0,0,0,0),
                                   maskColor = 0, imageList = imageList, cartoonID = cartoonID, inverse = inverseList[[e]])
          }
        }

      }
      if(transformRef[1] == 'meanshape'){

        imageTr <- maskOutline(imageTr, outline = maskOutline, refShape = 'mean',
                               IDlist = IDlist, landList = landList, imageList = imageList,
                               adjustCoords = TRUE, cartoonID = cartoonID)

        cropEx <- c(min(maskOutlineMean[,1]), max(maskOutlineMean[,1]), min(maskOutlineMean[,2]), max(maskOutlineMean[,2]))
      }

      imageTr <- raster::crop(imageTr, cropEx)
      raster::extent(imageTr) <- cropEx
    }

    # Plot transformed raster
    if(plotTransformed){

      # imageTr <- raster::flip(imageTr, 'y')
      imageTr2 <- imageTr
      imageTr2[is.na(imageTr2)] <- 0

      x <- as.array(imageTr2)/255
      cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
      uniqueCols <- unique(cols)
      x2 <- match(cols, uniqueCols)
      dim(x2) <- dim(x)[1:2]
      raster::image(t(x2), col=uniqueCols, yaxt='n', xaxt='n', main = paste(names(landList)[n],'transformed', sep=' '))

      # imageTr <- raster::flip(imageTr, 'y')
    }

    # set resolution
    inCols <- res[1]
    inRows <- res[2]

    resampledRaster <- raster::raster(ncol=inCols, nrow=inRows)

    raster::extent(resampledRaster) <- raster::extent(imageTr)

    resampled <- raster::resample(imageTr, resampledRaster, method = 'ngb')

    rasterList[[names(landList)[n]]] <- resampled

    print(paste('sample', names(landList)[n], 'processed', sep=' '))}

  return(rasterList)
}
