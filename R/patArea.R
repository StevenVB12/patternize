#' Calculate relative area of color pattern
#'
#' @param rList List of raster objects.
#' @param IDlist List of sample IDs.
#' @param refShape This can be 'target' in case the reference shape is a single sample (for registration analysis) or 'mean' if the images were transformed to a mean shape (only for meanshape when using landmark transformation)
#' @param outline xy coordinates that define outline.
#' @param landList Landmark landmarkList.
#' @param adjustCoords Adjust landmark coordinates.
#' @param cartoonID ID of the sample for which the cartoon was drawn.
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the original image used in landmark or registration analysis.
#' @param flipRaster Whether to flip raster along xy axis (in case there is an inconsistency between raster and outline coordinates).
#' @param flipOutline Whether to flip plot along x, y or xy axis.
#' @param imageList List of image should be given if one wants to flip the outline or adjust landmark coordinates.
#'
#' @examples
#' data(rasterList_lanRGB)
#' data(rasterList_regRGB)
#' data(rasterList_lanK)
#' data(rasterList_regK)
#'
#' data(imageList)
#'
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' outline_BC0077 <- read.table(paste(system.file("extdata",  package = 'patternize'), '/BC0077_outline.txt', sep=''), h= F)
#' prepath <- system.file("extdata", package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#' landmarkList <- makeList(IDlist, 'landmark', prepath, extension)
#' summedRaster_lanRGB <- sumRaster(rasterList_lanRGB, IDlist, type = 'RGB')
#'
#' patArea(rasterList_lanRGB, IDlist, refShape = 'mean', type = 'RGB', outline = outline_BC0077, landList = landmarkList, adjustCoords = TRUE, imageList = imageList, cartoonID = 'BC0077')
#' patArea(rasterList_regRGB, IDlist, refShape = 'target', type = 'RGB', outline = outline_BC0077,  crop = c(1000,4000,400,2500), adjustCoords = TRUE, imageList = imageList, cartoonID = 'BC0077', flipRaster = 'xy')
#' patArea(rasterList_lanK, IDlist, refShape = 'mean', type = 'k', outline = outline_BC0077, landList = landmarkList, adjustCoords = TRUE, imageList = imageList, cartoonID = 'BC0077')
#' patArea(rasterList_regK, IDlist, refShape = 'target', type = 'k', outline = outline_BC0077,  crop = c(1000,4000,400,2500), adjustCoords = TRUE, imageList = imageList, cartoonID = 'BC0077', flipRaster = 'xy')
#'
#' @export


patArea <-function(rList, IDlist, refShape, type, outline = NULL, landList = NULL, adjustCoords = FALSE, cartoonID = NULL, crop = NULL, flipRaster = NULL, flipOutline = NULL, imageList = NULL){



  if(!is.null(flipOutline) || !is.null(flipRaster)){

    imageEx <- raster::extent(imageList[[1]])

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

      extPicture <- extent(imageList[[indx]])
      outline[,2] <- extPicture[4]-outline[,2]
    }

    invisible(capture.output(transformed <- Morpho::procSym(landArray)))


    invisible(capture.output(cartoonLandTrans <- Morpho::computeTransform(transformed$mshape, as.matrix(landArray[,,indx]), type="tps")))
    outlineTrans <- Morpho::applyTransform(as.matrix(outline), cartoonLandTrans)


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

      }
    }
  }


  if(refShape == 'mean'){

    rasterEx <- raster::extent(min(outlineTrans[,1]),max(outlineTrans[,1]),min(outlineTrans[,2]),max(outlineTrans[,2]))
    rRe <- raster::raster(nrow=150,ncol=150)
    extent(rRe) <- rasterEx

    if(type == 'RGB'){
      newRaster <- raster::resample(rList[[IDlist[[1]]]], rRe)
    }
    if(type == 'k'){
      newRaster <- raster::resample(rList[[IDlist[[1]]]][[1]], rRe)
    }
    poly <- sp::Polygons(list(Polygon(outlineTrans)),paste("r"))

  }



  if(refShape == 'target'){

    rasterEx <- raster::extent(min(outline[,1]),max(outline[,1]),min(outline[,2]),max(outline[,2]))
    rRe <- raster::raster(nrow=150,ncol=150)
    extent(rRe) <- rasterEx

    if(type == 'RGB'){
      newRaster <- raster::resample(rList[[IDlist[[1]]]], rRe)
    }
    if(type == 'k'){
      newRaster <- raster::resample(rList[[IDlist[[1]]]][[1]], rRe)
    }

    poly <- sp::Polygons(list(Polygon(outline)),paste("r"))

  }

  polyList  <- c(poly)
  polyNames <- c(paste("r"))
  sr=sp::SpatialPolygons(polyList)
  srdf=sp::SpatialPolygonsDataFrame(sr, data.frame(1:length(polyNames), row.names=polyNames))

  r <- raster::raster(rasterEx, nrow=dim(newRaster)[1], ncol=dim(newRaster)[2])
  rr <-raster::rasterize(srdf, r)

  nrCellsOutline  <- raster::freq(rr, value=1)


  outDf <-c()

  if(type == 'RGB'){

    for(n in 1:length(IDlist)){

      rast <- rList[[IDlist[[n]]]]

      if(!is.null(flipRaster)){
        if(flipRaster == 'x'){
          rast <- raster::flip(rast,'x')
        }
        if(flipRaster == 'y'){
          rast <- raster::flip(rast,'y')
        }
        if(flipRaster == 'xy'){
          rast <- raster::flip(rast,'x')
          rast <- raster::flip(rast,'y')
        }
      }


      rast[rast == 0] <- NA
      rast <- raster::resample(rast,rr, method='ngb')

      nrCells  <- raster::freq(rast, value=1)

      outDf <- rbind(outDf, c(IDlist[[n]], nrCells/nrCellsOutline))
    }

    outDf <- as.data.frame(outDf)
    colnames(outDf) <- c('SampleId','Area')
    outDf$Area <- as.numeric(as.character(outDf$Area))

    return(outDf)
  }

  if(type == 'k'){

    kList <- list()

    for(n in 1:length(IDlist)){

      sRast <- rList[[IDlist[[n]]]]

      for(e in 1:length(sRast)){

        if(!is.null(flipRaster)){
          if(flipRaster == 'x'){
            sRast[[e]] <- raster::flip(sRast[[e]],'x')
          }
          if(flipRaster == 'y'){
            sRast[[e]] <- raster::flip(sRast[[e]],'y')
          }
          if(flipRaster == 'xy'){
            sRast[[e]] <- raster::flip(sRast[[e]],'x')
            sRast[[e]] <- raster::flip(sRast[[e]],'y')
          }
        }

        sRast[[e]][sRast[[e]] == 0] <- NA
        sRast[[e]] <- raster::resample(sRast[[e]], rr, method='ngb')

        nrCells  <- raster::freq(sRast[[e]], value=1)

        if(n == 1){
          kList[[e]] <- c(IDlist[[n]], nrCells/nrCellsOutline)
        }
        else{
          kList[[e]] <- rbind(kList[[e]], c(IDlist[[n]], nrCells/nrCellsOutline))
        }
      }
    }

    for(i in 1:length(kList)){
      kList[[i]] <- as.data.frame(kList[[i]])
      colnames(kList[[i]]) <- c('SampleId','Area')
      kList[[i]]$Area <- as.numeric(as.character(kList[[i]]$Area))
    }
    return(kList)
  }
}






