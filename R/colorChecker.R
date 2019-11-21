#' Calibrate images using ColorChecker.
#'
#' @param IDlist List of sample IDs.
#' @param prepath Prepath (default = NULL).
#' @param extension Extension (default = NULL).
#' @param colorCheckerType Type of colorChecker. Options are 'X-Rite ' and 'ColorGauge Micro
#' Analyzer' (default = 'X-Rite ').
#' @param fixedCorners Specify whether to set the coordinates of the colorChecker corners
#' for every image (default = FALSE).
#' @param patchSize Proportion of ColorChecker patch that will be used for observed RGB values
#' (default = 0.6).
#' @param colorCheckerXY Landmark list of colorChecker corners as returned
#'   by \code{\link[patternize]{makeList}}. The image will not be plotted.
#' @param fixedModel Precalculated model to adjust colors. Should be a listof a model for R, G
#'   and B (the colorChecker function gives as output such a list obtained from the last image
#'   in the analysis).
#' @param resampleFactor Integer for downsampling used by \code{\link{redRes}}.
#'
#' @return  Calibrated image(s) ('filename_calibrated.jpg')
#'
#' @export
#' @import raster
#' @importFrom stats lm
#' @importFrom graphics locator
#' @importFrom imager load.image save.image as.cimg width R G B resize width height
#' @importFrom sp Polygons SpatialPolygons SpatialPolygonsDataFrame

colorChecker <- function(IDlist,
                         prepath = NULL,
                         extension = NULL,
                         colorCheckerType = 'X-Rite',
                         fixedCorners = FALSE,
                         patchSize = 0.6,
                         colorCheckerXY = NULL,
                         fixedModel = NULL,
                         resampleFactor = NULL){

  prop <- 1- patchSize

  for(n in 1:length(IDlist)){

    print(paste('Processing sample', n,  IDlist[n], sep=' '))

    if(is.null(prepath)){
      im <- imager::load.image(paste(IDlist[n], extension, sep=''))
    }
    else{
      im <- imager::load.image(paste(prepath,'/',IDlist[n], extension, sep=''))
    }

    # Reduce resolution
    if(!is.null(resampleFactor)){

      imRed <- imager::resize(im,round(width(im)/resampleFactor),round(height(im)/resampleFactor))

    }
    else{
      imRed <- im
    }

    mR <- raster::as.matrix(imager::R(imRed))*255
    mG <- raster::as.matrix(imager::G(imRed))*255
    mB <- raster::as.matrix(imager::B(imRed))*255

    # mR_O <- raster::as.matrix(imager::R(im))*255
    # mG_O <- raster::as.matrix(imager::G(im))*255
    # mB_O <- raster::as.matrix(imager::B(im))*255

    rR <- raster::raster(mR)
    rG <- raster::raster(mG)
    rB <- raster::raster(mB)

    extent(rR) <- c(0, dim(imRed)[2], 0, dim(imRed)[1])
    extent(rG) <- c(0, dim(imRed)[2], 0, dim(imRed)[1])
    extent(rB) <- c(0, dim(imRed)[2], 0, dim(imRed)[1])

    rR <- flip(t(rR),'y')
    rG <- flip(t(rG),'y')
    rB <- flip(t(rB),'y')

    # X-Rite colorChecker
    if(colorCheckerType == 'X-Rite' & is.null(fixedModel)){
      if(is.null(colorCheckerXY)){
        print("Click on the outer corners of the X-Rite ColorChecker (1: ~brown, 2: ~cyan, 3: ~black, 4: ~white)")

        layout(matrix(c(1,1), 2, 1, byrow = TRUE))
        plot(imRed, xlim = c(0, width(imRed)))

        if(fixedCorners == FALSE){
          xy <- locator(n=4)
        }

        if(fixedCorners == TRUE & n ==1 ){
          xy <- locator(n=4)
        }
      }
      if(!is.null(colorCheckerXY)){

        xy <- list()
        xy$x <- colorCheckerXY[[n]][,1]
        xy$y <- colorCheckerXY[[n]][,2]

      }


      print('Matching ColorCecker patches...')

      if(is.null(colorCheckerXY)){
        polygon(xy$x, xy$y, border = "green")
      }

      xyDF <- as.data.frame(xy)

      line1 <- xyDF[1:2,]
      line2 <- xyDF[2:3,]
      line3 <- xyDF[3:4,]
      line4 <- xyDF[c(4,1),]

      xdiff <- (line1$x[2]-line1$x[1])/6
      ydiff <- (line1$y[2]-line1$y[1])/6

      xmin <- line1$x[1]
      ymin <- line1$y[1]

      xySubA <- list(x = c((xmin+xmin+xdiff)/2,(xmin+xdiff+xmin+xdiff*2)/2,(xmin+xdiff*2+xmin+xdiff*3)/2,(xmin+xdiff*3+xmin+xdiff*4)/2,(xmin+xdiff*4+xmin+xdiff*5)/2,(xmin+xdiff*5+xmin+xdiff*6)/2),
                     y = c((ymin+ymin+ydiff)/2,(ymin+ydiff+ymin+ydiff*2)/2,(ymin+ydiff*2+ymin+ydiff*3)/2,(ymin+ydiff*3+ymin+ydiff*4)/2,(ymin+ydiff*4+ymin+ydiff*5)/2,(ymin+ydiff*5+ymin+ydiff*6)/2))

      xySubDF_1A <- as.data.frame(xySubA)

      xySubBa <- list(x = c(xmin+xdiff-(xdiff/2)*prop,xmin+xdiff*2-(xdiff/2)*prop,xmin+xdiff*3-(xdiff/2)*prop,xmin+xdiff*4-(xdiff/2)*prop,xmin+xdiff*5-(xdiff/2)*prop,xmin+xdiff*6-(xdiff/2)*prop),
                      y = c(ymin+ydiff-(ydiff/2)*prop,ymin+ydiff*2-(ydiff/2)*prop,ymin+ydiff*3-(ydiff/2)*prop,ymin+ydiff*4-(ydiff/2)*prop,ymin+ydiff*5-(ydiff/2)*prop,ymin+ydiff*6-(ydiff/2)*prop))
      xySubBb <- list(x = c(xmin+(xdiff/2)*prop,xmin+xdiff+(xdiff/2)*prop,xmin+xdiff*2+(xdiff/2)*prop,xmin+xdiff*3+(xdiff/2)*prop,xmin+xdiff*4+(xdiff/2)*prop,xmin+xdiff*5+(xdiff/2)*prop),
                      y = c(ymin+(ydiff/2)*prop,ymin+ydiff+(ydiff/2)*prop,ymin+ydiff*2+(ydiff/2)*prop,ymin+ydiff*3+(ydiff/2)*prop,ymin+ydiff*4+(ydiff/2)*prop,ymin+ydiff*5+(ydiff/2)*prop))

      xySubDF_1Ba <- as.data.frame(xySubBa)
      xySubDF_1Bb <- as.data.frame(xySubBb)


      xdiff <- (line3$x[2]-line3$x[1])/6
      ydiff <- (line3$y[2]-line3$y[1])/6

      xmin <- line3$x[1]
      ymin <- line3$y[1]

      xySubA <- list(x = c((xmin+xmin+xdiff)/2,(xmin+xdiff+xmin+xdiff*2)/2,(xmin+xdiff*2+xmin+xdiff*3)/2,(xmin+xdiff*3+xmin+xdiff*4)/2,(xmin+xdiff*4+xmin+xdiff*5)/2,(xmin+xdiff*5+xmin+xdiff*6)/2),y = c((ymin+ymin+ydiff)/2,(ymin+ydiff+ymin+ydiff*2)/2,(ymin+ydiff*2+ymin+ydiff*3)/2,(ymin+ydiff*3+ymin+ydiff*4)/2,(ymin+ydiff*4+ymin+ydiff*5)/2,(ymin+ydiff*5+ymin+ydiff*6)/2))

      xySubDF_3A <- as.data.frame(xySubA)

      xySubBa <- list(x = c(xmin+xdiff-(xdiff/2)*prop,xmin+xdiff*2-(xdiff/2)*prop,xmin+xdiff*3-(xdiff/2)*prop,xmin+xdiff*4-(xdiff/2)*prop,xmin+xdiff*5-(xdiff/2)*prop,xmin+xdiff*6-(xdiff/2)*prop),
                      y = c(ymin+ydiff-(ydiff/2)*prop,ymin+ydiff*2-(ydiff/2)*prop,ymin+ydiff*3-(ydiff/2)*prop,ymin+ydiff*4-(ydiff/2)*prop,ymin+ydiff*5-(ydiff/2)*prop,ymin+ydiff*6-(ydiff/2)*prop))
      xySubBb <- list(x = c(xmin+(xdiff/2)*prop,xmin+xdiff+(xdiff/2)*prop,xmin+xdiff*2+(xdiff/2)*prop,xmin+xdiff*3+(xdiff/2)*prop,xmin+xdiff*4+(xdiff/2)*prop,xmin+xdiff*5+(xdiff/2)*prop),
                      y = c(ymin+(ydiff/2)*prop,ymin+ydiff+(ydiff/2)*prop,ymin+ydiff*2+(ydiff/2)*prop,ymin+ydiff*3+(ydiff/2)*prop,ymin+ydiff*4+(ydiff/2)*prop,ymin+ydiff*5+(ydiff/2)*prop))

      xySubDF_3Ba <- as.data.frame(xySubBa)
      xySubDF_3Bb <- as.data.frame(xySubBb)

      xdiff <- (line2$x[2]-line2$x[1])/4
      ydiff <- (line2$y[2]-line2$y[1])/4

      xmin <- line2$x[1]
      ymin <- line2$y[1]

      xySubA <- list(x = c((xmin+xmin+xdiff)/2,(xmin+xdiff+xmin+xdiff*2)/2,(xmin+xdiff*2+xmin+xdiff*3)/2,(xmin+xdiff*3+xmin+xdiff*4)/2),y = c((ymin+ymin+ydiff)/2,(ymin+ydiff+ymin+ydiff*2)/2,(ymin+ydiff*2+ymin+ydiff*3)/2,(ymin+ydiff*3+ymin+ydiff*4)/2))

      xySubDF_2A <- as.data.frame(xySubA)

      xySubBa <- list(x = c(xmin+xdiff-(xdiff/2)*prop,xmin+xdiff*2-(xdiff/2)*prop,xmin+xdiff*3-(xdiff/2)*prop,xmin+xdiff*4-(xdiff/2)*prop),
                      y = c(ymin+ydiff-(ydiff/2)*prop,ymin+ydiff*2-(ydiff/2)*prop,ymin+ydiff*3-(ydiff/2)*prop,ymin+ydiff*4-(ydiff/2)*prop))
      xySubBb <- list(x = c(xmin+(xdiff/2)*prop,xmin+xdiff+(xdiff/2)*prop,xmin+xdiff*2+(xdiff/2)*prop,xmin+xdiff*3+(xdiff/2)*prop),
                      y = c(ymin+(ydiff/2)*prop,ymin+ydiff+(ydiff/2)*prop,ymin+ydiff*2+(ydiff/2)*prop,ymin+ydiff*3+(ydiff/2)*prop))

      xySubDF_2Ba <- as.data.frame(xySubBa)
      xySubDF_2Bb <- as.data.frame(xySubBb)

      xdiff <- (line4$x[2]-line4$x[1])/4
      ydiff <- (line4$y[2]-line4$y[1])/4

      xmin <- line4$x[1]
      ymin <- line4$y[1]

      xySubA <- list(x = c((xmin+xmin+xdiff)/2,(xmin+xdiff+xmin+xdiff*2)/2,(xmin+xdiff*2+xmin+xdiff*3)/2,(xmin+xdiff*3+xmin+xdiff*4)/2),
                     y = c((ymin+ymin+ydiff)/2,(ymin+ydiff+ymin+ydiff*2)/2,(ymin+ydiff*2+ymin+ydiff*3)/2,(ymin+ydiff*3+ymin+ydiff*4)/2))

      xySubDF_4A <- as.data.frame(xySubA)

      xySubBa <- list(x = c(xmin+xdiff-(xdiff/2)*prop,xmin+xdiff*2-(xdiff/2)*prop,xmin+xdiff*3-(xdiff/2)*prop,xmin+xdiff*4-(xdiff/2)*prop),
                      y = c(ymin+ydiff-(ydiff/2)*prop,ymin+ydiff*2-(ydiff/2)*prop,ymin+ydiff*3-(ydiff/2)*prop,ymin+ydiff*4-(ydiff/2)*prop))
      xySubBb <- list(x = c(xmin+(xdiff/2)*prop,xmin+xdiff+(xdiff/2)*prop,xmin+xdiff*2+(xdiff/2)*prop,xmin+xdiff*3+(xdiff/2)*prop),
                      y = c(ymin+(ydiff/2)*prop,ymin+ydiff+(ydiff/2)*prop,ymin+ydiff*2+(ydiff/2)*prop,ymin+ydiff*3+(ydiff/2)*prop))

      xySubDF_4Ba <- as.data.frame(xySubBa)
      xySubDF_4Bb <- as.data.frame(xySubBb)

      if(is.null(colorCheckerXY)){
        points(xySubDF_1A, pch=20, col = 'green')
        points(xySubDF_1Ba, pch=20, col = 'red')
        points(xySubDF_1Bb, pch=20, col = 'red')

        points(xySubDF_2A, pch=20, col = 'green')
        points(xySubDF_2Ba, pch=20, col = 'red')
        points(xySubDF_2Bb, pch=20, col = 'red')

        points(xySubDF_3A, pch=20, col = 'green')
        points(xySubDF_3Ba, pch=20, col = 'red')
        points(xySubDF_3Bb, pch=20, col = 'red')

        points(xySubDF_4A, pch=20, col = 'green')
        points(xySubDF_4Ba, pch=20, col = 'red')
        points(xySubDF_4Bb, pch=20, col = 'red')
      }

      labels <- list(c(1,7,13,19),
                     c(2,8,14,20),
                     c(3,9,15,21),
                     c(4,10,16,22),
                     c(5,11,17,23),
                     c(6,12,18,24))


      xyTot <- c()
      for(e in 1:nrow(xySubDF_1A)){

        xyLine1 <- xySubDF_1A[e,]
        xyLine2 <- xySubDF_3A[6:1,][e,]

        xdiff <- (xyLine2$x-xyLine1$x)/4
        ydiff <- (xyLine2$y-xyLine1$y)/4

        xmin <- xyLine1$x
        ymin <- xyLine1$y

        xySub <- list(x = c((xmin+xmin+xdiff)/2,(xmin+xdiff+xmin+xdiff*2)/2,(xmin+xdiff*2+xmin+xdiff*3)/2,(xmin+xdiff*3+xmin+xdiff*4)/2),
                      y = c((ymin+ymin+ydiff)/2,(ymin+ydiff+ymin+ydiff*2)/2,(ymin+ydiff*2+ymin+ydiff*3)/2,(ymin+ydiff*3+ymin+ydiff*4)/2))

        xySubDF <- as.data.frame(xySub)

        if(is.null(colorCheckerXY)){
          text(xySubDF, label=labels[[e]], col = 'green')
        }

        xySubDFLabel <- cbind(xySubDF, label=labels[[e]])



        xyLine1a <- xySubDF_1Ba[e,]
        xyLine1b <- xySubDF_1Bb[e,]

        xyLine3a <- xySubDF_3Bb[6:1,][e,]
        xyLine3b <- xySubDF_3Ba[6:1,][e,]

        xdiffa <- (xyLine3a$x-xyLine1a$x)/4
        ydiffa <- (xyLine3a$y-xyLine1a$y)/4

        xdiffb <- (xyLine3b$x-xyLine1b$x)/4
        ydiffb <- (xyLine3b$y-xyLine1b$y)/4

        xmina <- xyLine1a$x
        ymina <- xyLine1a$y

        xminb <- xyLine1b$x
        yminb <- xyLine1b$y


        xySubAa <- list(x = c(xmina+xdiffa-(xdiffa/2)*prop,xmina+xdiffa*2-(xdiffa/2)*prop,xmina+xdiffa*3-(xdiffa/2)*prop,xmina+xdiffa*4-(xdiffa/2)*prop),
                        y = c(ymina+ydiffa-(ydiffa/2)*prop,ymina+ydiffa*2-(ydiffa/2)*prop,ymina+ydiffa*3-(ydiffa/2)*prop,ymina+ydiffa*4-(ydiffa/2)*prop))
        xySubAb <- list(x = c(xmina+(xdiffa/2)*prop,xmina+xdiffa+(xdiffa/2)*prop,xmina+xdiffa*2+(xdiffa/2)*prop,xmina+xdiffa*3+(xdiffa/2)*prop),
                        y = c(ymina+(ydiffa/2)*prop,ymina+ydiffa+(ydiffa/2)*prop,ymina+ydiffa*2+(ydiffa/2)*prop,ymina+ydiffa*3+(ydiffa/2)*prop))

        xySubBa <- list(x = c(xminb+xdiffb-(xdiffb/2)*prop,xminb+xdiffb*2-(xdiffb/2)*prop,xminb+xdiffb*3-(xdiffb/2)*prop,xminb+xdiffb*4-(xdiffb/2)*prop),
                        y = c(yminb+ydiffb-(ydiffb/2)*prop,yminb+ydiffb*2-(ydiffb/2)*prop,yminb+ydiffb*3-(ydiffb/2)*prop,yminb+ydiffb*4-(ydiffb/2)*prop))
        xySubBb <- list(x = c(xminb+(xdiffb/2)*prop,xminb+xdiffb+(xdiffb/2)*prop,xminb+xdiffb*2+(xdiffb/2)*prop,xminb+xdiffb*3+(xdiffb/2)*prop),
                        y = c(yminb+(ydiffb/2)*prop,yminb+ydiffb+(ydiffb/2)*prop,yminb+ydiffb*2+(ydiffb/2)*prop,yminb+ydiffb*3+(ydiffb/2)*prop))

        xySubAaDF <- as.data.frame(xySubAa)
        xySubAbDF <- as.data.frame(xySubAb)
        xySubBaDF <- as.data.frame(xySubBa)
        xySubBbDF <- as.data.frame(xySubBb)

        xySubRow <- cbind(xySubAaDF,xySubAbDF,xySubBaDF,xySubBbDF, label=labels[[e]])
        colnames(xySubRow) <- c('x1','y1','x2','y2','x3','y3','x4','y4','label')

        xyTot <- rbind(xyTot, xySubRow)

      }


      xyTot$imR <- NA
      xyTot$imG <- NA
      xyTot$imB <- NA

      xyTot <- xyTot[order(xyTot$label),]

      for(e in 1:nrow(xyTot)){

        print(paste('Calculating observed RGB values for patch', e, sep = ' '))
        if(is.null(colorCheckerXY)){
          polygon(c(xyTot$x1[e],xyTot$x2[e],xyTot$x4[e],xyTot$x3[e]),
                  c(xyTot$y1[e],xyTot$y2[e],xyTot$y4[e],xyTot$y3[e]), border = 'red')
        }

        outline <- rbind(c(xyTot$x1[e],xyTot$y1[e]),
                         c(xyTot$x2[e],xyTot$y2[e]),
                         c(xyTot$x4[e],xyTot$y4[e]),
                         c(xyTot$x3[e],xyTot$y3[e]))

        poly <- sp::Polygons(list(sp::Polygon(outline)),paste("r"))
        polyList  <- c(poly)
        polyNames <- c(paste("r"))
        sr <- sp::SpatialPolygons(polyList)
        srdf <- sp::SpatialPolygonsDataFrame(sr, data.frame(1:length(polyNames), row.names=polyNames))

        extrR <- raster::extract(rR, srdf)
        extrG <- raster::extract(rG, srdf)
        extrB <- raster::extract(rB, srdf)

        xyTot$imR[e] <- mean(extrR[[1]])
        xyTot$imG[e] <- mean(extrG[[1]])
        xyTot$imB[e] <- mean(extrB[[1]])
      }



      # Colorimetric values for ColorCheker targets
      l1 <- c(1, 115, 82, 68)
      l2 <- c(2, 194, 150, 130)
      l3 <- c(3, 98, 122, 157)
      l4 <- c(4, 87, 108, 67)
      l5 <- c(5, 133, 128, 177)
      l6 <- c(6, 103, 189, 170)
      l7 <- c(7, 214, 126, 44)
      l8 <- c(8, 80, 91, 166)
      l9 <- c(9, 193, 90, 99)
      l10 <- c(10, 94, 60, 108)
      l11 <- c(11, 157, 188, 64)
      l12 <- c(12, 224, 163, 46)
      l13 <- c(13, 56, 61, 150)
      l14 <- c(14, 70, 148, 73)
      l15 <- c(15, 175, 54, 60)
      l16 <- c(16, 231, 199, 31)
      l17 <- c(17, 187, 86, 149)
      l18 <- c(18, 8, 133, 161)
      l19 <- c(19, 243, 243, 242)
      l20 <- c(20, 200, 200, 200)
      l21 <- c(21, 160, 160, 160)
      l22 <- c(22, 122, 122, 121)
      l23 <- c(23, 85, 85, 85)
      l24 <- c(24, 52, 52, 52)

      ColorCheckerRGB <- as.data.frame(rbind(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24))
      colnames(ColorCheckerRGB) <- c('label', 'sR', 'sG', 'sB')

      dat <- merge(xyTot, ColorCheckerRGB, by = 'label')


      print('Calculating polynomial regression...')

      sR <- dat$sR
      sG <- dat$sG
      sB <- dat$sB

      imR <- dat$imR
      imG <- dat$imG
      imB <- dat$imB

      modelR <- lm(sR ~ imR + imG + imB + imR^2 + imG^2 + imB^2)
      modelG <- lm(sG ~ imR + imG + imB + imR^2 + imG^2 + imB^2)
      modelB <- lm(sB ~ imR + imG + imB + imR^2 + imG^2 + imB^2)


      dfIm = data.frame(
        imR = matrix(mR, ncol=1),
        imG = matrix(mG, ncol=1),
        imB = matrix(mB, ncol=1)
      )

      print('Calibrating colors...')
      prR <- predict(modelR, dfIm)
      prG <- predict(modelG, dfIm)
      prB <- predict(modelB, dfIm)

      dfCal <- as.data.frame(cbind(prR,prG,prB))

      print('Rebuilding image...')
      Ri = matrix(dfCal$prR, nrow=dim(imRed)[1])
      Gi = matrix(dfCal$prG, nrow=dim(imRed)[1])
      Bi = matrix(dfCal$prB, nrow=dim(imRed)[1])

      imCal = array(dim=dim(imRed))
      imCal[,,,1] = Ri
      imCal[,,,2] = Gi
      imCal[,,,3] = Bi

      imCal <- imager::as.cimg(imCal)

      if(is.null(colorCheckerXY)){
        layout(matrix(c(1,2), 2, 1, byrow = TRUE))

          plot(imRed)
          plot(imCal)




        if(fixedCorners == FALSE & is.null(fixedModel)){
          todo <- readline(prompt="Press [enter] to continue and save image >>> ")
        }
      }
    }


    # ColorGauge Micro Analyzer colorChecker
    if(colorCheckerType == 'ColorGauge Micro Analyzer' & is.null(fixedModel)){

      if(is.null(colorCheckerXY)){
        print("Click on the outer corners of the ColorGauge Micro Analyzer ColorChecker (1: ~brown, 2: ~cyan, 3: ~purple, 4: ~red)")

        layout(matrix(c(1,1), 2, 1, byrow = TRUE))
        plot(imRed, xlim = c(0, width(imRed)))

        if(fixedCorners == FALSE){
          xy <- locator(n=4)
        }

        if(fixedCorners == TRUE & n ==1){
          xy <- locator(n=4)
        }
      }

      if(!is.null(colorCheckerXY)){

        xy <- list()
        xy$x <- colorCheckerXY[[n]][,1]
        xy$y <- colorCheckerXY[[n]][,2]

      }

      print('Matching ColorCecker patches...')

      if(is.null(colorCheckerXY)){
        polygon(xy$x, xy$y, border = "green")
      }

      xyDF <- as.data.frame(xy)

      line1 <- xyDF[1:2,]
      line2 <- xyDF[2:3,]
      line3 <- xyDF[3:4,]
      line4 <- xyDF[c(4,1),]

      xdiff <- (line1$x[2]-line1$x[1])/6
      ydiff <- (line1$y[2]-line1$y[1])/6

      xmin <- line1$x[1]
      ymin <- line1$y[1]

      xySubA <- list(x = c((xmin+xmin+xdiff)/2,(xmin+xdiff+xmin+xdiff*2)/2,(xmin+xdiff*2+xmin+xdiff*3)/2,(xmin+xdiff*3+xmin+xdiff*4)/2,(xmin+xdiff*4+xmin+xdiff*5)/2,(xmin+xdiff*5+xmin+xdiff*6)/2),
                     y = c((ymin+ymin+ydiff)/2,(ymin+ydiff+ymin+ydiff*2)/2,(ymin+ydiff*2+ymin+ydiff*3)/2,(ymin+ydiff*3+ymin+ydiff*4)/2,(ymin+ydiff*4+ymin+ydiff*5)/2,(ymin+ydiff*5+ymin+ydiff*6)/2))

      xySubDF_1A <- as.data.frame(xySubA)

      xySubBa <- list(x = c(xmin+xdiff-(xdiff/2)*prop,xmin+xdiff*2-(xdiff/2)*prop,xmin+xdiff*3-(xdiff/2)*prop,xmin+xdiff*4-(xdiff/2)*prop,xmin+xdiff*5-(xdiff/2)*prop,xmin+xdiff*6-(xdiff/2)*prop),
                      y = c(ymin+ydiff-(ydiff/2)*prop,ymin+ydiff*2-(ydiff/2)*prop,ymin+ydiff*3-(ydiff/2)*prop,ymin+ydiff*4-(ydiff/2)*prop,ymin+ydiff*5-(ydiff/2)*prop,ymin+ydiff*6-(ydiff/2)*prop))
      xySubBb <- list(x = c(xmin+(xdiff/2)*prop,xmin+xdiff+(xdiff/2)*prop,xmin+xdiff*2+(xdiff/2)*prop,xmin+xdiff*3+(xdiff/2)*prop,xmin+xdiff*4+(xdiff/2)*prop,xmin+xdiff*5+(xdiff/2)*prop),
                      y = c(ymin+(ydiff/2)*prop,ymin+ydiff+(ydiff/2)*prop,ymin+ydiff*2+(ydiff/2)*prop,ymin+ydiff*3+(ydiff/2)*prop,ymin+ydiff*4+(ydiff/2)*prop,ymin+ydiff*5+(ydiff/2)*prop))

      xySubDF_1Ba <- as.data.frame(xySubBa)
      xySubDF_1Bb <- as.data.frame(xySubBb)


      xdiff <- (line3$x[2]-line3$x[1])/6
      ydiff <- (line3$y[2]-line3$y[1])/6

      xmin <- line3$x[1]
      ymin <- line3$y[1]

      xySubA <- list(x = c((xmin+xmin+xdiff)/2,(xmin+xdiff+xmin+xdiff*2)/2,(xmin+xdiff*2+xmin+xdiff*3)/2,(xmin+xdiff*3+xmin+xdiff*4)/2,(xmin+xdiff*4+xmin+xdiff*5)/2,(xmin+xdiff*5+xmin+xdiff*6)/2),y = c((ymin+ymin+ydiff)/2,(ymin+ydiff+ymin+ydiff*2)/2,(ymin+ydiff*2+ymin+ydiff*3)/2,(ymin+ydiff*3+ymin+ydiff*4)/2,(ymin+ydiff*4+ymin+ydiff*5)/2,(ymin+ydiff*5+ymin+ydiff*6)/2))

      xySubDF_3A <- as.data.frame(xySubA)

      xySubBa <- list(x = c(xmin+xdiff-(xdiff/2)*prop,xmin+xdiff*2-(xdiff/2)*prop,xmin+xdiff*3-(xdiff/2)*prop,xmin+xdiff*4-(xdiff/2)*prop,xmin+xdiff*5-(xdiff/2)*prop,xmin+xdiff*6-(xdiff/2)*prop),
                      y = c(ymin+ydiff-(ydiff/2)*prop,ymin+ydiff*2-(ydiff/2)*prop,ymin+ydiff*3-(ydiff/2)*prop,ymin+ydiff*4-(ydiff/2)*prop,ymin+ydiff*5-(ydiff/2)*prop,ymin+ydiff*6-(ydiff/2)*prop))
      xySubBb <- list(x = c(xmin+(xdiff/2)*prop,xmin+xdiff+(xdiff/2)*prop,xmin+xdiff*2+(xdiff/2)*prop,xmin+xdiff*3+(xdiff/2)*prop,xmin+xdiff*4+(xdiff/2)*prop,xmin+xdiff*5+(xdiff/2)*prop),
                      y = c(ymin+(ydiff/2)*prop,ymin+ydiff+(ydiff/2)*prop,ymin+ydiff*2+(ydiff/2)*prop,ymin+ydiff*3+(ydiff/2)*prop,ymin+ydiff*4+(ydiff/2)*prop,ymin+ydiff*5+(ydiff/2)*prop))

      xySubDF_3Ba <- as.data.frame(xySubBa)
      xySubDF_3Bb <- as.data.frame(xySubBb)



      xdiff <- (line2$x[2]-line2$x[1])/5
      ydiff <- (line2$y[2]-line2$y[1])/5

      xmin <- line2$x[1]
      ymin <- line2$y[1]

      xySubA <- list(x = c((xmin+xmin+xdiff)/2,(xmin+xdiff+xmin+xdiff*2)/2,(xmin+xdiff*2+xmin+xdiff*3)/2,(xmin+xdiff*3+xmin+xdiff*4)/2,(xmin+xdiff*4+xmin+xdiff*5)/2),
                     y = c((ymin+ymin+ydiff)/2,(ymin+ydiff+ymin+ydiff*2)/2,(ymin+ydiff*2+ymin+ydiff*3)/2,(ymin+ydiff*3+ymin+ydiff*4)/2,(ymin+ydiff*4+ymin+ydiff*5)/2))

      xySubDF_2A <- as.data.frame(xySubA)

      xySubBa <- list(x = c(xmin+xdiff-(xdiff/2)*prop,xmin+xdiff*2-(xdiff/2)*prop,xmin+xdiff*3-(xdiff/2)*prop,xmin+xdiff*4-(xdiff/2)*prop,xmin+xdiff*5-(xdiff/2)*prop),
                      y = c(ymin+ydiff-(ydiff/2)*prop,ymin+ydiff*2-(ydiff/2)*prop,ymin+ydiff*3-(ydiff/2)*prop,ymin+ydiff*4-(ydiff/2)*prop,ymin+ydiff*5-(ydiff/2)*prop))
      xySubBb <- list(x = c(xmin+(xdiff/2)*prop,xmin+xdiff+(xdiff/2)*prop,xmin+xdiff*2+(xdiff/2)*prop,xmin+xdiff*3+(xdiff/2)*prop,xmin+xdiff*4+(xdiff/2)*prop),
                      y = c(ymin+(ydiff/2)*prop,ymin+ydiff+(ydiff/2)*prop,ymin+ydiff*2+(ydiff/2)*prop,ymin+ydiff*3+(ydiff/2)*prop,ymin+ydiff*4+(ydiff/2)*prop))

      xySubDF_2Ba <- as.data.frame(xySubBa)
      xySubDF_2Bb <- as.data.frame(xySubBb)

      xdiff <- (line4$x[2]-line4$x[1])/5
      ydiff <- (line4$y[2]-line4$y[1])/5

      xmin <- line4$x[1]
      ymin <- line4$y[1]

      xySubA <- list(x = c((xmin+xmin+xdiff)/2,(xmin+xdiff+xmin+xdiff*2)/2,(xmin+xdiff*2+xmin+xdiff*3)/2,(xmin+xdiff*3+xmin+xdiff*4)/2,(xmin+xdiff*4+xmin+xdiff*5)/2),
                     y = c((ymin+ymin+ydiff)/2,(ymin+ydiff+ymin+ydiff*2)/2,(ymin+ydiff*2+ymin+ydiff*3)/2,(ymin+ydiff*3+ymin+ydiff*4)/2,(ymin+ydiff*4+ymin+ydiff*5)/2))

      xySubDF_4A <- as.data.frame(xySubA)

      xySubBa <- list(x = c(xmin+xdiff-(xdiff/2)*prop,xmin+xdiff*2-(xdiff/2)*prop,xmin+xdiff*3-(xdiff/2)*prop,xmin+xdiff*4-(xdiff/2)*prop,xmin+xdiff*5-(xdiff/2)*prop),
                      y = c(ymin+ydiff-(ydiff/2)*prop,ymin+ydiff*2-(ydiff/2)*prop,ymin+ydiff*3-(ydiff/2)*prop,ymin+ydiff*4-(ydiff/2)*prop,ymin+ydiff*5-(ydiff/2)*prop))
      xySubBb <- list(x = c(xmin+(xdiff/2)*prop,xmin+xdiff+(xdiff/2)*prop,xmin+xdiff*2+(xdiff/2)*prop,xmin+xdiff*3+(xdiff/2)*prop,xmin+xdiff*4+(xdiff/2)*prop),
                      y = c(ymin+(ydiff/2)*prop,ymin+ydiff+(ydiff/2)*prop,ymin+ydiff*2+(ydiff/2)*prop,ymin+ydiff*3+(ydiff/2)*prop,ymin+ydiff*4+(ydiff/2)*prop))

      xySubDF_4Ba <- as.data.frame(xySubBa)
      xySubDF_4Bb <- as.data.frame(xySubBb)

      if(is.null(colorCheckerXY)){
        points(xySubDF_1A, pch=20, col = 'green')
        points(xySubDF_1Ba, pch=20, col = 'red')
        points(xySubDF_1Bb, pch=20, col = 'red')

        points(xySubDF_2A, pch=20, col = 'green')
        points(xySubDF_2Ba, pch=20, col = 'red')
        points(xySubDF_2Bb, pch=20, col = 'red')

        points(xySubDF_3A, pch=20, col = 'green')
        points(xySubDF_3Ba, pch=20, col = 'red')
        points(xySubDF_3Bb, pch=20, col = 'red')

        points(xySubDF_4A, pch=20, col = 'green')
        points(xySubDF_4Ba, pch=20, col = 'red')
        points(xySubDF_4Bb, pch=20, col = 'red')
      }

      labels <- list(c(1,7,13,19,25),
                     c(2,8,14,20,26),
                     c(3,9,15,21,27),
                     c(4,10,16,22,28),
                     c(5,11,17,23,29),
                     c(6,12,18,24,30))

      # mR <- raster::as.matrix(imager::R(im))*255
      # mG <- raster::as.matrix(imager::G(im))*255
      # mB <- raster::as.matrix(imager::B(im))*255
      #
      # rR <- raster::raster(mR)
      # rG <- raster::raster(mG)
      # rB <- raster::raster(mB)
      #
      # extent(rR) <- c(0, dim(im)[2], 0, dim(im)[1])
      # extent(rG) <- c(0, dim(im)[2], 0, dim(im)[1])
      # extent(rB) <- c(0, dim(im)[2], 0, dim(im)[1])
      #
      # rR <- flip(t(rR),'y')
      # rG <- flip(t(rG),'y')
      # rB <- flip(t(rB),'y')

      xyTot <- c()
      for(e in 1:nrow(xySubDF_1A)){

        xyLine1 <- xySubDF_1A[e,]
        xyLine2 <- xySubDF_3A[6:1,][e,]

        xdiff <- (xyLine2$x-xyLine1$x)/5
        ydiff <- (xyLine2$y-xyLine1$y)/5

        xmin <- xyLine1$x
        ymin <- xyLine1$y

        xySub <- list(x = c((xmin+xmin+xdiff)/2,(xmin+xdiff+xmin+xdiff*2)/2,(xmin+xdiff*2+xmin+xdiff*3)/2,(xmin+xdiff*3+xmin+xdiff*4)/2,(xmin+xdiff*4+xmin+xdiff*5)/2),
                      y = c((ymin+ymin+ydiff)/2,(ymin+ydiff+ymin+ydiff*2)/2,(ymin+ydiff*2+ymin+ydiff*3)/2,(ymin+ydiff*3+ymin+ydiff*4)/2,(ymin+ydiff*4+ymin+ydiff*5)/2))

        xySubDF <- as.data.frame(xySub)

        if(is.null(colorCheckerXY)){
          text(xySubDF, label=labels[[e]], col = 'green')
        }

        xySubDFLabel <- cbind(xySubDF, label=labels[[e]])



        xyLine1a <- xySubDF_1Ba[e,]
        xyLine1b <- xySubDF_1Bb[e,]

        xyLine3a <- xySubDF_3Bb[6:1,][e,]
        xyLine3b <- xySubDF_3Ba[6:1,][e,]

        xdiffa <- (xyLine3a$x-xyLine1a$x)/5
        ydiffa <- (xyLine3a$y-xyLine1a$y)/5

        xdiffb <- (xyLine3b$x-xyLine1b$x)/5
        ydiffb <- (xyLine3b$y-xyLine1b$y)/5

        xmina <- xyLine1a$x
        ymina <- xyLine1a$y

        xminb <- xyLine1b$x
        yminb <- xyLine1b$y


        xySubAa <- list(x = c(xmina+xdiffa-(xdiffa/2)*prop,xmina+xdiffa*2-(xdiffa/2)*prop,xmina+xdiffa*3-(xdiffa/2)*prop,xmina+xdiffa*4-(xdiffa/2)*prop,xmina+xdiffa*5-(xdiffa/2)*prop),
                        y = c(ymina+ydiffa-(ydiffa/2)*prop,ymina+ydiffa*2-(ydiffa/2)*prop,ymina+ydiffa*3-(ydiffa/2)*prop,ymina+ydiffa*4-(ydiffa/2)*prop,ymina+ydiffa*5-(ydiffa/2)*prop))
        xySubAb <- list(x = c(xmina+(xdiffa/2)*prop,xmina+xdiffa+(xdiffa/2)*prop,xmina+xdiffa*2+(xdiffa/2)*prop,xmina+xdiffa*3+(xdiffa/2)*prop,xmina+xdiffa*4+(xdiffa/2)*prop),
                        y = c(ymina+(ydiffa/2)*prop,ymina+ydiffa+(ydiffa/2)*prop,ymina+ydiffa*2+(ydiffa/2)*prop,ymina+ydiffa*3+(ydiffa/2)*prop,ymina+ydiffa*4+(ydiffa/2)*prop))

        xySubBa <- list(x = c(xminb+xdiffb-(xdiffb/2)*prop,xminb+xdiffb*2-(xdiffb/2)*prop,xminb+xdiffb*3-(xdiffb/2)*prop,xminb+xdiffb*4-(xdiffb/2)*prop,xminb+xdiffb*5-(xdiffb/2)*prop),
                        y = c(yminb+ydiffb-(ydiffb/2)*prop,yminb+ydiffb*2-(ydiffb/2)*prop,yminb+ydiffb*3-(ydiffb/2)*prop,yminb+ydiffb*4-(ydiffb/2)*prop,yminb+ydiffb*5-(ydiffb/2)*prop))
        xySubBb <- list(x = c(xminb+(xdiffb/2)*prop,xminb+xdiffb+(xdiffb/2)*prop,xminb+xdiffb*2+(xdiffb/2)*prop,xminb+xdiffb*3+(xdiffb/2)*prop,xminb+xdiffb*4+(xdiffb/2)*prop),
                        y = c(yminb+(ydiffb/2)*prop,yminb+ydiffb+(ydiffb/2)*prop,yminb+ydiffb*2+(ydiffb/2)*prop,yminb+ydiffb*3+(ydiffb/2)*prop,yminb+ydiffb*4+(ydiffb/2)*prop))

        xySubAaDF <- as.data.frame(xySubAa)
        xySubAbDF <- as.data.frame(xySubAb)
        xySubBaDF <- as.data.frame(xySubBa)
        xySubBbDF <- as.data.frame(xySubBb)

        xySubRow <- cbind(xySubAaDF,xySubAbDF,xySubBaDF,xySubBbDF, label=labels[[e]])
        colnames(xySubRow) <- c('x1','y1','x2','y2','x3','y3','x4','y4','label')

        xyTot <- rbind(xyTot, xySubRow)

      }

      xyTot$imR <- NA
      xyTot$imG <- NA
      xyTot$imB <- NA

      xyTot <- xyTot[order(xyTot$label),]

      for(e in 1:nrow(xyTot)){

        print(paste('Calculating observed RGB values for patch', e, sep = ' '))

        if(is.null(colorCheckerXY)){
          polygon(c(xyTot$x1[e],xyTot$x2[e],xyTot$x4[e],xyTot$x3[e]),
                  c(xyTot$y1[e],xyTot$y2[e],xyTot$y4[e],xyTot$y3[e]), border = 'red')
        }

        outline <- rbind(c(xyTot$x1[e],xyTot$y1[e]),
                         c(xyTot$x2[e],xyTot$y2[e]),
                         c(xyTot$x4[e],xyTot$y4[e]),
                         c(xyTot$x3[e],xyTot$y3[e]))

        poly <- sp::Polygons(list(sp::Polygon(outline)),paste("r"))
        polyList  <- c(poly)
        polyNames <- c(paste("r"))
        sr <- sp::SpatialPolygons(polyList)
        srdf <- sp::SpatialPolygonsDataFrame(sr, data.frame(1:length(polyNames), row.names=polyNames))

        extrR <- raster::extract(rR, srdf)
        extrG <- raster::extract(rG, srdf)
        extrB <- raster::extract(rB, srdf)

        xyTot$imR[e] <- mean(extrR[[1]])
        xyTot$imG[e] <- mean(extrG[[1]])
        xyTot$imB[e] <- mean(extrB[[1]])
      }

      # Colorimetric values for ColorCheker targets
      l1 <- c(1, 116, 91, 76)
      l2 <- c(2, 192, 158, 141)
      l3 <- c(3, 113, 128, 160)
      l4 <- c(4, 104, 118, 72)
      l5 <- c(5, 141, 135, 181)
      l6 <- c(6, 144, 195, 179)
      l7 <- c(7, 79, 141, 173)
      l8 <- c(8, 243, 244, 240)
      l9 <- c(9, 234, 235, 232)
      l10 <- c(10, 221, 222, 220)
      l11 <- c(11, 209, 211, 209)
      l12 <- c(12, 203, 138, 65)
      l13 <- c(13, 178, 93, 158)
      l14 <- c(14, 183, 185, 184)
      l15 <- c(15, 158, 160, 159)
      l16 <- c(16, 125, 127, 127)
      l17 <- c(17, 99, 99, 99)
      l18 <- c(18, 89, 94, 173)
      l19 <- c(19, 227, 209, 65)
      l20 <- c(20, 72, 74, 73)
      l21 <- c(21, 46, 46, 46)
      l22 <- c(22, 25, 24, 24)
      l23 <- c(23, 11, 11, 12)
      l24 <- c(24, 179, 95, 108)
      l25 <- c(25, 162, 66, 68)
      l26 <- c(26, 110, 159, 84)
      l27 <- c(27, 61, 63, 150)
      l28 <- c(28, 214, 174, 66)
      l29 <- c(29, 175, 197, 84)
      l30 <- c(30, 93, 61, 108)

      ColorCheckerRGB <- as.data.frame(rbind(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,l23,l24, l25, l26, l27, l28, l29, l30))
      colnames(ColorCheckerRGB) <- c('label', 'sR', 'sG', 'sB')

      dat <- merge(xyTot, ColorCheckerRGB, by = 'label')


      print('Calculating polynomial regression...')

      sR <- dat$sR
      sG <- dat$sG
      sB <- dat$sB

      imR <- dat$imR
      imG <- dat$imG
      imB <- dat$imB

      modelR <- lm(sR ~ imR + imG + imB + imR^2 + imG^2 + imB^2)
      modelG <- lm(sG ~ imR + imG + imB + imR^2 + imG^2 + imB^2)
      modelB <- lm(sB ~ imR + imG + imB + imR^2 + imG^2 + imB^2)


      dfIm = data.frame(
        imR = matrix(mR, ncol=1),
        imG = matrix(mG, ncol=1),
        imB = matrix(mB, ncol=1)
      )

      print('Calibrating colors...')
      prR <- predict(modelR, dfIm)
      prG <- predict(modelG, dfIm)
      prB <- predict(modelB, dfIm)

      dfCal <- as.data.frame(cbind(prR,prG,prB))

      print('Rebuilding image...')
      Ri = matrix(dfCal$prR, nrow=dim(im)[1])
      Gi = matrix(dfCal$prG, nrow=dim(im)[1])
      Bi = matrix(dfCal$prB, nrow=dim(im)[1])

      imCal = array(dim=dim(im))
      imCal[,,,1] = Ri
      imCal[,,,2] = Gi
      imCal[,,,3] = Bi

      imCal <- imager::as.cimg(imCal)

      if(is.null(colorCheckerXY)){
        layout(matrix(c(1,2), 2, 1, byrow = TRUE))

        plot(imRed)
        plot(imCal)


        if(fixedCorners == FALSE & is.null(fixedModel)){
          todo <- readline(prompt="Press [enter] to continue and save image >>> ")
        }
      }
    }


    if(!is.null(fixedModel)){

      modelR <- fixedModel[[1]]
      modelG <- fixedModel[[2]]
      modelB <- fixedModel[[3]]

      dfIm = data.frame(
        imR = matrix(mR, ncol=1),
        imG = matrix(mG, ncol=1),
        imB = matrix(mB, ncol=1)
      )

      print('Calibrating colors...')
      prR <- predict(modelR, dfIm)
      prG <- predict(modelG, dfIm)
      prB <- predict(modelB, dfIm)

      dfCal <- as.data.frame(cbind(prR,prG,prB))

      print('Rebuilding image...')
      Ri = matrix(dfCal$prR, nrow=dim(imRed)[1])
      Gi = matrix(dfCal$prG, nrow=dim(imRed)[1])
      Bi = matrix(dfCal$prB, nrow=dim(imRed)[1])

      imCal = array(dim=dim(imRed))
      imCal[,,,1] = Ri
      imCal[,,,2] = Gi
      imCal[,,,3] = Bi

      imCal <- imager::as.cimg(imCal)

    }

    if(is.null(prepath)){
      im <- imager::save.image(imCal, paste(IDlist[n], '_calibrated', extension, sep=''), quality = 1)
    }
    else{
      im <- imager::save.image(imCal, paste(prepath,'/',IDlist[n], '_calibrated', extension, sep=''), quality = 1)
    }
  }

  return(list(modelR,modelG,modelB))
}
