#' Aligns images using \code{\link[RNiftyReg]{niftyreg}} utilities for automated image registration
#' and extracts color pattern using watershed segmentation. This function works interactively by
#' allowing to pick a starting pixel within each pattern element from which the watershed will
#' extract the pattern. This function works best for patterns with sharp boundaries.
#'
#' @param sampleList List of RasterStack objects.
#' @param target Image imported as RasterStack used as target for registration.
#' @param resampleFactor Integer for downsampling image used by \code{\link{redRes}}.
#' @param crop Vector c(xmin, xmax, ymin, ymax) that specifies the pixel coordinates to crop the
#'    original image.
#' @param useBlockPercentage Block percentage as used in \code{\link[RNiftyReg]{niftyreg}}
#'    (default = 75).
#' @param removebgR Integer indicating the range RGB treshold to remove from image (e.g. 100 removes
#'    pixels with average RGB > 100; default = NULL) for registration analysis. This works only to
#'    remove a white background.
#' @param maskOutline When outline is specified, everything outside of the outline will be masked for
#'    the color extraction (default = NULL).
#' @param cartoonID ID of the sample for which the cartoon was drawn and will be used for masking.
#' @param correct Correct image illumination using a linear model (default = FALSE).
#' @param blur Blur image for priority map extraction (default = TRUE).
#' @param sigma Size of sigma for Gaussian blurring (default = 5).
#' @param bucketfill Use a bucket fill on the background to fill holes (default = TRUE).
#' @param cleanP Integer to remove spurious areas with width smaller than cleanP (default = NULL).
#' @param splitC Integer to split selected patterns into connected components and remove ones with
#'    areas smaller than splitC (default = NULL).
#' @param plotTransformed Plot transformed image (default = FALSE).
#' @param plotCorrect Plot corrected image, corrected for illumination using a linear model
#'    (default = FALSE).
#' @param plotEdges Plot image gradient (default = FALSE).
#' @param plotPriority Plot priority map (default = FALSE).
#' @param plotWS Plot watershed result (default = FALSE).
#' @param plotBF Plot bucketfill (default = FALSE).
#' @param plotFinal Plot extracted patterns (default = FALSE).
#'
#' @return List of raster objects.
#'
#' @examples
#' \dontrun{
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '.jpg'
#'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#'
#' target <- imageList[[1]]
#'
#' outline_BC0077 <- read.table(paste(system.file("extdata",  package = 'patternize'),
#' '/BC0077_outline.txt', sep=''), header = FALSE)
#'
#' rasterList_regW <- patRegW(imageList, target, plotTransformed = FALSE, cartoonID = 'BC0077',
#'                           correct = TRUE, plotCorrect = FALSE, blur = FALSE, sigma = 2,
#'                           bucketfill = FALSE, cleanP = 0, splitC = 10, plotPriority = TRUE,
#'                           plotWS = FALSE, plotBF = FALSE, plotFinal = TRUE, removebgR = 100,
#'                           maskOutline = outline_BC0077)
#' }
#'
#' @export
#'
#' @import raster
#' @importFrom imager as.cimg imfill isoblur imgradient imsplit watershed bucketfill clean split_connected highlight enorm parany add
#' @importFrom purrr discard
#' @importFrom dplyr sample_n
#' @importFrom magrittr %>%
#' @importFrom graphics locator
#' @importFrom stats lm

patRegW <- function(sampleList,
                    target,
                    resampleFactor = NULL,
                    useBlockPercentage = 75,
                    crop = c(0,0,0,0),
                    removebgR = NULL,
                    maskOutline = NULL,
                    cartoonID = NULL,
                    correct = FALSE,
                    blur = TRUE,
                    sigma = 3,
                    bucketfill = TRUE,
                    cleanP = NULL,
                    splitC = NULL,
                    plotTransformed = FALSE,
                    plotCorrect = FALSE,
                    plotEdges = FALSE,
                    plotPriority = FALSE,
                    plotWS = FALSE,
                    plotBF = FALSE,
                    plotFinal = FALSE){

  declare( . <- as_nonstandard())
  rasterList <- list()

  # Crop target image
  if(!identical(crop, c(0,0,0,0))){

    targetExtRaster <- crop
    target <- raster::crop(target, targetExtRaster)
  }

  # Reduce resolution of target image
  if(!is.null(resampleFactor)){
    target <- redRes(target, resampleFactor)
  }

  # Average color channels to gray scale
  targetA <- apply(raster::as.array(target), 1:2, mean)

  # Remove bacground of target image
  if(is.numeric(removebgR)){

    targetA <- apply(targetA, 1:2, function(x) ifelse(x > removebgR, 0, x))
  }

  if(!is.null(maskOutline)){

    if(is.null(cartoonID)){
      print('Please specify cartoonID')
    }

    indx <- which(names(sampleList) == cartoonID)

    maskOutlineNew <- maskOutline
    extPicture <- raster::extent(sampleList[[indx]])
    maskOutlineNew[,2] <- extPicture[4]-maskOutlineNew[,2]
  }


  n <- 1
  todo <- 0

  # Run the loop for each sample
  while(n <= length(sampleList)){

    if(todo != 'r' && todo != 'a'){

      sStack <- sampleList[[n]]
      extRaster <- raster::extent(sStack)

      # Crop image
      if(!identical(crop, c(0,0,0,0))){

        extRaster <- crop
        sStack <- crop(sStack, extRaster)
      }

      sourceRaster <- redRes(sStack, 1)

      # Reduce resolution
      if(!is.null(resampleFactor)){
        sourceRaster <- redRes(sStack, resampleFactor)
      }

      sourceRasterK <- sourceRaster

      # Average color channels to gray scale
      sourceRaster <- apply(raster::as.array(sourceRaster), 1:2, mean)

      # Remove bacground for registration
      if(is.numeric(removebgR)){

        sourceRaster <- apply(sourceRaster, 1:2, function(x) ifelse(x > removebgR, 0, x))
      }

      # Calculate transformation
      result <- RNiftyReg::niftyreg(sourceRaster, targetA, useBlockPercentage=useBlockPercentage)

      # Apply transformation
      transformedMap <- RNiftyReg::applyTransform(RNiftyReg::forward(result), raster::as.array(sourceRasterK), interpolation=0)

      r1 <- raster::raster(transformedMap[1:nrow(transformedMap),ncol(transformedMap):1,1])
      # r1 <- focal(r1, w=matrix(1,nrow=3,ncol=3), fun=fill.na, pad = TRUE, na.rm = FALSE)

      r2 <- raster::raster(transformedMap[1:nrow(transformedMap),ncol(transformedMap):1,2])
      # r2 <- focal(r2, w=matrix(1,nrow=3,ncol=3), fun=fill.na, pad = TRUE, na.rm = FALSE)

      r3 <- raster::raster(transformedMap[1:nrow(transformedMap),ncol(transformedMap):1,3])
      # r3 <- focal(r3, w=matrix(1,nrow=3,ncol=3), fun=fill.na, pad = TRUE, na.rm = FALSE)


      imageTr <-raster::stack(r1,r2,r3)
      imageTr <- raster::flip(imageTr,'x')

      raster::extent(imageTr) <- raster::extent(sourceRasterK)

      imageTr[imageTr == 0] <- NA

      # Plot transformed raster
      if(plotTransformed){

        imageTr <- raster::flip(imageTr, 'y')

        x <- as.array(imageTr)/255
        cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
        uniqueCols <- unique(cols)
        x2 <- match(cols, uniqueCols)
        dim(x2) <- dim(x)[1:2]
        raster::image(t(x2), col=uniqueCols, yaxt='n', xaxt='n', main = paste(names(sampleList)[n],'transformed', sep=' '))

        imageTr <- raster::flip(imageTr, 'y')
      }

      # MaskOutline
      if(!is.null(maskOutline)){

        imageTr <- maskOutline(imageTr, maskOutline, refShape = 'target', flipOutline = 'y', crop = crop,
                               maskColor = 255, imageList = sampleList)

        cropEx <- c(min(maskOutlineNew[,1]), max(maskOutlineNew[,1]), min(maskOutlineNew[,2]), max(maskOutlineNew[,2]))

        imageTr <- raster::crop(imageTr, cropEx)
        raster::extent(imageTr) <- cropEx
      }

      imageTr[is.na(imageTr)] <- 0

      # Transform to imager format
      imA <- raster::as.array(imageTr)
      imR <- as.raster(imA, nrow = dim(imageTr)[1], ncol = dim(imageTr)[2], max = 255)
      im <- imager::as.cimg(imR)

      if(n == 1){
        im1 <- imager::as.cimg(imR)
      }

      # Correct illumination using linear model
      if(correct){

        if(plotCorrect){

          layout(t(1:2))
          par(mar=c(2, 2, 2, 2))
          plot(im, main = paste(names(sampleList)[n], 'original', sep=' '))
        }

        d <- as.data.frame(im)

        # Subsamble, fit a linear model
        m <- sample_n(d, 1e4) %>% lm(value ~ x*y, data=.)

        # Correct by removing the trend
        im <- im-predict(m,d)

        if(plotCorrect) plot(im, main = paste(names(sampleList)[n], 'corrected', sep=' '))
      }


      # # Blur image for color extraction
      # if(focal){
      #
      #   gf <- focalWeight(image, sigma, "Gauss")
      #
      #   rrr1 <- raster::focal(image[[1]], gf)
      #   rrr2 <- raster::focal(image[[2]], gf)
      #   rrr3 <- raster::focal(image[[3]], gf)
      #
      #   image <- raster::stack(rrr1, rrr2, rrr3)
      #
      #   image[is.na(image)] <- 255
      # }


      # Blur image for priority map extraction
      if(blur){

        edges <- isoblur(im, sigma) %>% imgradient("xy") %>% enorm %>% imsplit("c") %>% add %>% sqrt
      }

      else{

        edges <- imgradient(im,"xy") %>% enorm
      }

      layout(t(1:1))

      if(plotEdges) plot(edges)

      priority <- 1/(255+edges)

      if(plotPriority) plot(priority)


      # Extract seed points from image
      seed <- imager::imfill(dim=dim(im)) #Empty image
    }

    layout(t(1:1))
    par(mar=c(2, 2, 2, 2))

    plot(im)

    print("Choose points to identify patterns. Click outside image area to stop.")

    while(1){

      xy <- locator(n=1)

      if(as.numeric(xy)[1] > dim(im)[1] || as.numeric(xy)[1] < 0 || as.numeric(xy)[2] > dim(im)[2] || as.numeric(xy)[2] < 0){
        break
      }

      print(paste('x: ', as.character(xy)[1], 'y: ', as.character(xy)[2]))

      seed[as.numeric(xy[1]),as.numeric(xy[2]),1,1] <- 2
    }

    print("Choose at leat one point to identify background. Click outside image area to stop.")

    while(1){

      xy <- locator(n=1)

      if(as.numeric(xy)[1] > dim(im)[1] || as.numeric(xy)[1] < 0 || as.numeric(xy)[2] > dim(im)[2] || as.numeric(xy)[2] < 0){
        break
      }

      print(paste('x: ', as.character(xy)[1], 'y: ', as.character(xy)[2]))

      xyBF <- xy

      seed[as.numeric(xy[1]),as.numeric(xy[2]),1,1] <- 1
    }


    # Run watershed transform
    ws <- (imager::watershed(seed, priority) == 2) # multiply by 255 if RGB extracted

    if(plotWS) plot(ws)


    # Use a bucket fill on the background to fill holes (starts from last selected bg pixel)
    if(bucketfill){

      ws <- imager::bucketfill(ws, as.numeric(xyBF[1]), as.numeric(xyBF[2]), color=2) %>% {!( . == 2)}
    }

    if(plotBF) plot(ws)


    # Remove spurious areas with width < cleanP
    if(!is.null(cleanP)){

      ws <- imager::clean(ws,cleanP)
    }


    # Split into connected components and remove ones with areas < splitC
    if(!is.null(splitC)){

      ws <- imager::split_connected(ws) %>% purrr::discard(~ sum(.) < splitC) %>% parany
    }


    # Transform to raster
    wsR <- as.raster(ws)
    wsM <- matrix(t(col2rgb(wsR))[,1]/255, nrow = dim(im)[1], ncol = dim(im)[2])
    patternRaster <- raster::raster(t(wsM))
    raster::extent(patternRaster) <- raster::extent(0, dim(im)[1], 0, dim(im)[2])

    if(plotFinal){

      layout(t(1:2))

      plot(im, main="Watershed")
      highlight(ws, col="red")

      plot(patternRaster, legend = FALSE, bty = 'n', main = 'Extracted pattern')
    }

    todo <- readline(prompt="Press [enter] to continue, 'a' to add points 'r' to redo or 'q' to quit >>> ")

    if(todo == 'q') break

    if(todo == 'a') next

    if(todo == 'r'){

      seed <- imager::imfill(dim=dim(im))
      next
    }

    else{

      y <- raster::raster(ncol = dim(im1)[2], nrow = dim(im1)[1])
      extent(y) <- raster::extent(0, dim(im)[1], 0, dim(im)[2])

      patternRaster <- resample(patternRaster, y, method = 'ngb')

      raster::extent(patternRaster) <- raster::extent(target)

      rasterList[[names(sampleList)[n]]] <- patternRaster

      print(paste('sample', names(sampleList)[n], 'done and added to rasterList', sep=' '))

      n = n + 1
    }
  }

  return(rasterList)

}

###
# The following code is added to suppress raising a R CMD check note resulting from magrittr '.'
###
# Placeholder function.
as_nonstandard <- function()
{
  stop("This is a placeholder function, and should not be called directly.")
}

# Declare symbols for nonstandard evalution.
declare <- function(declaration)
{
  d <- substitute(declaration)
  p <- parent.frame()

  if (!identical(d[[1]], quote(`<-`)))
    stop("Declarations follow the <- syntax.", call. = FALSE)

  if (!is.symbol(d[[2]]))
    stop("LHS needs to be a symbol", call. = FALSE)

  sym <- as.character(d[[2]])

  msg <- sprintf("Symbol '%s' is declared for non-standard use.", sym)

  eval(
    substitute(
      delayedAssign(x, stop(m, call. = FALSE), environment(), environment()),
      list(m = msg, x = sym)), p, p)
}

