#' Extracts color pattern from landmark transformed image using watershed segmentation. This function
#' works interactively by allowing to pick a starting pixel within each pattern element from which the
#' watershed will extract the pattern. This function works best for patterns with sharp boundaries.
#'
#' @param sampleList List of RasterStack objects.
#' @param landList Landmark list as returned by \code{\link[patternize]{makeList}}.
#' @param IDlist List of sample IDs should be specified when masking outline and transformRef
#'    is 'meanshape'.
#' @param adjustCoords Adjust landmark coordinates in case they are reversed compared to pixel
#'    coordinates (default = FALSE).
#' @param transformRef ID or landmark matrix of reference sample for shape to which color patterns
#'    will be transformed to. Can be 'meanshape' for transforming to mean shape of Procrustes
#'    analysis.
#' @param resampleFactor Integer for downsampling image used by \code{\link{redRes}}.
#' @param transformType Transformation type as used by \code{\link[Morpho]{computeTransform}}
#'    (default ='tps').
#' @param maskOutline When outline is specified, everything outside of the outline will be masked for
#'    the color extraction (default = NULL).
#' @param cartoonID ID of the sample for which the cartoon was drawn and will be used for masking
#'    (should be set when transformRef = 'meanShape').
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
#'
#' \dontrun{
#' IDlist <- c('BC0077','BC0071','BC0050','BC0049','BC0004')
#' prepath <- system.file("extdata",  package = 'patternize')
#' extension <- '_landmarks_LFW.txt'
#'
#' landmarkList <- makeList(IDlist, 'landmark', prepath, extension)
#'
#' extension <- '.jpg'
#' imageList <- makeList(IDlist, 'image', prepath, extension)
#'
#' outline_BC0077 <- read.table(paste(system.file("extdata",  package = 'patternize'),
#' '/BC0077_outline.txt', sep=''), header = FALSE)
#'
#' rasterList_W <- patLanW(imageList, landmarkList, IDlist, transformRef = 'meanshape',
#' adjustCoords = TRUE, plotTransformed = FALSE, correct = TRUE, plotCorrect = FALSE, blur = FALSE,
#' sigma = 2, bucketfill = FALSE, cleanP = 0, splitC = 10, plotPriority = TRUE, plotWS = TRUE,
#' plotBF = TRUE, plotFinal = TRUE, maskOutline = outline_BC0077, cartoonID = 'BC0077')
#' }
#'
#' @export
#' @import raster
#' @importFrom imager as.cimg imfill isoblur imgradient imsplit watershed bucketfill clean split_connected highlight enorm parany add
#' @importFrom purrr discard
#' @importFrom dplyr sample_n
#' @importFrom magrittr %>%
#' @importFrom graphics locator
#' @importFrom stats lm
#' @importFrom Morpho procSym computeTransform applyTransform

patLanW <- function(sampleList,
                    landList,
                    IDlist = NULL,
                    adjustCoords = FALSE,
                    transformRef = 'meanshape',
                    resampleFactor = NULL,
                    transformType = 'tps',
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


  # Check whether sampleList and landList have the same length
  if(length(sampleList) != length(landList)){
    stop("sampleList is not of the same length as lanArray")
  }

  for(n in 1:length(sampleList)){
    if(names(sampleList)[n] != names(landList)[n]){
      stop("samples are not in the same order in sampleList and lanArray")
    }
  }


  # Make landmark array
  lanArray <- lanArray(landList, adjustCoords, sampleList)


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

    indx <- which(names(sampleList) == cartoonID)
    maskOutlineNew <- maskOutline
    extPicture <- raster::extent(sampleList[[indx]])
    maskOutlineNew[,2] <- extPicture[4]-maskOutlineNew[,2]
    }
  }

  if(is.null(cartoonID)){
    maskOutlineNew <- maskOutline
  }

  if(all(c(!is.null(maskOutline), transformRef == 'meanshape'))){

    invisible(capture.output(cartoonLandTrans <- Morpho::computeTransform(refShape,
                                                                          as.matrix(lanArray[,,indx]),
                                                                          type="tps")))

    maskOutlineMean <- Morpho::applyTransform(as.matrix(maskOutlineNew), cartoonLandTrans)
  }


  n <- 1
  todo <- 0

  # Run the loop for each sample
  while(n <= length(sampleList)){

    if(todo != 'r' && todo != 'a'){

      image <- sampleList[[n]]
      extRaster <- raster::extent(image)

      # Reduce resolution
      if(!is.null(resampleFactor)){
        image <- redRes(image, resampleFactor)
      }


      # Transform image using landmarks
      invisible(capture.output(transMatrix <- Morpho::computeTransform(refShape,
                                                                       as.matrix(lanArray[,,n]),
                                                                       type = transformType)))

      imageDF1 <- raster::as.data.frame(image[[1]], xy = TRUE)
      imageDF2 <- raster::as.data.frame(image[[2]], xy = TRUE)
      imageDF3 <- raster::as.data.frame(image[[3]], xy = TRUE)

      invisible(capture.output(imageT <- Morpho::applyTransform(as.matrix(imageDF1)[,1:2], transMatrix)))

      r <- raster::raster(nrow = dim(image)[1], ncol = dim(image)[2])

      raster::extent(r) <- c(min(imageT[,1]),max(imageT[,1]),min(imageT[,2]),max(imageT[,2]))

      # Rasterize the transformed image and fill in NA values using
      imageT1r <- raster::rasterize(imageT, field = imageDF1[,3], r, fun = mean)
      imageT1rf <- focal(imageT1r, w=matrix(1,nrow=3,ncol=3), fun=fill.na, pad = TRUE, na.rm = FALSE)

      imageT2r <- raster::rasterize(imageT, field = imageDF2[,3], r, fun = mean)
      imageT2rf <- focal(imageT2r, w=matrix(1,nrow=3,ncol=3), fun=fill.na, pad = TRUE, na.rm = FALSE)

      imageT3r <- raster::rasterize(imageT, field = imageDF3[,3], r, fun = mean)
      imageT3rf <- focal(imageT3r, w=matrix(1,nrow=3,ncol=3), fun=fill.na, pad = TRUE, na.rm = FALSE)

      imageTr <- raster::stack(imageT1rf, imageT2rf, imageT3rf)

      imageTr[is.na(imageTr)] <- 255


      # Plot transformed raster
      if(plotTransformed){

        # imageTr <- raster::flip(imageTr, 'y')

        x <- as.array(imageTr)/255
        cols <- rgb(x[,,1], x[,,2], x[,,3], maxColorValue=1)
        uniqueCols <- unique(cols)
        x2 <- match(cols, uniqueCols)
        dim(x2) <- dim(x)[1:2]
        raster::image(t(x2), col=uniqueCols, yaxt='n', xaxt='n', main = paste(names(landList)[n],'transformed', sep=' '))

        # imageTr <- raster::flip(imageTr, 'y')
      }

      # MaskOutline
      if(!is.null(maskOutline)){

        if(transformRef[1] != 'meanshape'){
          imageTr <- maskOutline(imageTr, maskOutlineNew, refShape = 'target', crop = c(0,0,0,0),
                                maskColor = 255, imageList = sampleList)

          cropEx <- c(min(maskOutlineNew[,1]), max(maskOutlineNew[,1]), min(maskOutlineNew[,2]), max(maskOutlineNew[,2]))

        }
        if(transformRef[1] == 'meanshape'){

          imageTr <- maskOutline(imageTr, outline = maskOutline, refShape = 'mean',
                      IDlist = IDlist, landList = landList, imageList = sampleList,
                      adjustCoords = TRUE, cartoonID = cartoonID)

          cropEx <- c(min(maskOutlineMean[,1]), max(maskOutlineMean[,1]), min(maskOutlineMean[,2]), max(maskOutlineMean[,2]))
        }

        imageTr <- raster::crop(imageTr, cropEx)
        raster::extent(imageTr) <- cropEx
      }


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
          plot(im, main = paste(names(landList)[n], 'original', sep=' '))
        }

        d <- as.data.frame(im)

        # Subsamble, fit a linear model
        m <- sample_n(d, 1e4) %>% lm(value ~ x*y, data=.)

        # Correct by removing the trend
        im <- im-predict(m,d)

        if(plotCorrect) plot(im, main = paste(names(landList)[n], 'corrected', sep=' '))
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

      if(any(c(as.numeric(xy)[1] > dim(im)[1], as.numeric(xy)[1] < 0, as.numeric(xy)[2] > dim(im)[2], as.numeric(xy)[2] < 0))){
        break
      }

      print(paste('x: ', as.character(xy)[1], 'y: ', as.character(xy)[2]))

      seed[as.numeric(xy[1]),as.numeric(xy[2]),1,1] <- 2
    }

    print("Choose at least one point to identify background. Click outside image area to stop.")

    while(1){

      xy <- locator(n=1)

      if(any(c(as.numeric(xy)[1] > dim(im)[1], as.numeric(xy)[1] < 0, as.numeric(xy)[2] > dim(im)[2], as.numeric(xy)[2] < 0))){
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
      highlight(ws, col="blue")

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

      raster::extent(patternRaster) <- c(min(refShape[,1]), max(refShape[,1]), min(refShape[,2]), max((refShape[,2])))

      rasterList[[names(landList)[n]]] <- patternRaster

      print(paste('sample', names(landList)[n], 'done and added to rasterList', sep=' '))

      n = n + 1
    }
  }

  return(rasterList)
}

# Function used in focal to calculate average for NA values in tranformed matrix
fill.na <- function(x, i=5) {
  if( is.na(x)[i] ) {
    return( round(mean(x, na.rm=TRUE),0) )
  } else {
    return( round(x[i],0) )
  }
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
