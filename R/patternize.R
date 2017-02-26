#' patternize - An R package for quantifying color pattern variation.
#'
#' Quantifying variation in color patterns to study and compare the consistency of their expression necessitates the homologous alignment and color-based segmentation of images. Patternize is an R package that quantifies variation in color patterns as obtained from image data. Patternize defines homology between pattern positions across specimens either through fixed landmarks or image registration. Pattern identification is performed by categorizing the distribution of colors using either an RGB threshold or an unsupervised image segmentation. The quantification of the color patterns can be visualized as heat maps and compared between sets of samples.
#'
#' @author Steven M. Van Belleghem
#'
#' @section patternize main functions:
#'
#' The package has four main functions depending on how you want the alignment of the iamges and the color extraction to be performed.
#'
#' \code{patLanRGB} \cr
#'    Aligns images usings transformations obtained from fixed landmarks and extracts colors using a predefined RGB values and cutoff value.
#'
#' \code{patLanK} \cr
#'    Aligns images usings transformations obtained from fixed landmarks and extracts colors using k-means clustering.
#'
#' \code{patRegRGB} \cr
#'    Aligns images using \code{\link[RNiftyReg]{RNiftyReg}} utilities for automated image registration and extracts colors using a predefined RGB values and cutoff value.
#'
#' \code{patRegK} \cr
#'    Aligns images using \code{\link[RNiftyReg]{RNiftyReg}} utilities for automated image registration and extracts colors using k-means clustering.
#'
#'
#'
#' @section patternize preprocessing functions:
#'
#' The input for the main patternize functions are \code{\link[raster]{RasterStack}} objects and when landmark transformation is used, landmark arrays.
#'
#' \code{makeList} \cr
#'    This function returns a list of RasterStacks or a list of landmarks depending on the input provided.
#'
#' \code{lanArray} \cr
#'    This function creates a landmark array as used by \code{\link[Morpho]{procSym}} in the package \code{Morpho}.
#'
#'
#'
#' @section patternize postprocessing functions:
#'
#' \code{sumRaster} \cr
#'    This function sums the individual color pattern rasters as obtained by the main patternize functions.
#'
#' \code{plotHeat} \cr
#'    Plots the color pattern heatmaps. Uses \code{sumRaster} output.
#'
#' \code{patPCA} \cr
#'    This function transforms the individual color pattern rasters as obtained by the main patternize functions to a dataframe of 0 and 1 values that can be used for Principal Component Analysis (\code{\link[stats]{prcomp}}).
#'
#' \code{patArea} \cr
#'    This fucntion calculates the area in which the color pattern is expressed in each sample as the relative proportion using the provided outline of the considered trait or structure.
#'
#'
#'
#' @section patternize miscellaneous functions:
#'
#' \code{redRes} \cr
#'    Reduces the resolution of the \code{\link[raster]{RasterStack}} objects to speed up analysis.
#'
#' \code{kImage} \cr
#'    Performs k-means clustering of images.
#'
#' \code{createTarget} \cr
#'    Creates an artificial target images using a provided outline that can be used for image registration (experimantal).
#'
#'
#' @docType package
#' @name patternize
NULL
