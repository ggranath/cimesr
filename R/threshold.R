# This file is a part of cimesr

##' @title Greyscale to black and white
##'
##' @description Convert a greyscale image (matrix) into black and white using a set threshold.
##'
##' @details A threshold value can be estimated from imaging software or by using the /code{AutoThres()} function.
##'
##' @param mat numeric matrix 
##' @param thresValue numeric; a value defining the threshold between white and black. Must be between max and min value of /code{mat}.
##' 
##' @return matrix with cell values 0 (white) or 1 (black).
##'
##' @author  Gustaf Granath, Alemu Gonsamo, Jean-Michel Walter
##.'
##' @export 
##'
##' @examples
##'
##' # Read image and convert to black and white
##' blue <- readbitmap::read.bitmap(system.file("pictures", "fisheye.jpg", 
##'                                 package="cimesr"), channel = 3)
##' bw <- threshold(blue, thresValue = 0.73)

threshold <- function (mat, thresValue) {
  if(thresValue > 1| thresValue < 0) stop("A Threshold value between 0 and 1 must be provided")
  thresholdR(mat, thresValue, nrow(mat), ncol(mat))
}
