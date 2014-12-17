# This file is a part of cimesr

##' @title Gap size extractor
##'
##' @description Extract gap sizes along 'circular transects' at a range of zenith angles. 
##'
##' @details Sequences of black and white pixels start at the first edge from gap to non-gap or vice-
##' versa, for each zenith angle. This is called foliage clumping. Black pixels
##' represent foliage or sky obstruction, whereas white pixels represent sky.
##'
##' @param mat a matrix with 0s and 1s
##' @param MinZenith numeric; minimum zenith angle in degrees (0 to 90)
##' @param MaxZenith numeric;  maximum zenith angle in degrees (0 to 90)
##' @param step  numeric; angle step (1 degree is recomended)
##' @return a list
##' 
##' \code{ZenithWhite} zenith angle, degrees, corresponding to the white sequence 
##' 
##' \code{WhiteSeq} lengths (in pixel units) of sequences of white pixels
##' 
##' \code{ZenithBlack} zenith angle, degrees, corresponding to the black sequence 
##' 
##' \code{BlackSeq} lengths (in pixel units) of sequences of black pixels
##' 
##' @author  Gustaf Granath, Alemu Gonsamo, Jean-Michel Walter
##.'
##' @export 
##'
##' @keywords gap extraction
##'
##' @examples
##'
##' # Define three points on the fish-eye ring
##'pointA <- c(x = 140, y = 102)
##'pointB <- c(x = 401.5, y = 25.9)
##'pointC <- c(x = 506.6, y = 238.6)
##'DFpoints <- data.frame(A = pointA, B = pointB, C = pointC)
##' # Read image and convert to black and white
##' blue <- ReadImg(system.file("pictures", "fisheye.jpg", package="cimesr"), 
##'                 channel = 3, coord = DFpoints, crop = TRUE)
##' bw <- threshold(blue, thresValue = 0.73)
##' # Extract gap size
##' siz <- size(bw,  MinZenith = 0, MaxZenith = 90, step = 1)


size <- function (mat, MinZenith = 0, MaxZenith = 90, step = 1) {
  matVec = array(t(mat)) # need to convert matrix to vector by row for the C++ function
  #revNum = c((Azimut - 1):0)
  siz <- BMPPolarGroupingR(matVec, nrow(mat), ncol(mat), MinZenith, MaxZenith, step)
  class(siz) <- "cimesr"
  return(siz)
}