# This file is a part of cimesr

##' @title Gap fraction extractor
##'
##' @description Gap fraction extractor (white pixels) for each sky segment which are defined by number of zenith annuli (vertical plane) and azimuth sectors (horizontal plane).
##'
##' @details The number of segments are controlled by the \code{Zenith}  and \code{Azimut} arguments. A recomended division is \code{Zenith}=18 and \code{Azimut}=24 which are the default settings.
##'
##' @param mat a matrix with 0s and 1s
##' @param Zenith integer; the number of divisions in the vertical plane
##' @param Azimut integer; the number of sectors in the horizontal plant
##' @param declination the magnetic declination in degrees at the site 
##' @return a list
##' (1) zenith angles, in degrees, mid-point values. 
##' (2) azimuth angles, degrees, mid-point values 
##' (3) gap fractions 
##' (4) total number of pixels in each sky division. 
##'
##' @author  Gustaf Granath, Alemu Gonsamo, Jean-Michel Walter
##.'
##' @export 
##'
##' @keywords gap extraction
##'
##' @seealso gap size extraction will be added 
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
##' # Extract gap fractions
##' frac <- fraction(bw, Zenith = 18, Azimut = 24, declination = 0)


fraction <- function (mat, Zenith = 18, Azimut = 24, declination  = 0) {
  matVec = array(t(mat)) # need to convert matrix to vector by row for the C++ function
  revNum = c((Azimut - 1):0)
  frac <- BMPPolarProjectionR(matVec, nrow(mat), ncol(mat), Zenith, Azimut, declination, revNum)
  class(frac) <- "cimesr"
  return(frac)
}