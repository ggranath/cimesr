# This file is a part of cimesr

##' @title Canopy openness
##'
##' @description Calculate canopy geometry indices with the possibility to correct for slope and aspect.
##'
##' @details Openness indices are calculated for each zenith angle.
##' The difference between unweighted and weighted canopy openness (CO) 
##' is that unweighted is a simple average gap projected onto the image plane 
##' for a given solid angle, whereas weighted corresponds to the average gap 
##' back projected to the object hemisphere.
##'
##' ISF: Anderson's "diffuse (indirect) site factor", Chazdon and Field "weighted canopy openness". Cosine correction is applied to both ISF. 
##' 
##' UOC: "uniform overcast" sky luminance distribution model.
##' 
##' SOC: "standard overcast" sky luminance distribution model.
##' 
##' 
##' A clumping index (CI) is calculated for a zenith view interval 30|D60| (Leblanc and Chen, 2001),
##' as:
##' CI = mean log(gap fraction) / log (mean gap fraction)
##'
##'
##' @param fracList a list returned by the \code{fraction()} function
##' @param slope numeric; in degrees (0-90 deg)
##' @param aspect numeric; in degrees (0-360 deg)
##' @param imgId the name of the image/plot 
##' @return a list
##' @return (1) Zenith angle 
##' @return (2) mean gap (arithmetically, or linearly, averaged over azimuth), i.e. gap profile over the whole zenith view angles. 
##' @return (3) cumulated percent gap (gap total, GT, or unweighted canopy openness)
##' @return (4) cumulated percent weighted canopy openness (visible sky, CO)
##' @return (5) cumulated percent canopy cover (CC, canopy/crown closure)
##' @return (6) cumulated percent relative diffuse light (ISF, UOC model)
##' @return (7) cumulated percent relative diffuse light (ISF, SOC model)
##' @return (8) Image id
##'
##' @references Leblanc S & Chen JM (2001) A practical scheme for correcting multiple scattering effects 
##' on optical LAI measurements. Agricultural and Forest Meteorology 110, 125-139.  
##' 
##' Gonsamo A, Walter J-M and Pellikka P  (2011) CIMES: A package of programs for determining 
##' canopy geometry and solar radiation regimes through hemispherical photographs. 
##' Computers and Electronics in Agriculture 79, 207-215. 
##' 
##' Walter J-MN (2009) CIMES (C) A package of programs for the Assessment of Canopy Geometry 
##' through Hemispherical Photographs. Manuel. Universite Louis Pasteur, Institut de Botanique, 
##' Strasbourg I.
##' 
##' @author Gustaf Granath, Alemu Gonsamo, Jean-Michel Walter
##.'
##' @export 
##'
##' @keywords openness
##'
##' @examples
##'
##' # Read image and convert to black and white
##' # Extract fish-eye area
##' # Define three points on the fish-eye ring
##' pointA <- c(x = 140, y = 102)
##' pointB <- c(x = 401.5, y = 25.9)
##' pointC <- c(x = 506.6, y = 238.6)
##' DFpoints <- data.frame(A = pointA, B = pointB, C = pointC)
##' blue <- ReadImg(system.file("pictures", "fisheye.jpg", package="cimesr"), 
##'                 channel = 3, coord = DFpoints, crop = TRUE)
##' bw <- threshold(blue, thresValue = 0.73)
##' frac <- fraction(bw, Zenith = 18, Azimut = 24, declination = 0)
##' gap.open <- openness(frac, slope = 0, aspect = 0, imgId = "plot1")
##' 
##' # display result
##' round(data.frame(gap.open[1:8]), 3) 


openness <- function (fracList, slope = 0, aspect = 0, imgId = "1") {
  opennessR(fracList, slope, aspect, imgId)
}
