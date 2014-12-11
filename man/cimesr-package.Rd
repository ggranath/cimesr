\name{cimesr-package}
\alias{cimesr-package}
\alias{cimesr}
\docType{package}
\title{
cimesr: analysing canopy geometry and solar radiation regimes using
hemispherical photographs
}
\description{
This package is for analyzing structural properties of 
forest canopies using digital hemispherical photographs (fish-eye images). At the moment only 
functions to estimate openness are included but eventually you will be able to estimate canopy attributes like 
gap size, morphology and distribution, and vegetation indices like canopy and site openness, leaf area index (LAI), 
leaf angle distribution and clumping factors. Cimesr is build on the software CIMES-FISHEYE which was mainly developed 
by Jean-Michel Walter.
}
\details{
\tabular{ll}{
Package: \tab cimesr\cr
Type: \tab Package\cr
Version: \tab 0.1-0\cr
Date: \tab 2014-11-17\cr
License: \tab GPL 3\cr
}
Currently only a few base functions are available in cimesr. The user can load an image (jpg or bmp) and extract the blue channel of the fish-eye image. The image can then be transformed into a black and white image (a matrix with 0s and 1s) setting a threshold. A sugested threshold can be retrived with the AutoThres function. The number of segments, which is a functions of The number of azimuth and zenith sections in a polar projection, can be defined and gap fractions are calcuated for each sky segment. The declination at the site can also be defined. Openness can then be derived taking slope into account.
}
\author{
Gustaf Granath[aut, cre], Alemu Gonsamo [aut],  Jean-Michel Walter [aut]

Maintainer: Gustaf Granath <gustaf.granath@gmail.com>
}
\references{
Gonsamo, A., Walter, J. M. N., & Pellikka, P. (2011). CIMES: A package of programs for determining canopy geometry and solar radiation regimes through hemispherical photographs. Computers and Electronics in Agriculture, 79(2), 207-215.

Walter J-MN (2009) CIMES-FISHEYE. Hemispherical Photography of Forest Canopies. A package of
Programs for the Determination of Canopy Geometry and Solar Radiation Regimes through Hemispherical
Photographs. Manual. Universite de Strasbourg, France. http://jmnw.free.fr/
}

\keyword{Canopy openness}
\keyword{Gap fraction}
\keyword{Hemispherical photography}
\keyword{LAI}
\keyword{Solar radiation}
