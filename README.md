cimesr: analysing canopy geometry and solar radiation regimes using hemispherical photographs
===========

Version: 0.1-dev

## Development

The cimesr package is now in development. It is an R implementation of the software bundle CIMES-FISHEYE which was developed by Jean-Michel Walter. 

## Features

At the moment only the canopy openness function is available. This function calulcates gap distribution, canopy openness, canopy closure, fraction diffuse light and can correct for aspect and slope.

## Installation

Development version from Github:
```
library("devtools"); install_github("cimesr", "ggranath", dependencies=TRUE)
```
This requires `devtools` (`install.packages("devtools")`)and the approach builds the package from source. Adding `build_vignettes=FALSE` can help if you have trouble installing.

## Citation
To cite `cimesr`, run `citation("cimesr")`

## Usage
```
# Define three points on the fish-eye ring
pointA <- c(x = 140, y = 102)
pointB <- c(x = 401.5, y = 25.9)
pointC <- c(x = 506.6, y = 238.6)
DFpoints <- data.frame(A = pointA, B = pointB, C = pointC)

# Read fish-eye image. Blue channel is extracted (greyscale). 
blue <- ReadImg(system.file("pictures", "fisheye.jpg", package="cimesr"), 
                 channel = 3, coord = DFpoints, crop = TRUE)
# Convert to black and white. Here using threshold value 0.73
bw <- threshold(blue, thresValue = 0.73)

# Extract gap fraction and estimate openness attributes
frac <- fraction(bw, Zenith = 18, Azimut = 24, declination = 0)
gap.open <- openness(frac, slope = 0, aspect = 0, imgId = "plot1")
 
# display result
round(data.frame(gap.open[1:8]), 3) 
```