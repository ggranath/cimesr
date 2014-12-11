# The imgExtract function is based on C++ code in the cimes software. GNU Copyright (C) 1982-2014 by Jean-Michel Walter. Contributors: Alemu Gonsamo  
# This file is a part of cimesr

##' @title Read jpg/bmp image
##'
##' @description Reads a hemispherical (fish-eye) jpg/bmp file. Returns a specified color channel and the image area is cropped out.
##'
##' @details This function calls the \code{read.bitmap()} function in the \code{readbitmap} package. The recomended channel to use for canopy analyses is the blue channel. This function crop the image so a qudratic image of the fish-eye area is returned. Although cropping is needed for any gap analyses in cimesr, it can be disabled by setting \code{crop = FALSE}. The whole image is then returned. The output is always a matrix with cell values from 0 to 1 and represents a greyscale.
##'
##' @param filename character; a jpg or bmp file
##' @param channel integer (1-3); which RGB channel to store, 3 = blue.
##' @param coord data frame or matrix with three pixel coordinates. Three columns where row one has x coordinates and row two y coordinates. See example.
##' @param crop logical; if TRUE, the hemispherical image is cropped out. 
##' 
##' @return a matrix with cell values between 0 and 1.
##'
##' @author Gustaf Granath
##.'
##' @export 
##'
##' @keywords read image
##'
##' @examples
##' # Define three points on the fish-eye ring
##' pointA <- c(x = 128.4, y = 92.9)
##' pointB <- c(x = 381.7, y = 17.3)
##' pointC <- c(x = 494.9, y = 225.4)
##' DFpoints <- data.frame(A = pointA, B = pointB, C = pointC)
##' blue <- ReadImg(system.file("pictures", "fisheye.jpg", package="cimesr"), 
##'                 channel = 3, coord = DFpoints, crop = TRUE)


ReadImg <- function(filename, channel = 3, coord = NULL, crop = TRUE) {
  if(crop == TRUE) if(!(all(dim(coord) == c(2, 3)))) stop("Three x, y coordinates are not given")
  mat <- readbitmap::read.bitmap(filename, channel = channel)
  if(crop == FALSE) return(mat)
  
  # imgExtract function  
  xa = coord[1, 1]
  xb = coord[1, 2]
  xc = coord[1, 3]
  ya = coord[2, 1]
  yb = coord[2, 2]
  yc = coord[2, 3]
  
  a = xa-xb
  b = ya-yb
  c = (xb*xb-xa*xa+yb*yb-ya*ya) / 2
  d = xa-xc
  e = ya-yc
  f = (xc*xc-xa*xa+yc*yc-ya*ya) / 2
  
  dd = a*e - b*d
  
  xm = (b*f-c*e) / dd
  ym = (c*d-a*f) / dd
  radius = sqrt((xc-xm)^2 + (yc-ym)^2)
  
  cat("Center: ", "x =", xm," y =", ym, "\n")
  cat("Radius: ", radius)
  
  # get new image measurements
  centerPoint <- c(xm, ym)
  xStart <- floor(centerPoint[1] - radius)
  xEnd <- floor(centerPoint[1] + radius)
  yStart <- floor(centerPoint[2] - radius)
  yEnd <- floor(centerPoint[2] + radius)
  # adjust the image if radius falls outside the orignal image
  # here black rows/columns are added if needed
  if(yStart < 0) {
    mat = rbind(mat, matrix(0, nrow = -1*yStart, ncol = ncol(mat)))
    yEnd = yEnd + (-1*yStart)
    yStart = 0
  }
  if(yEnd > dim(mat)[1]) {
    mat = rbind(mat, matrix(0, nrow = (yEnd-dim(mat)[1]), ncol = ncol(mat)))
    #yEnd = nrow(mat)
  }
  if(xStart < 0) {
    mat = cbind(mat, matrix(0, nrow = nrow(mat), ncol = -1*xStart))
    xEnd = xEnd + -1*xStart
    xStart = 0
  }
  if(xEnd > dim(mat)[2]) {
    mat = cbind(mat, matrix(0, nrow = nrow(mat), ncol = (xEnd-dim(mat)[2])))
    #xEnd = ncol(mat)
  }
  print(c(xStart, xEnd, yStart, yEnd))
  print(dim(mat))
  imgCrop <- mat[yStart:yEnd, xStart:xEnd]
  return(imgCrop)
}