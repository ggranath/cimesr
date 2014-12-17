# Auto thresholding function to get threshold value
# Copyright (C) 2014 Gustaf Granath
# This file is a part of cimesr


##' @title Find greyscale threshold value
##'
##' @description Function to automaticlly find a threshold value  to turn greyscale image into black and white. EXPERIMENTAL!
##'
##' @details A threshold value is determined by assuming a bimodal distribution of the greyscale values. One peak towards the white end and one towards the black end of the scale. The function divides the greyscale into 100 steps and starts in between the peaks. The mean of three steps ahead is then compared with the previous three steps. At a specific change (50 percent increase by default) the second peak is assumed to begin and the value at this point is returned.
##'
##' @param mat numeric matrix; can also be a greyscale image of the class "pixmapGrey" 
##' @param start numeric; set starting point of the search for a second peak. A value between 0 and 100 indicating the scale from black to white.
##' @param limit numeric; the increase required to define a second peak. Default value is 1.5 which is a 50 percent increase 
##' @param plot logical; show a histogram of the greyscale whith the suggested threshold 
##' 
##' @return threshold value
##'
##' @author Gustaf Granath
##'
##' @export 
##'
##' @examples
##'
##' # Read image and find threshold
##' blue <- readbitmap::read.bitmap(system.file("pictures", "fisheye.jpg", 
##'                                 package="cimesr"), channel = 3)
##' AutoThres(blue)

AutoThres <- function (mat, start = 30, limit = 1.5, plot = FALSE) {
  if (class(mat) == "pixmapGrey") mat = mat@grey
  dens <- hist(mat, breaks = 100, plot = FALSE)
  
  change = 0 
  while (change < limit) { # threshold (limit) is set to 50% increase by default
    m1 <- mean(dens$density[start:(start + 2)]) #mean over 3 values
    m2 <- mean(dens$density[(start + 3):(start + 5)]) # and the mean over the next 3 values
    start = start +3
    print(start)
    if(m2 == 0 | m1==0) { 
                          if(start > 90) {change = 2}
                                          else
                                          {change = 1}
                          } 
          else 
            {change <- m2/m1} # calc change and then check if its more than a 50% increase    
  }
  threshold <- start/100
  if(plot == TRUE) {
    hist(mat, breaks = 100)
    abline(v  = threshold, col = "red")
  }
  return(threshold)
}
