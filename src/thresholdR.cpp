// thresholdR.cpp: threshold function to turn greyscale into black and white 
// GNU Copyright (C) 1982-2014 by Jean-Michel Walter. Contributors: Alemu Gonsamo
// Modifications for R implementation by Gustaf Granath 2014
//
// This file is a part of cimesr

#include <algorithm> 
#include <vector> 
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

SEXP thresholdR(SEXP Vs, SEXP Th, SEXP h, SEXP w){
  
  //read in R data
  NumericVector Vr(Vs), threshold(Th), height(h), width(w);
  
  //set up out vector
  NumericVector outVr = clone(Vr);
  
  //make C++ variables
  double thresholdC = Rcpp::as<double>(threshold);
  int heightC = Rcpp::as<int>(height);
  int widthC = Rcpp::as<int>(width);
  int i,j;

  //loop to to do the thresholding.Returns 0 or 1
  for (j=0; j<heightC; j++) 
  	for (i=0; i<widthC; i++)
			outVr[j*widthC+i]=int(Vr[j*widthC+i]>thresholdC);

  return outVr;
}
