// gfa_fraction.cpp: cimes function to extract white/black pixels in sky segments 
// GNU Copyright (C) 1982-2014 by Jean-Michel Walter. Contributors: Alemu Gonsamo
// Modifications for R implementation by Gustaf Granath 2014
//
// This file is a part of cimesr

#include <algorithm> 
#include <vector> 
#include <Rcpp.h>
using namespace Rcpp;

const double pi=3.14159265358979323846;

// [[Rcpp::export]]

SEXP BMPPolarProjectionR(SEXP Vs, SEXP h, SEXP w, SEXP nZenithR, SEXP nAzimutR, SEXP a0R, SEXP revNum)
{
  
  // int BMPPolarProjection(t_Bitmap *B, int nZenith, int nAzimut, double a0, char *name)
  // {

  //FILE *r;
  
  //read in R data
  NumericVector Vr(Vs);
  NumericVector rev(revNum);

  //make C++ variables
  int heightC = Rcpp::as<int>(h);
  int widthC = Rcpp::as<int>(w);
  int nZenith = Rcpp::as<int>(nZenithR);
  int nAzimut = Rcpp::as<int>(nAzimutR);
  double a0 = Rcpp::as<double>(a0R);
  int i,j;
  double d, x, y, a;
  long Total[90][360], Ciel[90][360];
  int centerX=widthC/2, centerY=heightC/2, radius=widthC/2;
	int omega, alpha;

	for (j=0; j<nZenith; j++)
		for (i=0; i<nAzimut; i++)
			Total[j][i]=Ciel[j][i]=0;

	for (j=0; j<heightC; j++)
		for (i=0; i<widthC; i++)
		{
			x=i-centerX;
			y=centerY-j;
			d=hypot(x,y);
			if (d<=radius)
			{
			omega=(d/radius)*nZenith;
                if (!(y || x)) a=0; else a=atan2(y,x);
		if (a<0) a+=2*pi;
                a-=pi/2+pi*a0/180.0;
                if (a<0) a+=2*pi;
			alpha=(a/(2*pi))*nAzimut;
  
		if (Vr[j*widthC+i])
				Ciel[omega][alpha]++;
		Total[omega][alpha]++;
			}
		}
      
  NumericVector  zen(nZenith*nAzimut), az(nZenith*nAzimut), cpix(nZenith*nAzimut), tpix(nZenith*nAzimut);
  NumericVector seg = NumericVector::create(nZenith, nAzimut);
  
  for (j=0; j<nZenith; j++)
		for (i=nAzimut-1; i>=0; i--) {    
      zen[j*nAzimut + rev[i]] = (double)90*j/nZenith+(double)45/nZenith;
      az[j*nAzimut + rev[i]] = (double)360*i/nAzimut+(double)180/nAzimut;
      cpix[j*nAzimut + rev[i]] = (Total[j][i]?(double)Ciel[j][i]/Total[j][i]:(double)0);
      tpix[j*nAzimut + rev[i]] = Total[j][i];
		}
  
  //export to R as list
  SEXP res = R_NilValue;
  res = Rcpp::List::create(Rcpp::Named("segments")  = seg,
                           Rcpp::Named("zenith")  = zen,
                           Rcpp::Named("azimut")  = az,
                           Rcpp::Named("countPix")  = cpix,
                           Rcpp::Named("totalPix")  = tpix);
  
	return res;
}
