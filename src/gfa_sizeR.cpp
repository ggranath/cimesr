// gfa_size.cpp: cimes function to extract size 
// --- GNU Copyright (C) 1982-2014 by Jean-Michel Walter. Contributors: Alemu Gonsamo ---
// Modifications for R implementation by Gustaf Granath 2014
//
// This file is a part of cimesr

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <vector> //GG

#include <Rcpp.h> //GG
using namespace Rcpp; //GG

const double pi=3.14159265358979323846;

// [[Rcpp::export]] 

SEXP BMPPolarGroupingR(SEXP Vs, SEXP height, SEXP width, SEXP ZinfR, SEXP ZsupR, SEXP ZpasR)
  {

    // int BMPPolarGrouping(t_Bitmap *B, double Zinf, double Zsup, double Zpas, char *name)
    // {

    // FILE *histo;
    
  // GG section added 
  //read in R data
  NumericVector Vr(Vs);

  //make C++ variables
  int heightC = Rcpp::as<int>(height);
  int widthC = Rcpp::as<int>(width);
  double Zinf = Rcpp::as<double>(ZinfR);
  double Zsup = Rcpp::as<double>(ZsupR);
  double Zpas = Rcpp::as<double>(ZpasR);
  int centerX=widthC/2, centerY=heightC/2, radius=widthC/2;
  // end GG section
  
  // B->Header.biWidth replaced with widthC and B->Header.biHeight replaced by heightC
    
  int sequence[65536], seq1[8192], seq2[8192], seq3[8192], seq4[8192], seq5[8192], seq6[8192], seq7[8192], seq8[8192],cc,cp; //GG
  //unsigned char sequence[65536], seq1[8192], seq2[8192], seq3[8192], seq4[8192], seq5[8192], seq6[8192], seq7[8192], seq8[8192],cc,cp;
  
  int histoB[32768], histoW[32768];
  int seqsize, seq1size, seq2size, seq3size, seq4size, seq5size, seq6size, seq7size, seq8size;
	int x1,y1,x3,y3,x5,y5,x7,y7;
	// int centerX=B->Header.biWidth/2, centerY=B->Header.biHeight/2, radius=B->Header.biWidth/2;
	double mB,mW,sB2,sW2,pq,Sm=0;
	double Z;
  int r,h;
	int nB,nW,sumB,sumW,nSm=0,l;
    int npB,npW,ntp;
	int x,y,e,xback,yback,eback,j;

    //cout << "Polar Grouping Algorith: " << Zinf << " degree -> " << Zsup << " degree, +" << Zpas << " degree" << endl;
	Zinf=pi*Zinf/180; Zsup=pi*Zsup/180.0;
    if (Zpas)   Zpas=pi*Zpas/180.0;
    else        Zpas=pi*(90.0/radius)/180.0;

    //results=fopen("Grouping.dat", "wt");
    //fprintf(results,"Theta\t#Seq\tMean_L_White\tMean_L_Black\tp+q\t\t95%% CI\t\t#PixW\t#PixB\t#Pix\n\n");
    //histo=fopen("gapsize.txt", "wt"); GG removed
    
    //GG add vectors for R
    std::vector<double> colzenW;
    std::vector<double> colzenB;
    std::vector<int> colW;
    std::vector<int> colB;
    //

	for (Z=Zinf; Z<=(Zsup+(pi/2.0)/radius); Z+=Zpas)
	{
    	seqsize=seq1size=seq2size=seq3size=seq4size=seq5size=seq6size=seq7size=seq8size=0;
        histoB[0]=histoW[0]=0;

    	r=(double)radius*Z/(pi/2.0);
    	x=0;
    	y=r;

      e=3-2*r;
    	while (x<=y)
	    {
		    seq1[seq1size++]=Vr[widthC*(centerY-y)+centerX+x]; // 1

            xback=x; yback=y;

	    	if (e<0)
            {
                eback=e;
                e+=4*x+++6;
            }
    		else
            {
                eback=e;
	    		e+=4*(x++-y--)+10;
            }
    	}
    	x1=centerX+xback; y1=centerY-yback;

    	x=0;
    	y=r;
    	e=3-2*r;
    	while (x<=y)
    	{
	    	seq2[seq2size++]=Vr[widthC*(centerY-x)+centerX+y]; // 2
            xback=x; yback=y;

	    	if (e<0)
            {
                eback=e;
		    	e+=4*x+++6;
            }
    		else
            {
                eback=e;
	    		e+=4*(x++-y--)+10;
            }
    	}
    	if (eback<0)
        {
            if ((centerY-xback==y1) && (centerX+yback==x1)) seq2size--;
        }
    	else
        {
    		if ((centerY-xback==y1) && (centerX+yback==x1)) seq2size--;
        }

    	e=3-2*r;
    	if (e<0)
         	{ x=1; y=r; }
    	else
        	{ x=1; y=r-1; }
    	while (x<=y)
    	{
    		seq3[seq3size++]=Vr[widthC*(centerY+x)+centerX+y]; // 3
            xback=x; yback=y;

    		if (e<0)
            {
                eback=e;
    			e+=4*x+++6;
            }
    		else
            {
                eback=e;
    			e+=4*(x++-y--)+10;
            }
    	}
    	x3=centerX+yback; y3=centerY+xback;

    	x=0;
    	y=r;
    	e=3-2*r;
    	while (x<=y)
    	{
    		seq4[seq4size++]=Vr[widthC*(centerY+y)+centerX+x]; // 4
            xback=x; yback=y;

	    	if (e<0)
            {
                eback=e;
	    		e+=4*x+++6;
            }
	    	else
            {
                eback=e;
	    		e+=4*(x++-y--)+10;
            }
	    }

    	if (eback<0)
        {
    		if ((centerY+yback==y3) && (centerX+xback==x3)) seq4size--;
        }
	    else
        {
	    	if ((centerY+yback==y3) && (centerX+xback==x3)) seq4size--;
        }

    	e=3-2*r;
    	if (e<0)
        	{ x=1; y=r; }
	    else
        	{ x=1; y=r-1; }
    	while (x<=y)
    	{
    		seq5[seq5size++]=Vr[widthC*(centerY+y)+centerX-x]; // 5
            xback=x; yback=y;

    		if (e<0)
            {
                eback=e;
    			e+=4*x+++6;
            }
    		else
            {
                eback=e;
    			e+=4*(x++-y--)+10;
            }
    	}
    	x5=centerX-xback; y5=centerY+yback;

	    x=0;
    	y=r;
    	e=3-2*r;
	    while (x<=y)
    	{
	    	seq6[seq6size++]=Vr[widthC*(centerY+x)+centerX-y]; // 6
            xback=x; yback=y;

		    if (e<0)
            {
                eback=e;
    			e+=4*x+++6;
            }
    		else
            {
                eback=e;
	    		e+=4*(x++-y--)+10;
            }
    	}
    	if (eback<0)
        {
    		if ((centerY+xback==y5) && (centerX-yback==x5)) seq6size--;
        }
    	else
        {
    		if ((centerY+xback==y5) && (centerX-yback==x5)) seq6size--;
        }

    	e=3-2*r;
    	if (e<0)
        	{ x=1; y=r; }
    	else
        	{ x=1; y=r-1; }
    	while (x<=y)
    	{
        	seq7[seq7size++]=Vr[widthC*(centerY-x)+centerX-y]; // 7
            xback=x; yback=y;

	    	if (e<0)
            {
                eback=e;
	    		e+=4*x+++6;
            }
	    	else
            {
                eback=e;
		    	e+=4*(x++-y--)+10;
            }
    	}
	    x7=centerX-yback; y7=centerY-xback;

    	e=3-2*r;
    	if (e<0)
	        { x=1; y=r; }
    	else
	        { x=1; y=r-1; }
    	while (x<=y)
	    {
		    seq8[seq8size++]=Vr[widthC*(centerY-y)+centerX-x]; // 8
            xback=x; yback=y;

	    	if (e<0)
            {
                eback=e;
		    	e+=4*x+++6;
            }
    		else
            {
                eback=e;
	    		e+=4*(x++-y--)+10;
            }
    	}

    	if (eback<0)
        {
    		if ((centerY-yback==y7) && (centerX-xback==x7)) seq8size--;
        }
    	else
        {
    		if ((centerY-yback==y7) && (centerX-xback==x7)) seq8size--;
        }

    	for (j=0; j<seq1size; j++)
    		sequence[seqsize++]=seq1[j];
    	for (j=seq2size-1; j>=0; j--)
    		sequence[seqsize++]=seq2[j];
	    for (j=0; j<seq3size; j++)
	    	sequence[seqsize++]=seq3[j];
    	for (j=seq4size-1; j>=0; j--)
	    	sequence[seqsize++]=seq4[j];
    	for (j=0; j<seq5size; j++)
	    	sequence[seqsize++]=seq5[j];
    	for (j=seq6size-1; j>=0; j--)
	    	sequence[seqsize++]=seq6[j];
    	for (j=0; j<seq7size; j++)
	    	sequence[seqsize++]=seq7[j];
    	for (j=seq8size-1; j>=0; j--)
		    sequence[seqsize++]=seq8[j];

    	nW=nB=0;
    	sumW=sumB=0;
    	l=1;
    	cc=sequence[0];
    	for (j=1; j<seqsize; j++)
    	{
    		cp=cc;
    		cc=sequence[j];
    		if (cc==cp)
    			l++;
    		else
                if (cp) { histoW[nW]=l; nW++; sumW+=l; l=1; }
                else    { histoB[nB]=l; nB++; sumB+=l; l=1; }
    	}
    	if (cc)	sumW+=l;
    	else	sumB+=l;
    	if (cc!=sequence[0])
                if (cc) nW++;
                else    nB++;
        else    if (cc) histoW[0]+=l;
                else    histoB[0]+=l;

    	mB=(nB?(double)sumB/nB:0);
    	mW=(nW?(double)sumW/nW:0);
        sB2=(nB?(1.0/nB)*(mB-1)/(mB*mB*mB):0);
        sW2=(nW?(1.0/nW)*(mW-1)/(mW*mW*mW):0);
    	if (nB && nW)
	    {
	    	pq=1.0/mB+1.0/mW;
	    	Sm+=pq;
	    	nSm++;
    	}

    npB=npW=ntp=0;
      
    //fprintf(histo,"%2.1lf\n",180.0*Z/pi); GG remove
    double zz=180.0*Z/pi; //GG add zenith angle
    
    for (h=0; h<nW; h++)
    {
        if (histoW[h]>0)
		{
			//fprintf(histo,"%i ",histoW[h]); GG remove
      colW.push_back(histoW[h]); //GG add
      colzenW.push_back(zz); //GG add

      npW+=histoW[h];
		}
    }
    //fprintf(histo,"\n");
    for (h=0; h<nB; h++)
    {

        if (histoB[h]>0)
		{
			//fprintf(histo,"%i ",histoB[h]);
      colB.push_back(histoB[h]); //GG add
      colzenB.push_back(180.0*Z/pi); //GG add
      
      npB+=histoB[h];
		}			
    }
    //fprintf(histo,"\n"); GG remove
  
    ntp=npW+npB;
    if (nB && nW)
    {
        //cout << (180.0*Z/pi) << " degree : " << "p+q= " << pq << "  [" << (1.0-1.96*sqrt(sB2+sW2)) << " ," << (1.0+1.96*sqrt(sB2+sW2)) << "]" << endl;
        //fprintf(results,"%2.1lf\t%4i\t%lf\t%lf\t%lf\t%lf\t%4i\t%4i\t%4i\n",180.0*Z/pi,nB,mW,mB,pq,1.96*sqrt(sB2+sW2),npW,npB,ntp);
    }
    else
    {
        //cout << 180.0*Z/pi << " degree : " << "One color" << endl;
        //fprintf(results,"%2.1lf\t*\n",180.0*Z/pi);
    }

    } // au debut de la fonction

    if (nSm)
    {
        Sm/=nSm;

        //cout << endl << "Sampled circles: " << nSm << endl;
        //cout << "Mean: " << Sm << endl;

        //fprintf(results,"\n\nSampled circles: %i   Mean: %lf\n",nSm,Sm);
    }
    else
    {
        //cout << "Can't compute Grouping algorithm !" << endl;

        //fprintf(results,"\n\nCan't compute Grouping algorithm !\n");
    }

    //fclose(results);
    //fclose(histo); GG remove

	//return 0; GG remove
  
  //GG add export to R as list        
  SEXP res = R_NilValue;
  res = Rcpp::List::create(Rcpp::Named("ZenithWhite")  = wrap(colzenW),
                           Rcpp::Named("WhiteSeq")  = wrap(colW),
                           Rcpp::Named("ZenithBlack")  = wrap(colzenB),
                           Rcpp::Named("BlackSeq")  = wrap(colB));
  
  return res;
}
