// opennessR.cpp: OPENNESS function in cimes ported to R 
// GNU Copyright (C) 1982-2014 by Jean-Michel Walter. Contributors: Alemu Gonsamo
// Modifications for R implementation by Gustaf Granath 2014
//
// This file is a part of cimesr

/******************************************************************************
CIMES - CANOPY GEOMETRY                                                                        
                                < OPENNESS >  (C)                                    
                       Programmed by Jean-Michel Walter                          
                                                                             
  This program computes the:

        (1) mean gap (arithmetically, or linearly, averaged over azimuth)
      i.e. gap profile over the whole zenith view angles,
        (2) cumulated percent gap (gap total, GT, or unweighted canopy openness),
        (3) cumulated percent weighted canopy openness (visible sky, CO),
        (4) cumulated percent canopy cover (CC, canopy/crown closure),
        (5) cumulated percent relative diffuse light (ISF, UOC model),
        (6) cumulated percent relative diffuse light (ISF, SOC model),

  for each zenith angle. The difference between unweighted and weighted CO is
  that unweighted is a simple average gap projected onto the image plane for a given
  solid angle, whereas weighted corresponds to the average gap back projected to the objet
  hemisphere.

  ISF: Anderson's "diffuse (indirect) site factor", Chazdon and Field "weighted canopy
  openness". Cosine correction is applied to both ISF. 
  UOC: "uniform overcast" sky luminance distribution model.
  SOC: "standard overcast" sky luminance distribution model.
  
  A clumping index (CI) is calculated for a zenith view interval 30¡Ð60¡ (Leblanc and Chen, 2001),
  as:
            CI = mean log(gap fraction) / log (mean gap fraction)


  Usage:    --> OPENNESS [parameter file][gap fraction data file][output file] <--
  ------

  Parameter file:
  --------------
        (1) Slope, degrees (double), e.g. 25.0
        (2) Aspect, degrees (double), e.g. 180.0

  Note:
  -----
        Zenith and horizon ring masks are not useful in this program since CO
        calculated for the whole hemisfere.

  
  References:
  ----------
  - Leblanc S & Chen JM (2001) A practical scheme for correcting multiple scattering effects 
    on optical LAI measurements. Agricultural and Forest Meteorology 110, 125Ð139.  
  - Gonsamo A, Walter J-M and Pellikka P  (2011) CIMES: A package of programs for determining 
    canopy geometry and solar radiation regimes through hemispherical photographs. 
    Computers and Electronics in Agriculture 79, 207Ð215. 
  - Walter J-MN (2009) CIMES (C) A package of programs for the Assessment of Canopy Geometry 
    through Hemispherical Photographs. Manuel. Universite Louis Pasteur, Institut de Botanique, 
    Strasbourg I.

  1st version : Feb 18 1999
  2nd    "    : Mar 15 1999
  3rd    "    : Jul 31 2000
  4th    "    : Sep 22 2000
  5th    "    : Jul 01 2008

*****************************************************************************/
#include <Rcpp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
using namespace Rcpp;

#define DIM         8     /* array dimension */
#define COL         4     /* number of columns (gap fraction file)    */
#define LIN      6000     /* max. number of lines (gap fraction file) */
#define MAX_ZEN    40     /* max. number of zenithal rings            */
#define MAX_AZI   150     /* max. number of azimuthal sectors         */

#define DEG_TO_RAD (PI/180.0)     /* conversion degrees to radians         */
#define RAD_TO_DEG (180.0/PI)     /* conversion radians to degrees         */
#define OVERCAST   2.0            /* value of "b" for sky radiance distribution
                                  under standard overcast sky (SOC), where
                                  (1+b) is the ratio of radiance at the zenith
                                  to that at the horizon. For British Isles,
                                  b=1.23 (Stevens and Unsworth, 1980); after
                                  Moon & Spencer (1942), b=2.0 (default), for sky luminance,
                                  useful in the tropics (Mitchell & Whitmore (1993).
                                  */
#define ZEN_MIN		30.0		/* minimum zenith angle for clumping index, degrees */
#define ZEN_MAX		60.0		/* maximum zenith angle for clumping index, degrees */

// int main(int argc, char *argv[])

// [[Rcpp::export]]

SEXP opennessR(SEXP fracR, SEXP slopeR, SEXP aspectR, SEXP imgId) 
{
        Rcpp::List fracC(fracR);
        Rcpp::NumericVector seg =
            Rcpp::as<Rcpp::NumericVector>(fracC["segments"]); //extract zen and azi
        Rcpp::NumericVector frac =
          Rcpp::as<Rcpp::NumericVector>(fracC["countPix"]); //extract fraction column
        CharacterVector imgIdC(imgId);
        //double open(List frac, )  
        //FILE *fp_prm, *fp_gap, *fpout;

        int i, j, k, compt, compt_plus, compt_gap, compt_log;
        int ZENITH, AZIMUTH;

        double sum_mean_log, sum_gap, mean_log, log_mean;
	  	  double openness, tot_azim, zen, ann_zen, az, gap;
        double opn_sky_can, sky_hem, theta1, theta2;
        double inclin, aspect, cosinclin, sininclin;
        double cosaltit, sinaltit, cosgam, coszen;
        double uoc_int=0.0, uoc_ext=0.0, soc_int=0.0, soc_ext=0.0;
        double lum, lum_uoc, lum_soc;
        double STP_ZN, STP_AZ, STRTD_ZN, STRT_AZ;

        // static double tableau[LIN+1][COL+1];
        /* angular and gap fraction data: the 1st and 2nd columns are mid-points
        of zenith and azimuth angles, respectively, the 3rd column is gap fraction
        data, the 4th column is number of pixels
        */

        static double fract_open[MAX_ZEN+1][MAX_AZI+1];
        /* work array for gap fraction data
        */
       // static double fract_sum[MAX_ZEN][DIM];
        /* work array for gap, CO, OS, CC:
                fract_sum[i][0] --> zenith angles, degrees
                fract_sum[i][1] --> average gap fraction
                fract_sum[i][2] --> fraction cumulated gap
                fract_sum[i][3] --> fraction cumulated CO
                fract_sum[i][4] --> fraction cumulated CC
                fract_sum[i][5] --> fraction cumulated UOC
                fract_sum[i][6] --> fraction cumulated SOC
        */
      NumericVector zenang(seg[0]), zengap(seg[0]), cumgap(seg[0]), ope(seg[0]), uoc(seg[0]), soc(seg[0]);

        /*------------------------ tests command line ----------------------*/
        //if (argc != 4)  {
          //      fprintf(stderr,
            //    "\n--->Usage: < %s > [parameter file][gap fraction file][output file]<---\n\n", argv[0]);
              //  exit(1);
        //}

        /*------------------------ opens parameter (configutation) file ----*/
        //if ((fp_prm=fopen(argv[1],"r+")) == NULL)        {
          //      fprintf(stderr,"Cannot open the file %s...\n",argv[1]);
            //    exit(1);
        //}

        //fscanf(fp_prm,"%lf",&inclin);         /* slope  */
        inclin = Rcpp::as<double>(slopeR);

        //fscanf(fp_prm,"%lf",&aspect);         /* aspect */
        aspect = Rcpp::as<double>(aspectR);


        /*------------------------ opens gap fraction data file ------------*/
        //if ((fp_gap=fopen(argv[2],"r+")) == NULL)     {
          //      fprintf(stderr,"Cannot open the file %s...\n",argv[2]);
            //    exit(1);
        //}

        /*------- adjusts for sky divisions --------------------------------*/
        //fscanf(fp_gap, "%d %d", &z0, &a0);
        
        ZENITH=seg[0];
        STP_ZN=90.0/(double)ZENITH;
        STRTD_ZN=STP_ZN/2.0;
        AZIMUTH=seg[1];
        STP_AZ=360.0/(double)AZIMUTH;
        STRT_AZ=STP_AZ/2.0;

        /* reading data --------------------------------------------------*/
        //for (i = 0; i < LIN; i++)
           // for (j = 0; j < COL; j++)
             //   fscanf(fp_gap,"%lf",&tableau[i][j]);

        for (compt_plus = 0, i = 0; i < ZENITH; i++)
            for (j = 0; j < AZIMUTH; j++)    {
                fract_open[i][j] = frac[compt_plus]; //tableau[compt_plus][2];
                compt_plus++;
            }

        /*----------- computations of gap, GT, CO, and CC -------------------*/
        inclin *= DEG_TO_RAD;
        aspect *= DEG_TO_RAD;
        cosinclin = cos(inclin);
        sininclin = sin(inclin);

	  	sum_gap = sum_mean_log = 0.0;
	  	compt_gap = compt_log = 0;
        openness = opn_sky_can = 0.0;
	 	 uoc_ext = soc_ext = uoc_int =soc_int = 0.0;
        for (zen = STRTD_ZN, k = 0, i = 0; i < ZENITH; i++)    {
			coszen = cos(zen * DEG_TO_RAD);
			cosaltit = cos(M_PI_2 - zen * DEG_TO_RAD);
			sinaltit = sin(M_PI_2 - zen * DEG_TO_RAD);
			lum = (1.0 + OVERCAST * sinaltit) / (1.0 + OVERCAST);
			lum_uoc = sinaltit * cosaltit;
			lum_soc = lum * lum_uoc;
			theta1 = (zen - (STP_ZN / 2.0)) * DEG_TO_RAD;
			theta2 = (zen + (STP_ZN / 2.0)) * DEG_TO_RAD;
			ann_zen = cos(theta1) - cos(theta2);
			sky_hem = 1.0 - cos(theta2);
			for (compt = AZIMUTH, az = STRT_AZ, tot_azim = 0.0, j = 0; j < AZIMUTH; j++)   {
					cosgam = cosinclin * sinaltit + sininclin * cosaltit
								* cos(az * DEG_TO_RAD - aspect);
					/* hemisphere hidden (orographic mask) */
					if (cosgam < 0.0)       {
						compt--;
						sky_hem -= (ann_zen / (double)AZIMUTH);
					}
					/* hemisphere open ------------------- */
					if (cosgam >= 0.0) {
						gap = fract_open[i][j];
						tot_azim += gap;
						openness += gap;
						opn_sky_can += gap * ann_zen / (double)AZIMUTH;
						uoc_ext += lum_uoc;
						soc_ext += lum_soc;
						uoc_int += gap * lum_uoc;
						soc_int += gap * lum_soc;
						k++;
						if (zen >= ZEN_MIN && zen <= ZEN_MAX)	{
							sum_gap += gap;
							compt_gap++;
						}
					}
					az += STP_AZ;
			}
			if (sum_gap)	{
				sum_mean_log += log(sum_gap/(double)compt_gap);
				compt_log++;
			}
    	zenang[i] = zen;                       /* zenith angle */
			zengap[i] = tot_azim / (double)compt;  /* zenithal gap profile */
			cumgap[i] = openness / (double)k;      /* cumulated gap */
			ope[i] = opn_sky_can / sky_hem;     /* canopy openness */
			uoc[i] = uoc_int / uoc_ext;         /* UOC */
			soc[i] = soc_int / soc_ext;         /* SOC */
      //fract_sum[i][0] = zen;                       /* zenith angle */
			//fract_sum[i][1] = tot_azim / (double)compt;  /* zenithal gap profile */
			//fract_sum[i][2] = openness / (double)k;      /* cumulated gap */
			//fract_sum[i][3] = opn_sky_can / sky_hem;     /* canopy openness */
			//fract_sum[i][4] = uoc_int / uoc_ext;         /* UOC */
			//fract_sum[i][5] = soc_int / soc_ext;         /* SOC */
			zen += STP_ZN;
		}
		mean_log = sum_mean_log/(double)compt_log;
		log_mean = log(sum_gap/(double)compt_gap);

        
        
        
        
        
        /*~~~~~~~ open output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
//        if ((fpout = fopen(argv[3], "a+")) == NULL)      {
//			fprintf(stderr, "Cannot open the file %s...\n", argv[3]);
//			exit(1);
//		}
//
//        fprintf(fpout,"< Program OPENNESS > --- %s ---\n", argv[3]);
//	  	fprintf(fpout, "CIMES (C) 1989-2009 Jean-Michel Walter\n\n");
//        fprintf(fpout,"Photosite: %s\n",argv[2]);
//        fprintf(fpout, "\nNbr of Zenith Rings = %d (Width = %4.1lf deg)", ZENITH, STP_ZN);
//        fprintf(fpout, "\nNbr of Azimuth Sectors = %d (Width = %4.1lf deg)\n",AZIMUTH, STP_AZ);
//	  	fprintf(fpout, "Clumping Index (30 deg --- 60 deg) = %5.3lf\n", mean_log/log_mean);
//        if (inclin)
//                fprintf(fpout,"\nSlope = %2.0lf deg   Aspect = %3.0lf deg\n",
//                               RAD_TO_DEG * inclin, RAD_TO_DEG * aspect);
//        else if (inclin == 0.0)
//                fprintf(fpout, "\nHorizontal terrain\n");
//        fprintf(fpout,"\nGap and Canopy Morphology --- Diffuse Sky Light\n");
//        fprintf(fpout,"------------------------------------------------\n");
//        fprintf(fpout,"Zenith Gap    Gap    Canopy Canopy Diff.  Diff.  \n");
//        fprintf(fpout,"Angle  Fract. Total  Open.  Closr. UOC    SOC    \n");
//        fprintf(fpout,"------------------------------------------------\n");
//        for (i= 0; i < ZENITH; i++)
//                fprintf(fpout,"%-7.1lf%-7.3lf%-7.3lf%-7.3lf%-7.3lf%-7.3lf%-7.3lf\n",
//                        fract_sum[i][0],
//                        fract_sum[i][1],
//                        fract_sum[i][2],
//                        fract_sum[i][3],
//                        1.0 - fract_sum[i][3],
//                        fract_sum[i][4],
//                        fract_sum[i][5]);
//        fprintf(fpout,"------------------------------------------------\n\n\n");

SEXP res = R_NilValue;
  res = Rcpp::List::create(Rcpp::Named("zenith_angle")  = zenang,
                           Rcpp::Named("gap_fract")  = zengap,
                           Rcpp::Named("gap_total")  = cumgap,
                           Rcpp::Named("canopy_open")  = ope,
                           Rcpp::Named("canopy_closr")  = (1.0 - ope),
                           Rcpp::Named("diff_UOC")  = uoc,
                           Rcpp::Named("diff_SOC")  = soc,
                           Rcpp::Named("zenith")  = ZENITH,
                           Rcpp::Named("azimuth")  = AZIMUTH,
                           Rcpp::Named("zenith_width")  = STP_ZN,
                           Rcpp::Named("azimuth_width")  = STP_AZ,
                           Rcpp::Named("clumping")  = (mean_log/log_mean),
                           Rcpp::Named("slope")  = inclin,
                           Rcpp::Named("aspect")  = aspect,
                           Rcpp::Named("imgId")  = imgIdC);

        //fclose(fp_prm);
        //fclose(fp_gap);
        //fclose(fpout);
		
		return res;
}
