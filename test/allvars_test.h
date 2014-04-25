/**
 * @file allvars_test.h
 * @brief Header file with includes, definitions and global variables for 
 * use with rdtest.c
 * @author Tomas E. Tecce
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#include "../lib/nrsrc/nr.h"
#include "../lib/nrsrc/nrsag.h"
#include "../lib/nrsrc/nrutil.h"

/** Nr. of points for polynomial interpolation */
#define POLINT_POINTS 4    


extern int NumOutputTimes;
extern int NumCbins;
extern int NumOutputTimes;
extern int NumMdiscBins;
extern int NumLambdaBins;
extern int NumMbulgeBins;

extern double *Time_rscale, *Conc0_rscale, *Conc1_rscale;
extern double *Mdisc_rscale, *Lambda_rscale, *Mbulge_rscale;
extern double ****Table0_rscale, ****Table1_rscale;
extern double ****VcTable0_rscale, ****VcTable1_rscale;
