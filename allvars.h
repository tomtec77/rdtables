#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>

#include "lib/nrsrc/nrutil.h"
#include "lib/nrsrc/nrsag.h"
#include "lib/myutils/myutils.h"

/** Solar mass in g */
#define SOLARMASSG 1.98892e33 

/** Hubble constant in km/s/Mpc */
#define HUBBLE 100

/** Gravitational constant in Mpc*(km/s)^2/Msun */
#define GGRAV 4.3e-9 

/** 1 Mpc in cm */
#define MPCINCM 3.08568025e24 

/** Maximum number of redshift outputs */
#define MAXTIMEOUTS 200 

/** Upper limit for redshift */
#define CVIRZSAFETY 30

/** Upper limit for spin parameter */
#define LAMBDASAFETY 2 

/** Upper limit for c200 */
#define MINC200 0

/** Lower limit for c200 */
#define MAXC200 50

/** Log min/max disc mass in Msun */
#define MINDISCMASS 6.0  

#define NPOINTSGL 4

#define OUTPUTS 102
#define NMASS 6
#define NMBULGE 9

#define MAXDISCMASS 13.0
#define MINMD -6
#define MAXMD -0.6
#define MDSTEP 0.1

/** Log of minimum spin parameter */
#define SPINMIN -2.3

/** Log of maximum spin parameter */
#define SPINMAX 0.2     

#define SPINSTEP 0.1 

//#define STEPS 10000          /* Number of points in r to calculate Vc */

#define ZBRAC_FACTOR 1.6
#define ZBRAC_NTRY 50

#define RDUNDEFINED -99        


extern double *Aexp, *ZZ;

extern double *Mbulge_rscale, *Conc_rscale, *Mdisc_rscale;
extern double *Lambda_rscale;

extern double *Cmin_conc, *Cmax_conc;
extern double *Time_conc, *Mass_conc;
extern double **Table_conc;

extern int NumTableElements;

extern double *Xgl, *Wgl;

extern double MaxCvirRedshift;
extern unsigned int Nmd, Nc, Nl, Num;
extern unsigned int OutputListOn;

extern double TimeOfFirstSnapshot, TimeBetSnapshot;

extern double Mvir[NMASS];
extern double Mbulge[NMBULGE];
extern double *Mdisc, *Lambda, *Frtable, *Vctable, *Rotable;

/** Path to output tables */
extern char path[FILENAME_MAX]; 

/** Path to file with output times */
extern char path3[FILENAME_MAX];

extern unsigned int Count_div, Count_itmax;
