#include "allvars.h"

double *Aexp, *ZZ;

double *Mbulge_rscale, *Conc_rscale, *Mdisc_rscale;
double *Lambda_rscale;

double *Cmin_conc, *Cmax_conc;
double *Time_conc, *Mass_conc;
double **Table_conc;

int NumTableElements;

double *Xgl, *Wgl;

double MaxCvirRedshift;
unsigned int Nmd, Nc, Nl, Num;
//unsigned int Snapshot, StartSnapshot, MaxSnapshot;
unsigned int OutputListOn;
double Mvir[NMASS]= {1e10, 1e11, 1e12, 1e13, 1e14, 1e15};
double Mbulge[NMBULGE] = {0.0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15};

double *Mdisc, *Lambda, *Frtable, *Vctable, *Rotable;

double TimeOfFirstSnapshot, TimeBetSnapshot;

char path[FILENAME_MAX];  /* Path to output tables */
char path3[FILENAME_MAX];

unsigned int Count_div, Count_itmax;


