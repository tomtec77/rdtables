#include "allvars.h"
#include "proto.h"
#include "readparameters.h"

/*
 * Load file with general information on the rscale tables, and allocate
 * memory for the tables
 *    char *pathname:    Path to table files
 */
void load_tables_header(char *pathname)
{
  FILE *fp;
  int i, j;
  char filename[FILENAME_MAX];
  double zcurr, c200, mbulge;
  
  printf("Reading tables header file from directory '%s'...\n",
	 ntables, pathname);

  Time_rscale = dvector(1, ntime);
  Conc_rscale = dvector(1, nc);

  sprintf(filename, "%s/tables_header.dat", pathname);
  fp = fopen(filename, "r");
  if (fp==NULL) {
    fprintf(stderr, "Error: cannot open table header file '%s': '%s' "
	    "(%u)\n", filename, strerror(errno), errno); fflush(stderr);
    exit(EXIT_FAILURE);
  }
      
  fscanf(fp, "%d %d %d %d %d\n", &NumOutputTimes, &NumCbins, &NumMdiscBins, 
	 &NumLambdaBins, &NumMbulgeBins);
  
  /* Allocate memory for the table axes */
  Time_rscale   = dvector(1, NumOutputTimes);
  Conc_rscale   = dvector(1, NumCbins);
  Mdisc_rscale  = dvector(1, NumMdiscBins);
  Lambda_rscale = dvector(1, NumLambdaBins);
  Mbulge_rscale = dvector(1, NumMbulgeBins);

  for (k=1; k<=NumOutputTimes; k++)
    fscanf(fp, "%le ", &Time_rscale[k]);
  for (k=1; k<=NumCbins; k++)
    fscanf(fp, "%le ", &Conc_rscale[k]);
  for (k=1; k<=NumMdiscBins; k++)
    fscanf(fp, "%le ", &Mdisc_rscale[k]);
  for (k=1; k<=NumLambdaBins; k++)
    fscanf(fp, "%le ", &Lambda_rscale[k]);
  for (k=1; k<=NumMbulgeBins; k++)
    fscanf(fp, "%le ", &Mbulge_rscale[k]);

  return;
}

	 
/*
 * For a given redshift, load onto memory the two sets of tables whose
 * redshifts bracket the selected time. Function load_tables_header must be
 * called first to initialise values.
 *     double zcurr:    Selected output time
 *     char *pathname:  Path to table files
 */
void load_rscale_tables(double zcurr, char *pathname)
{
  FILE *fp;
  int i, j, k, p, q;
  unsigned long itime;
  char filename[FILENAME_MAX];
  double zcurr, c200, mbulge;


  /* Locate the corresponding outputs in time */
  dlocate(Time_rscale, NumOutputTimes, zcurr, &itime);
  if (itime == 0 || itime == NumOutputTimes) {
    fprintf(stderr, "Error (load_rscale_tables): redshift out of range - "
	    "Exit\n"); fflush(stdout);
    exit(EXIT_FAILURE);
  }

  /* Each tensor stores an entire set of data files for a given
     redshift. The first dimension labels concentrations, the second values of
     Mbulge and the data table for those values is stored linearly in the
     third dimension, row-major ordered. */
  Table0_rscale = d3tensor(1, NumCbins, 1, NumMbulgeBins, 
			   1, NumMdiscBins*NumLambdaBins);
  Table1_rscale = d3tensor(1, NumCbins, 1, NumMbulgeBins, 
			   1, NumMdiscBins*NumLambdaBins);
  

  for (k=0; k<2; k++) {
    for (i=0; i<NumCbins; i++) {
      if (k == 0) 
	sprintf(filename, "%s/rdtable_%03d_%03d.dat", 
		pathname, itime-1, i);
      else
	sprintf(filename, "%s/rdtable_%03d_%03d.dat", 
		pathname, itime, i);
    
      if (!(fp = fopen(filename, "r"))) {
	fprintf(stderr, "Error: cannot open table file '%s': %s (%u)\n",
		filename0, strerror(errno), errno);
	fflush(stderr);
	exit(EXIT_FAILURE);
      }

      fscanf(fp, "%le %le\n", zcurr, c200);
      if (zcurr != Time_rscale[itime] || c200 != Conc_rscale[i]) {
	fprintf(stderr, "Error (load_rscale_tables): read mismatch (1) - "
		"Exit\n"); fflush(stderr);
	exit(EXIT_FAILURE);
      }

      for (j=1; j<=NumMbulgeBins; j++) {
	fscanf(fp, "%le\n", &mbulge);
	if (mbulge != Mbulge_rscale[j]) {
	  fprintf(stderr, "Error (load_rscale_tables): read mismatch (2) - "
		  "Exit\n"); fflush(stderr);
	  exit(EXIT_FAILURE);
	} 

	for (p=0; p<NumMdiscBins; p++) {
	  for (q=0; q<NumLambdaBins; q++) {
	    offset = indexrm(p, q, NumMdiscBins, NumLambdaBins);
	    
	    if (k==0)
	      fscanf(fp, "%le ", &Table0_rscale[i][j][offset]);
	    else
	      fscanf(fp, "%le ", &Table1_rscale[i][j][offset]);
	  }
	}
      }
    }
  }
	 
  return;
}


/*
 * Calculate the scale radius interpolating on the tables.
 * For the given values, find the corresponding values in the two tables
 * that bracket the current redshift, and perform a simple linear
 * interpolation.
 */
double get_scale_radius(double chalo, double md, double lambda, double mb)
{
  unsigned long ic, imb, imd, il;
  double 

  /* Find the indices that enclose the selected values */
  dlocate(Conc_rscale, NumCbins, chalo, &ic);
  dlocate(Mbulge_rscale, NumMbulgeBins, mb, &imb);
  dlocate(Mdisc_rscale, NumMdiscBins, md, &imd);
  dlocate(Lambda_rscale, NumLambdaBins, lambda, &il);

  

  

  return;
}
