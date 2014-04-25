/**
 * @file init.c
 * @brief Initialisation functions
 * @author Tomas E. Tecce
 */
#include "allvars.h"
#include "proto.h"
#include "readparameters.h"


void allocate_arrays(void)
{
  /*FILE *fp;*/
  int bytes;
  

  if (!(Aexp = malloc(bytes=NumOutputTimes*sizeof(double)))) {
    fprintf(stderr, "Error: failed to allocate memory for array Aexp: "
	    "%s (%u)\n", strerror(errno), errno); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  printf("Allocated memory %f MB for Aexp\n", 
	 ((float)bytes)/(1024.0*1024.0)); fflush(stdout);

  if (!(ZZ = malloc(bytes=NumOutputTimes*sizeof(double)))) {
    fprintf(stderr, "Error: failed to allocate memory for array ZZ: "
	    "%s (%u)\n", strerror(errno), errno); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  printf("Allocated memory %f MB for ZZ\n", 
	 ((float)bytes)/(1024.0*1024.0)); fflush(stdout);

  if (!(Mbulge_rscale = malloc(bytes=NumMbulgeBins*sizeof(double)))) {
    fprintf(stderr, "Error: failed to allocate memory for array "
	    "Mbulge_rscale: %s (%u)\n", strerror(errno), errno);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  printf("Allocated memory %f MB for Mbulge\n", 
	 ((float)bytes)/(1024.0*1024.0)); fflush(stdout);

  if (!(Conc_rscale = malloc(bytes=NumCbins*sizeof(double)))) {
    fprintf(stderr, "Error: failed to allocate memory for array "
	    "Conc_rscale: %s (%u)\n", strerror(errno), errno);  
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  printf("Allocated memory %f MB for Conc_rscale\n", 
	 ((float)bytes)/(1024.0*1024.0)); fflush(stdout);

  if (!(Mdisc_rscale = malloc(bytes=NumMdiscBins*sizeof(double)))) {
    fprintf(stderr, "Error: failed to allocate memory for array "
	    "Mdisc_rscale: %s (%u)\n", strerror(errno), errno);  
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  printf("Allocated memory %f MB for Mdisc_rscale\n", 
	 ((float)bytes)/(1024.0*1024.0)); fflush(stdout);

  if (!(Lambda_rscale = malloc(bytes=NumLambdaBins*sizeof(double)))) {
    fprintf(stderr, "Error: failed to allocate memory for array "
	    "Lambda_rscale: %s (%u)\n", strerror(errno), errno); 
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  printf("Allocated memory %f MB for Lambda_rscale\n", 
	 ((float)bytes)/(1024.0*1024.0)); fflush(stdout);
  
  NumTableElements = NumMdiscBins*NumLambdaBins;
  if (!(Frtable = malloc(bytes=NumTableElements*sizeof(double)))) {
    fprintf(stderr, "Error: failed to allocate memory for Frtable: "
	    "%s (%u)\n", strerror(errno), errno); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  printf("Allocated memory %f MB for Frtable\n\n", 
	 ((float)bytes)/(1024.0*1024.0)); fflush(stdout);

  /*if (!(Rotable = malloc(bytes=NumTableElements*sizeof(double)))) {
    fprintf(stderr, "Error: failed to allocate memory for Rotable: "
	    "%s (%u)\n", strerror(errno), errno); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  printf("Allocated memory %f MB for Rotable\n\n", 
  ((float)bytes)/(1024.0*1024.0)); fflush(stdout);*/

  if (!(Vctable = malloc(bytes=NumTableElements*sizeof(double)))) {
    fprintf(stderr, "Error: failed to allocate memory for Vctable: "
	    "%s (%u)\n", strerror(errno), errno); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  printf("Allocated memory %f MB for Vctable\n\n", 
	 ((float)bytes)/(1024.0*1024.0)); fflush(stdout);

  /* Abscissas and weights for Gauss-Laguerre integration */
  Xgl = dvector(1, NPOINTSGL);
  Wgl = dvector(1, NPOINTSGL);

  return;
}


void free_arrays(void)
{
  free(Aexp);
  free(ZZ);
  free(Mdisc_rscale);
  free(Lambda_rscale);
  free(Mbulge_rscale);
  free(Frtable);
  free(Vctable);
  /*free(Rotable);*/

  free_dvector(Xgl, 1, NPOINTSGL);
  free_dvector(Wgl, 1, NPOINTSGL);

  return;
}
