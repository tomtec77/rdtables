/*
  axes.c
  
  Calculate the values of disc mass fraction, spin parameter and bulge mass
  fraction which label each table entry, or read them from files

  Link with rdtables.c
*/
#include "allvars.h"
#include "proto.h"
#include "readparameters.h"


void determine_axes_values(void)
{
  FILE *fp;
  int i, check;
  /*double a;*/

  
  for (i=0; i<NumMdiscBins; i++) {
    Mdisc_rscale[i] = pow(10.0, StartingMdisc+i*WidthOfMdiscBins);
  }
  printf("Mdisc range: [%g - %g]\n", Mdisc_rscale[0], 
	 Mdisc_rscale[NumMdiscBins-1]);
  if (Mdisc_rscale[NumMdiscBins-1] > 1.25*BaryonFrac) {
    printf("#############################################################\n");
    printf("# WARNING: maximum md exceeds 1.25 times the univ. f_baryon #\n");
    printf("#############################################################\n");
    check = validate_continuation();
    if (check == 0) {
      printf("Aborted on user input.\n");
      exit(0);
    }
    printf("\n");
  }
#ifdef DEBUG
  printf("\nDisc mass fractions:\n");
  for (i=0; i<NumMdiscBins; i++) {
    printf("  %.6f\n", Mdisc_rscale[i]);
  }
  printf("\n");
#endif


  for (i=0; i<NumLambdaBins; i++)
    Lambda_rscale[i] = pow(10.0, StartingLambda+i*WidthOfLambdaBins);
  printf("Spin parameter range: [%g - %g]\n", Lambda_rscale[0], 
	 Lambda_rscale[NumLambdaBins-1]);
  if (Lambda_rscale[NumLambdaBins-1] > LAMBDASAFETY) {
    printf("#############################################################\n");
    printf("# WARNING: values of lambda > %d are very rare           #\n",
	   LAMBDASAFETY);
    printf("#############################################################\n");
    check = validate_continuation();
    if (check == 0) {
      printf("Aborted on user input.\n");
      exit(0);
    }
    printf("\n");
  }
#ifdef DEBUG
  printf("\nSpin parameters:\n");
  for (i=0; i<NumLambdaBins; i++) {
    printf("  %.6f\n", Lambda_rscale[i]);
  }
  printf("\n");
#endif

  if (MbulgeListOn == 1) { /* Read bulge mass values from file */
    printf("Reading table bulge mass values from file '%s'...\n",
	   NameMbulgeList);
    if (!(fp=fopen(NameMbulgeList,"r"))) {
      fprintf(stderr, "Error: cannot open file '%s': %s (%u)\n", 
	      NameMbulgeList, strerror(errno), errno);
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
    for (i=0; i<NumMbulgeBins; i++) {
      fscanf(fp, " %lf\n", &Mbulge_rscale[i]);
    }
  }
  else {
    Mbulge_rscale[0] = StartingMbulge;
    
    for (i=1; i<NumMbulgeBins; i++)
      Mbulge_rscale[i] = Mbulge_rscale[i-1]+WidthOfMbulgeBins;
  }
  printf("Mbulge range: [%g - %g]\n", Mbulge_rscale[0], 
	 Mbulge_rscale[NumMbulgeBins-1]);
#ifdef DEBUG
  printf("\nBulge mass fractions:\n");
  for (i=0; i<NumMbulgeBins; i++) {
    printf("  %.6f\n", Mbulge_rscale[i]);
  }
  printf("\n");
#endif
    
  return;
}

