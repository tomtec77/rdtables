/**
 * @file output_times.c
 * @brief Determine the redshift of the disc scale length tables
 * @author Tomas E. Tecce
 */
#include "allvars.h"
#include "proto.h"
#include "readparameters.h"

void output_times(void)
{
  FILE *fp;
  int snap;
  double aend, a;
  unsigned int check;
#ifdef DEBUG
  int i;
#endif


  if (OutputListOn == 1) {    /* Reading output times from external file */

    if (NumOutputTimes >= MAXTIMEOUTS) {
      fprintf(stderr,"Error: number of time outputs exceeded - Exit\n");
      fflush(stderr);
      exit(EXIT_FAILURE);
    }
    
    printf("Reading output times from file '%s'...\n",
	   NameOutputList);
    if (!(fp = fopen(NameOutputList, "r"))) {
      fprintf(stderr, "Error: cannot open file '%s': %s (%u)\n",
	      NameOutputList, strerror(errno), errno);
      fflush(stderr);
      exit(EXIT_FAILURE);
    }

    for (snap=0; snap<NumOutputTimes; snap++) {
      fscanf(fp, " %lf ", &Aexp[snap]);
      
      ZZ[snap] = 1.0/Aexp[snap] - 1;
    }

    printf("Tables will be calculated from z = %g (a = %g) to "
	   "z = %g (a = %g)\n\n", ZZ[0], Aexp[0], ZZ[NumOutputTimes-1], 
	   Aexp[NumOutputTimes-1]);
  }
  else {
    /* Alternative: give a target final value for the expansion factor a,
       the interval between outputs in log a and the number of total output 
       times desired */
    aend = TimeOfFinalOutput/pow(TimeBetOutputs,NumOutputTimes);

    for (snap=NumOutputTimes-1, a=TimeOfFinalOutput; snap>=0; 
	 snap--, a/=TimeBetOutputs) {
      Aexp[snap] = a;
      ZZ[snap] = 1.0/a -1;
    }    

    printf("With the parameters selected, tables will be calculated from\n");
    printf("z = %g (a = %g) to z = %g (a = %g)\n\n", ZZ[0], Aexp[0], 
	   ZZ[NumOutputTimes-1], Aexp[NumOutputTimes-1]);
  }

  if (ZZ[0] > CVIRZSAFETY) {
    printf("#############################################################\n");
    printf("# WARNING: program will probably fail for z > %2d            #\n",
	   CVIRZSAFETY);
    printf("#############################################################\n");
    check = validate_continuation();
    if (check == 0) {
      printf("Aborted on user input.\n");
      exit(0);
    }
    printf("\n");
  }

#ifdef DEBUG
  printf("\nOutput times:\n");
  printf("  Aexp        Z\n");
  for (i=0; i<NumOutputTimes; i++) {
    printf("  %.6f    %.6f\n", Aexp[i], ZZ[i]);
  }
  printf("\n");
#endif
  
  return;
}
