/**
 * @file generate_rscale_table.c
 * @brief Create tables of scale radius and circular velocity
 * @author Tomas E. Tecce
 */
#include "allvars.h"
#include "proto.h"
#include "readparameters.h"


int generate_rscale_table(int isnap, double c200, char *rdfile, 
			   char *vcfile)
{
  int status=0;
  FILE *fp, *fv;
  int i, j, k, total;
  unsigned int offset;
  double buf, zcurr;
  int invalidFr, invalidVc;
  char headerfile[FILENAME_MAX];

  
  if (RdTablesOn == 0 && VcTablesOn == 0) {
    fprintf(stderr, "Error (generate_rscale_table): no tables will be "
	    "generated!\n");
    fprintf(stderr, "Check parameter file and try again - Exit\n");
    fflush(stderr);
    return -1;
  }
  
  zcurr = ZZ[isnap];

  /* Table header: redshift, concentration, number of rows and columns and
     number of data sets (values of mbulge), then label values of rows,
     then label values of columns */
  if (isnap == 0) {
    sprintf(headerfile, "%s/tables_header.dat", Path1);
    fp = fopen(headerfile, "w");
    if (fp==NULL) {
      fprintf(stderr, "Error (generate_rscale_table): "
	      "cannot open output file '%s': %s (%u)\n", 
	      headerfile, strerror(errno), errno); fflush(stdout);
      return -1;
    }

    printf("Creating header file %s/tables_header.dat... ", Path1);
    fprintf(fp, "%d %d %d %d %d\n", NumOutputTimes, NumCbins, NumMdiscBins, 
            NumLambdaBins, NumMbulgeBins);
    for (i=0; i<NumOutputTimes; i++)
      fprintf(fp, "%le ", ZZ[i]);
    fprintf(fp, "\n");
    for (i=0; i<NumMdiscBins; i++)
      fprintf(fp, "%le ", Mdisc_rscale[i]);
    fprintf(fp, "\n");
    for (i=0; i<NumLambdaBins; i++)
      fprintf(fp, "%le ", Lambda_rscale[i]);
    fprintf(fp, "\n");
    for (i=0; i<NumMbulgeBins; i++)
      fprintf(fp, "%le ", Mbulge_rscale[i]);
    fprintf(fp, "\n");

    fclose(fp);
    printf("done.\n");
  }
  
  if (RdTablesOn != 0) {
    fp = fopen(rdfile, "w");
    if (fp==NULL) {
      fprintf(stderr, "Error (generate_rscale_table): "
	      "cannot open output file '%s': %s (%u)\n", 
	      rdfile, strerror(errno), errno); fflush(stderr);
      return -1;
    }
    fprintf(fp, "%le %le\n", zcurr, c200);
  }

  if (VcTablesOn != 0) {
    fv = fopen(vcfile, "w");
    if (fv==NULL) {
      fprintf(stderr, "Error (generate_rscale_table): "
	      "cannot open output file '%s': %s (%u)\n", 
	      vcfile, strerror(errno), errno); fflush(stdout);
      return -1;
    }
    fprintf(fv, "%le %le\n", zcurr, c200);
  }

  invalidFr = invalidVc = 0;
  total = NumMbulgeBins*NumLambdaBins*NumMdiscBins;
  for (k=0; k<NumMbulgeBins; k++) {

    if (RdTablesOn != 0) fprintf(fp, "%le\n", Mbulge_rscale[k]);
    if (VcTablesOn != 0) fprintf(fv, "%le\n", Mbulge_rscale[k]);

    for (i=0; i<NumMdiscBins; i++) {
      for (j=0; j<NumLambdaBins; j++) {
	
	offset = indexrm(i,j,NumMdiscBins,NumLambdaBins);

	rditer(c200, Mdisc_rscale[i], Lambda_rscale[j], Mbulge_rscale[k],
	       zcurr, 1.0, &Frtable[offset], &Vctable[offset], &buf);
	printf("%5.1f %% \r", 
	       (float)100*(j + NumLambdaBins*(i + NumMdiscBins*k))/total); 
	fflush(stdout);

	/* Print data to file */
	if (RdTablesOn != 0) fprintf(fp, "%le ", Frtable[offset]);
	if (VcTablesOn != 0) {
	  fprintf(fv, "%le ", Vctable[offset]);
	  if (Vctable[offset] == 0.0) {
	    fprintf(stderr, "Error: found Vc = 0!\n");
	    fprintf(stderr, "i=%d, j=%d, rd=%g, vc=%g\n\n", i, j, 
		    Frtable[offset], Vctable[offset]);
	    exit(EXIT_FAILURE);
	  }
	}

	/* Data check */
	if (RdTablesOn != 0) {
	  if (Frtable[offset] == RDUNDEFINED) 
	    invalidFr++;
	}
	if (VcTablesOn != 0) {
	  if (Vctable[offset] == RDUNDEFINED) 
	    invalidVc++;
	}
      }
      if (RdTablesOn != 0) fprintf(fp, "\n");
      if (VcTablesOn != 0) fprintf(fv, "\n");
    }
  }
  printf("100.0 %%           \n"); fflush(stdout);

  if (RdTablesOn != 0) {
    if (invalidFr > 0)
      printf("Warning: %d undefined values in table '%s'\n", 
	     invalidFr, rdfile);
  }

  if (VcTablesOn != 0) {
    if (invalidVc > 0)
      printf("Warning: %d undefined values in table '%s'\n", 
	     invalidVc, vcfile);
  }

  if (RdTablesOn != 0) fclose(fp);
  if (VcTablesOn != 0) fclose(fv);

  return status;
}
