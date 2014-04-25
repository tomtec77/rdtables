/**
 * @file read_rscale_table.c
 * @brief Functions to load and read the new tables of disc scale radius
 * and circular velocity
 * @author Tomas E. Tecce
 */
#include "allvars_test.h"
#include "proto_test.h"


/**
 * @brief Load file with general information on the rscale tables, and 
 * allocate memory for the tables
 * @param *pathname Path to table files
 */
void load_tables_header(char *pathname)
{
  FILE *fp;
  int k;
  char filename[FILENAME_MAX];

  
  printf("Reading tables header file from directory '%s'...\n", pathname);

  sprintf(filename, "%s/tables_header.dat", pathname);
  fp = fopen(filename, "r");
  if (fp==NULL) {
    fprintf(stderr, "Error (load_tables_header): "
	    "cannot open table header file '%s': '%s' (%u)\n", 
	    filename, strerror(errno), errno); fflush(stderr);
    exit(EXIT_FAILURE);
  }
      
  fscanf(fp, "%d %d %d %d %d\n", &NumOutputTimes, &NumCbins, &NumMdiscBins, 
	 &NumLambdaBins, &NumMbulgeBins);
  
  /* Allocate memory for the table axes */
  /* The values stored in these vectors are the same for all tables */
  Time_rscale   = dvector(1, NumOutputTimes);
  Mdisc_rscale  = dvector(1, NumMdiscBins);
  Lambda_rscale = dvector(1, NumLambdaBins);
  Mbulge_rscale = dvector(1, NumMbulgeBins);

  for (k=1; k<=NumOutputTimes; k++)
    fscanf(fp, "%le ", &Time_rscale[k]);
  for (k=1; k<=NumMdiscBins; k++)
    fscanf(fp, "%le ", &Mdisc_rscale[k]);
  for (k=1; k<=NumLambdaBins; k++)
    fscanf(fp, "%le ", &Lambda_rscale[k]);
  for (k=1; k<=NumMbulgeBins; k++)
    fscanf(fp, "%le ", &Mbulge_rscale[k]);

  /* Allocate memory for the concentration tables selected */
  /* Each tensor stores an entire set of data files for a given
     redshift. The first dimension labels concentrations, the second values 
     of Mbulge and the data table for those values is stored in the
     third (Mdisc) and fourth (Lambda) dimensions. */
  Conc0_rscale    = dvector(1, NumCbins);
  Conc1_rscale    = dvector(1, NumCbins);
  Table0_rscale   = d4tensor(1,NumCbins, 1,NumMbulgeBins, 
			     1,NumMdiscBins, 1,NumLambdaBins);
  Table1_rscale   = d4tensor(1,NumCbins, 1,NumMbulgeBins, 
			     1,NumMdiscBins, 1,NumLambdaBins);
  VcTable0_rscale = d4tensor(1,NumCbins, 1,NumMbulgeBins, 
			     1,NumMdiscBins, 1,NumLambdaBins);
  VcTable1_rscale = d4tensor(1,NumCbins, 1,NumMbulgeBins, 
			     1,NumMdiscBins, 1,NumLambdaBins);

  printf("%d output times, %d concentration values\n\n", NumOutputTimes, 
	 NumCbins);
  
  return;
}

	 
/**
 * @brief Load into memory the two sets of tables whose redshifts bracket 
 * the selected time. 
 * @param zcurr Selected output time
 * @param *pathrscale Path to scale radius table files
 * @param *pathvc Path to circular velocity table files
 *
 * Function @a load_tables_header must be called first to initialise values
 */
void load_rscale_tables(double zcurr, char *pathrscale, char *pathvc)
{
  FILE *fp;
  int i, j, k, l, p;
  unsigned long itime;
  char filename[FILENAME_MAX];
  double zread, c200, mbulge;


  /* Locate the corresponding outputs in time */
  dlocate(Time_rscale, NumOutputTimes, zcurr, &itime);
  if (itime == 0 || itime == NumOutputTimes) {
    fprintf(stderr, "Error (load_rscale_tables): redshift out of range "
	    "[%f - %f] (%d) - Exit\n", Time_rscale[1], 
	    Time_rscale[NumOutputTimes], (int)itime); fflush(stdout);
    exit(EXIT_FAILURE);
  }

  for (p=0; p<2; p++) {
    for (i=0; i<NumCbins; i++) {
      /* Here we use itime-1 and itime because tables are numbered 
	 from 0 to N-1 */
      if (p == 0) {
	sprintf(filename, "%s/rdtable_%03d_%03d", 
		pathrscale, (int)(itime-1), i);
	if (i==NumCbins-1)
	  printf("Loading scale radius tables '%s/rdtable_%03d_000' to "
		 "'%s/rdtable_%03d_%03d' for z = %g\n", 
		 pathrscale, (int)(itime-1), pathrscale, (int)(itime-1), 
		 NumCbins-1, Time_rscale[itime]);
      }
      else {
	sprintf(filename, "%s/rdtable_%03d_%03d", 
		pathrscale, (int)itime, i);
	if (i==NumCbins-1)
	  printf("Loading scale radius tables '%s/rdtable_%03d_000' to "
		 "'%s/rdtable_%03d_%03d' for z = %g\n", 
		 pathrscale, (int)itime, pathrscale, (int)itime, NumCbins-1,
		 Time_rscale[itime+1]);
      }
    
      if (!(fp = fopen(filename, "r"))) {
	fprintf(stderr, "Error: cannot open table file '%s': %s (%u)\n",
		filename, strerror(errno), errno);
	fflush(stderr);
	exit(EXIT_FAILURE);
      }
      
      fscanf(fp, "%le %le\n", &zread, &c200);
      /*if (zread != Time_rscale[itime+p]) {
	fprintf(stderr, "Error: time mismatch when reading table '%s': "
		"z(table) = %f, z(outlist) = %f - Exit\n",
		filename, zread, Time_rscale[itime+p]); fflush(stdout);
	exit(EXIT_FAILURE);
	}*/
      if (p == 0)
	Conc0_rscale[i+1] = c200;
      else
	Conc1_rscale[i+1] = c200;

      for (j=1; j<=NumMbulgeBins; j++) {
	fscanf(fp, "%le\n", &mbulge);
	if (mbulge != Mbulge_rscale[j]) {
	  fprintf(stderr, 
		  "Error (load_rscale_tables, r): read mismatch (mbulge)"
		  " j=%d, %g, %g -  Exit\n", j, mbulge, Mbulge_rscale[j]); 
	  fflush(stderr);
	  exit(EXIT_FAILURE);
	} 

	for (k=1; k<=NumMdiscBins; k++) {
	  for (l=1; l<=NumLambdaBins; l++) {
	    if (p==0)
	      fscanf(fp, "%le ", &Table0_rscale[i+1][j][k][l]);
	    else
	      fscanf(fp, "%le ", &Table1_rscale[i+1][j][k][l]);
	  }
	}
      }
      fclose(fp);
    }
  }

  printf("z = %g concentration range: [%g - %g]\n", Time_rscale[itime],
	 Conc0_rscale[1], Conc0_rscale[NumCbins]);
  printf("z = %g concentration range: [%g - %g]\n", Time_rscale[itime+1],
	 Conc1_rscale[1], Conc1_rscale[NumCbins]);
  printf("\n");


  /*
   * Now load circular velocity tables
   */
  for (p=0; p<2; p++) {
    for (i=0; i<NumCbins; i++) {
      /* Here we use itime-1 and itime because tables are numbered 
	 from 0 to N-1 */
      if (p == 0) {
	sprintf(filename, "%s/vctable_%03d_%03d", 
		pathvc, (int)(itime-1), i);
	if (i==NumCbins-1)
	  printf("Loading circular velocity tables '%s/vctable_%03d_000' to"
		 " '%s/vctable_%03d_%03d' for z = %g\n", 
		 pathvc, (int)(itime-1), pathvc, (int)(itime-1), NumCbins-1,
		 Time_rscale[itime]);
      }
      else {
	sprintf(filename, "%s/vctable_%03d_%03d", pathvc, (int)itime, i);
	if (i==NumCbins-1)
	  printf("Loading circular velocity tables '%s/vctable_%03d_000' to "
		 "'%s/vctable_%03d_%03d' for z = %g\n", 
		 pathvc, (int)itime, pathvc, (int)itime, NumCbins-1,
		 Time_rscale[itime+1]);
      }
    
      if (!(fp = fopen(filename, "r"))) {
	fprintf(stderr, "Error (load_rscale_tables): "
		"cannot open table file '%s': %s (%u)\n",
		filename, strerror(errno), errno);
	fflush(stderr);
	exit(EXIT_FAILURE);
      }
      
      fscanf(fp, "%le %le\n", &zread, &c200);
      /*if (zread != Time_rscale[itime+p]) {
	fprintf(stderr, "Error: time mismatch when reading table '%s': "
		"z(table) = %f, z(outlist) = %f - Exit\n",
		filename, zread, Time_rscale[itime+p]); fflush(stdout);
	exit(EXIT_FAILURE);
	}*/

      for (j=1; j<=NumMbulgeBins; j++) {
	fscanf(fp, "%le\n", &mbulge);
	//printf("mbulge=%g\n", mbulge);
	if (mbulge != Mbulge_rscale[j]) {
	  fprintf(stderr, 
		  "Error (load_rscale_tables, v): read mismatch (mbulge)"
		  "j=%d, %g, %g -  Exit\n", j, mbulge, Mbulge_rscale[j]); 
	  fflush(stderr);
	  exit(EXIT_FAILURE);
	} 

	for (k=1; k<=NumMdiscBins; k++) {
	  for (l=1; l<=NumLambdaBins; l++) {
	    if (p==0)
	      fscanf(fp, "%le ", &VcTable0_rscale[i+1][j][k][l]);
	    else
	      fscanf(fp, "%le ", &VcTable1_rscale[i+1][j][k][l]);
	  }
	}
      }
      fclose(fp);
    }
  }
  printf("\n");

  return;
}


/**
 * @brief Calculate the scale radius or circular velocity interpolating on 
 * the corresponding tables
 * @param ****t0 First data table for the current redshift
 * @param ****t1 Second data table for the current redshift
 * @param zcurr Current redshift
 * @param chalo Host halo concentration @f$c_{200}@f$
 * @param mb Bulge mass fraction @f$M_b / M_{200}@f$
 * @param md Disc mass fraction @f$M_d / M_{200}@f$
 * @param lambda Spin parameter of the host halo
 *
 * For the given values, find the corresponding values in the two
 * tables that bracket the current redshift and interpolate to obtain the
 * result. Note that the value returned by this function is in the units of
 * the tables; it has to be multiplied by (virial mass/critical
 * density)^1/3 in the case of rscale and by the virial velocity in the
 * case of circular velocity to obtain a value in physical units.
 */
double get_data_from_table(double ****data0, double ****data1, 
			   double zcurr, double chalo, double mb, 
			   double md, double lambda)
{
  int i, j, kmd, kl;
  unsigned long ic[2], imb[2], iz[2];
  unsigned long imd, il;
  double **table;
  double ymb0[2][2], ymb1[2][2];
  double yc0[2], yc1[2];
  double yz0, yz1, buf;
#ifdef DEBUG
  int p, q;
#endif


  /* TOMAS 2012-12-06: patch to fix values that fall outside the range
     of the tables */
  if (chalo < Conc0_rscale[1]) chalo = Conc0_rscale[1];
  if (chalo > Conc0_rscale[NumCbins]) chalo = Conc0_rscale[NumCbins];
  
  if (mb < Mbulge_rscale[1]) mb = Mbulge_rscale[1];
  if (mb > Mbulge_rscale[NumMbulgeBins]) mb = Mbulge_rscale[NumMbulgeBins];
  
  if (md < Mdisc_rscale[1]) md = Mdisc_rscale[1];
  if (md > Mdisc_rscale[NumMdiscBins]) md = Mdisc_rscale[NumMdiscBins];
  
  if (lambda < Lambda_rscale[1]) lambda = Lambda_rscale[1];
  if (lambda > Lambda_rscale[NumLambdaBins]) 
    lambda = Lambda_rscale[NumLambdaBins];

  /* Find the indices that enclose the selected values */
  dlocate(Conc0_rscale, NumCbins, chalo, &ic[0]);
  ic[1] = ic[0] + 1;
#ifdef DEBUG
  printf("Conc0: %f < %f < %f (%d;%d)\n", Conc0_rscale[ic[0]], chalo, 
	 Conc0_rscale[ic[1]], ic[0], ic[1]);
#endif

  dlocate(Mbulge_rscale, NumMbulgeBins, mb, &imb[0]);
  imb[1] = imb[0] + 1;
#ifdef DEBUG
  printf("Mbulge: %f < %f < %f (%d;%d)\n", Mbulge_rscale[imb[0]], mb, 
	 Mbulge_rscale[imb[1]], imb[0], imb[1]);
#endif

  dlocate(Mdisc_rscale, NumMdiscBins, md, &imd);
#ifdef DEBUG
  printf("Mdisc: %f < %f < %f (%d;%d)\n", Mdisc_rscale[imd], md, 
	 Mdisc_rscale[imd+1], imd, imd+1);
#endif

  dlocate(Lambda_rscale, NumLambdaBins, lambda, &il);
#ifdef DEBUG
  printf("Lambda: %f < %f < %f (%d;%d)\n", Lambda_rscale[il], lambda, 
	 Lambda_rscale[il+1], il, il+1);
#endif
  
  /* We extract from the tables submatrices of size POLINT_POINTS x
     POLINT_POINTS, which are then provided to the 2D interpolation
     routine dpolin2 */
  kmd = IMIN(IMAX(imd-(POLINT_POINTS-1)/2,1), NumMdiscBins+1-POLINT_POINTS);
  kl = IMIN(IMAX(il-(POLINT_POINTS-1)/2,1), NumLambdaBins+1-POLINT_POINTS);
#ifdef DEBUG
  printf("kmd, kl, POLINT_POINTS: %d %d %d\n", kmd, kl, POLINT_POINTS);
  printf("Mdisc_rscale: ");
  for (p=0; p<POLINT_POINTS; p++)
    printf("%f ", Mdisc_rscale[kmd-1+p]);
  printf("\nLambda_rscale: ");
  for (p=0; p<POLINT_POINTS; p++)
    printf("%f ", Lambda_rscale[kl-1+p]);
  printf("\n\n");
#endif
  
  /* 
   * Eastern set of forking paths. 
   */
  /* For the first enclosing redshift, we get 4 values corresponding to 
     all combinations of bracketing concentration and bulge mass */
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      /*table = dsubmatrix(Table0_rscale[ic[i]][imb[j]], 
			 kmd-1,kmd-1+POLINT_POINTS,
			 kl-1,kl-1+POLINT_POINTS, 1,1);*/
      table = dsubmatrix(data0[ic[i]][imb[j]], 
			 kmd,kmd+POLINT_POINTS-1,
			 kl,kl+POLINT_POINTS-1, 1,1);
#ifdef DEBUG
      printf("table:\n");
      for (p=1; p<=POLINT_POINTS; p++) {
	for (q=1; q<=POLINT_POINTS; q++) {
	  printf("%8f ", table[p][q]);
	}
	printf("\n");
      }
      printf("\n"); fflush(stdout);
#endif
      dpolin2(&Mdisc_rscale[kmd-1], &Lambda_rscale[kl-1], table, 
	      POLINT_POINTS, POLINT_POINTS, md, lambda, &ymb0[i][j], &buf);
      free_dsubmatrix(table,1,POLINT_POINTS,1,POLINT_POINTS);
    }
  }
#ifdef DEBUG
  printf("ymb0: %f %f %f %f\n", ymb0[0][0], ymb0[0][1], ymb0[1][0], 
    ymb0[1][1]);
#endif

  /* First step backwards towards the garden's door: reduce mbulge */
  yc0[0] = ymb0[0][0] + (mb-Mbulge_rscale[imb[0]])*(ymb0[0][1]-ymb0[0][0])
    /(Mbulge_rscale[imb[1]]-Mbulge_rscale[imb[0]]);
  yc0[1] = ymb0[1][0] + (mb-Mbulge_rscale[imb[0]])*(ymb0[1][1]-ymb0[1][0])
    /(Mbulge_rscale[imb[1]]-Mbulge_rscale[imb[0]]);
#ifdef DEBUG
  printf("yc0: %f %f\n", yc0[0], yc0[1]);
#endif

  /* Second step towards the door: reduce c */
  yz0 = yc0[0] + (chalo-Conc0_rscale[ic[0]])*(yc0[1]-yc0[0])
    /(Conc0_rscale[ic[1]]-Conc0_rscale[ic[0]]);
#ifdef DEBUG
  printf("yz0 = %f\n", yz0);
#endif


  /* 
   * Western set of forking paths. 
   */
  dlocate(Conc1_rscale, NumCbins, chalo, &ic[0]);
  ic[1] = ic[0] + 1;
#ifdef DEBUG
  printf("Conc1: %f < %f < %f\n", Conc1_rscale[ic[0]], chalo, 
	 Conc1_rscale[ic[1]]);
#endif

  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      /*table = dsubmatrix(Table1_rscale[ic[i]][imb[j]], 
			 kmd-1,kmd-1+POLINT_POINTS,
			 kl-1,kl-1+POLINT_POINTS, 1,1);*/
      table = dsubmatrix(data1[ic[i]][imb[j]], 
			 kmd,kmd+POLINT_POINTS-1,
			 kl,kl+POLINT_POINTS-1, 1,1);
      dpolin2(&Mdisc_rscale[kmd-1], &Lambda_rscale[kl-1], table, 
	      POLINT_POINTS, POLINT_POINTS, md, lambda, &ymb1[i][j], &buf);
      free_dsubmatrix(table, 1,POLINT_POINTS, 1,POLINT_POINTS);
    }
  }
#ifdef DEBUG
  printf("ymb1: %f %f %f %f\n", ymb1[0][0], ymb1[0][1], ymb1[1][0], 
	 ymb1[1][1]);
#endif

  /* First step backwards reduce mbulge */
  yc1[0] = ymb1[0][0] + (mb-Mbulge_rscale[imb[0]])*(ymb1[0][1]-ymb1[0][0])
    /(Mbulge_rscale[imb[1]]-Mbulge_rscale[imb[0]]);
  yc1[1] = ymb1[1][0] + (mb-Mbulge_rscale[imb[0]])*(ymb1[1][1]-ymb1[1][0])
    /(Mbulge_rscale[imb[1]]-Mbulge_rscale[imb[0]]);
#ifdef DEBUG
  printf("yc1: %f %f\n", yc1[0], yc1[1]);
#endif

  /* Second step backwards: reduce c */
  yz1 = yc1[0] + (chalo-Conc1_rscale[ic[0]])*(yc1[1]-yc1[0])
    /(Conc1_rscale[ic[1]]-Conc1_rscale[ic[0]]);
#ifdef DEBUG
  printf("yz1 = %f\n", yz1);
#endif

  /* Final step to the garden's door. Interpolate between yz0 and yz1 */
  dlocate(Time_rscale, NumOutputTimes, zcurr, &iz[0]);
  iz[1] = iz[0] + 1;
  
  return yz0 + (zcurr-Time_rscale[iz[0]])*(yz1-yz0)
    /(Time_rscale[iz[1]]-Time_rscale[iz[0]]);
}


/**
 * @brief Free memory allocated to the scale radius tables
 */
void free_rscale_tables(void)
{
  free_dvector(Time_rscale, 1, NumOutputTimes);
  free_dvector(Mdisc_rscale, 1, NumMdiscBins);
  free_dvector(Lambda_rscale, 1, NumLambdaBins);
  free_dvector(Mbulge_rscale, 1, NumMbulgeBins);

  free_dvector(Conc0_rscale, 1, NumCbins);
  free_dvector(Conc1_rscale, 1, NumCbins);

  free_d4tensor(Table0_rscale, 1, NumCbins, 1, NumMbulgeBins, 
		1, NumMdiscBins, 1, NumLambdaBins);
  free_d4tensor(Table1_rscale, 1, NumCbins, 1, NumMbulgeBins, 
		1, NumMdiscBins, 1, NumLambdaBins);

  free_d4tensor(VcTable0_rscale, 1, NumCbins, 1, NumMbulgeBins, 
		1, NumMdiscBins, 1, NumLambdaBins);
  free_d4tensor(VcTable1_rscale, 1, NumCbins, 1, NumMbulgeBins, 
		1, NumMdiscBins, 1, NumLambdaBins);

  return;
}
