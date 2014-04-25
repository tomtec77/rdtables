/**
 * @file concentration.c
 * @brief Read concentration vs mass tables and find the range of 
 * concentration values for each output time
 * @author Tomas E. Tecce
 */
#include "allvars.h"
#include "proto.h"
#include "readparameters.h"

/** Total number of table files */
#define NTABLES_CONC 51

/** Number of data lines in table files */
#define NLINES_CONC 101
 
/** Number of points for polynomial interpolation */
#define POLINT_POINTS 4

/** Maximum halo mass in Msun */
#define MHALOSAFETY 5e15 


/** 
 * @brief Load tables of concentration vs virial mass into memory for
 * quick access
 * @param *path Path to concentration tables
 *
 * Each file starts with a line with a single number (the redshift at 
 * which the table is calculated), then nlines lines of data
 */
void load_concentration_tables(char *path)
{
  FILE *fp;
  int i, j, check;
  char filename[FILENAME_MAX];
  double buf1, buf2, buf3;
#ifndef TESTMODE
  char *error;
  char buffer[256];
#endif

  Cmin_conc = dvector(1, NTABLES_CONC);
  Cmax_conc = dvector(1, NTABLES_CONC);
  Time_conc = dvector(1, NTABLES_CONC);

#ifdef TESTMODE
  printf("Reading %d tables c200 vs m200 from directory '%s'...\n",
	 NTABLES_CONC, path);

  Mass_conc = dvector(1, NLINES_CONC);
  Table_conc = dmatrix(1,NTABLES_CONC,1,NLINES_CONC);

  for (i=1; i<=NTABLES_CONC; i++) {

  /* Table files are numbered from 0 to NTABLES_CONC-1 */
    sprintf(filename, "%s/cmtable_%03d", path, i-1);
    if (!(fp=fopen(filename,"r"))) {
      fprintf(stderr, "Error (load_concentration_tables): "
	      "cannot open concentration table '%s': %s "
	      "(%u)\n", filename, strerror(errno), errno); fflush(stderr);
      exit(EXIT_FAILURE);
    }

    fscanf(fp, " %lf \n", &Time_conc[i]);      /* Read table redshift */
    if (i == 1) {
      printf("Starting redshift: %g\n", Time_conc[i]);
      if (Time_conc[i] < ZZ[NumOutputTimes-1]) {
	printf("#############################################################\n");
	printf("# WARNING: max redshift of tables < max redshift selected   #\n");
	printf("#############################################################\n");
	check = validate_continuation();
	if (check == 0) {
	  printf("Aborted on user input.\n");
	  exit(0);
	}
      }
    }
    if (i == NTABLES_CONC) 
      printf("Final redshift: %g\n", Time_conc[i]);

    /* Read data from file */
    for (j=1; j<=NLINES_CONC; j++) {
      fscanf(fp, " %lf %lf %lf %lf %lf \n", &Mass_conc[j], 
	     &Table_conc[i][j], &buf1, &buf2, &buf3); 
      /*#ifdef DEBUG
      if (j==1 || j==NLINES_CONC) {
	printf("i, j, Mass_conc, Table_conc: %d %d %g %g\n",
	       i, j, Mass_conc[j], Table_conc[i][j]);
	}
	#endif*/
    }
    fclose(fp);
  }

  /* Find range of concentration values for each table */
  find_table_limits();

  printf("Done reading c200 vs m200 tables.\n\n");
#else
  if (!(fp=fopen(path ,"r"))) {
    fprintf(stderr, "Error (load_concentration_tables): "
	    "cannot open concentration range table '%s': %s "
	    "(%u)\n", filename, strerror(errno), errno); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  /* Skip the first two lines (comments) */
  for (i=0; i<2; i++)
    error = fgets(buffer,sizeof(buffer),fp);

  for (i=1; i<=NTABLES_CONC; i++) {
    fscanf(fp, "%lf %lf %lf %lf %lf\n", &Time_conc[i], &Cmin_conc[i],
	   &Cmax_conc[i], &buf1, &buf2);

#ifdef DEBUG
    printf("i, Time_conc, Cmin_conc, Cmax_conc: %3d %6.3f %6.3f %6.3f\n",
	   i, Time_conc[i], Cmin_conc[i], Cmax_conc[i]);
#endif
  }
  printf("Done reading concentration range table.\n\n");
  
#endif

  return;
}


/**
 * @brief Find range of concentration values at a given redshift
 *
 * Determine, for each table, the range of concentration values and store
 * the maximum and minimum for each redshift in the vectors Cmin_conc and
 * Cmax_conc, respectively.
 */
void find_table_limits(void)
{
  int i, j;
  unsigned int offset;

  for (i=1; i<=NTABLES_CONC; i++) {
    Cmin_conc[i] = MAXC200;
    Cmax_conc[i] = MINC200;
  }
  
  printf("Concentration range:\n");
  printf("    Time         cmin        cmax\n");

  /* Search the tables for the maximum and minimum values */
  for (i=1; i<=NTABLES_CONC; i++) {
    offset = indexrm(i-1, 0, NTABLES_CONC, NLINES_CONC) + 1;

    for (j=1; j<=NLINES_CONC; j++) {
      /*printf("cmin, cmax, table: %g %g %g\n", Cmin_conc[i], Cmax_conc[i],
	Table_conc[i][j]);*/
      Cmin_conc[i] = dmin(Cmin_conc[i],Table_conc[i][j]);
      Cmax_conc[i] = dmax(Cmax_conc[i],Table_conc[i][j]);
    }
    printf("    %9.6f    %f    %f\n", Time_conc[i], Cmin_conc[i], 
	   Cmax_conc[i]);
  }

  return;
}


/**
 * @brief Find the range of valid concentrations at redshift z
 * @param z Redshift
 * @param *cmin Pointer to store minimum halo concentration
 * @param *cmax Pointer to store maximum halo concentration
 *
 * For a given redshift, interpolate in vectors Cmin_conc and Cmax_conc and
 * return the range of valid concentrations [cmin, cmax].
 */
void concentration_range(double z, double *cmin, double *cmax)
{
  unsigned long itime;
  int k;
  double cerr;


  /* Find where to interpolate */
  dlocate(Time_conc, NTABLES_CONC, z, &itime);
  if (itime == 0 || itime == NTABLES_CONC) {
    fprintf(stderr, "Error (concentration_range): redshift %g out of "
	    "range [%g - %g] - Exit\n", 
	    z, Time_conc[1], Time_conc[NTABLES_CONC]); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  k = IMIN(IMAX(itime-(POLINT_POINTS-1)/2,1),
	   NTABLES_CONC+1-POLINT_POINTS);
  dpolint(&Time_conc[k-1], &Cmin_conc[k-1], POLINT_POINTS, z, cmin, &cerr);
  dpolint(&Time_conc[k-1], &Cmax_conc[k-1], POLINT_POINTS, z, cmax, &cerr);
  
  return;
}


/** 
 * @brief Interpolate in the concentration tables to get the value of halo 
 * mass M_200 corresponding to a given concentration
 * @param chalo Halo concentration c_200 (assuming a NFW profile)
 * @param z Redshift
 * @param *result Pointer to store the result
 */
void get_mass_from_concentration(double chalo, double zcurr, double *result)
{
  int i, j;
  int k[2];
  unsigned long itime[2], ic[2];
  double *conc0, *conc1;
  double m200[2], merr[2];


  conc0 = dvector(1, NLINES_CONC);
  conc1 = dvector(1, NLINES_CONC);

  /* Find where to interpolate in time */
  dlocate(Time_conc, NTABLES_CONC, zcurr, &itime[1]);
  itime[0] = itime[1]+1;

  /* If z is outside the range of the tables, exit with error */
  if (itime[1] == 0 || itime[1] == NTABLES_CONC) {
    fprintf(stderr, "Error: selected redshift %g is outside the range "
	    "of the tables - Exit\n", zcurr); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  /* Get the corresponding concentration values */
  for (i=1; i<=NTABLES_CONC; i++) {
    for (j=1; j<=NLINES_CONC; j++) {
      conc0[j] = Table_conc[itime[0]][j];
      conc1[j] = Table_conc[itime[1]][j];
    }
  }

  /* Find where to interpolate in concentration */
  /* Note that here we need to locate in two arrays, not just one as when
     doing the inverse calculation (concentration from mass) since the mass
     ranges are the same in all tables, but that is not the case for the
     concentrations */
  
  /* Find where to interpolate in concentration */
  dlocate(conc0, NLINES_CONC, chalo, &ic[0]);
  dlocate(conc1, NLINES_CONC, chalo, &ic[1]);

  for (i=0; i<2; i++)
    k[i] = IMIN(IMAX(ic[i]-(POLINT_POINTS-1)/2,1),
		NLINES_CONC+1-POLINT_POINTS);

  dpolint(&conc0[k[0]-1], &Mass_conc[k[0]-1], POLINT_POINTS, chalo,
	  &m200[0], &merr[0]);
  dpolint(&conc1[k[1]-1], &Mass_conc[k[1]-1], POLINT_POINTS, chalo,
	  &m200[1], &merr[1]);

  free_dvector(conc0, 1, NLINES_CONC);
  free_dvector(conc1, 1, NLINES_CONC);

  /* Return a linear interpolation between the results for the two
     redshifts */
  *result =  m200[0] + (zcurr-Time_conc[itime[0]])*(m200[1]-m200[0])/
    (Time_conc[itime[1]]-Time_conc[itime[0]]);

  return;
}


/**
 * @brief Free memory allocated to concentration tables
 */
void free_concentration_tables(void)
{
  free_dvector(Cmin_conc,1,NTABLES_CONC);
  free_dvector(Cmax_conc,1,NTABLES_CONC);
  free_dvector(Time_conc,1,NTABLES_CONC);

#ifdef TESTMODE
  free_dvector(Mass_conc,1,NLINES_CONC);
  free_dmatrix(Table_conc,1,NTABLES_CONC,1,NLINES_CONC);
#endif

  return;
}


