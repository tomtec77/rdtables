/**
 * @file polin2.c
 * @brief Function for polinomial interpolation
 *
 * Given arrays x1a[1..m] and x2a[1..n] of independent variables, and a 
 * submatrix of function values ya[1..m][1..n], tabulated at the grid points
 * defined by x1a and x2a; and given values x1 and x2 of the independent
 * variables; this routine returns an interpolated function value y, and
 * an accuracy indication dy (based only on the interpolation in the x1
 * direction, however).
 */
#include "nr.h"
#include "nrutil.h"
#include "nrsag.h"

void dpolin2(double x1a[], double x2a[], double **ya, int m, int n, 
	     double x1, double x2, double *y, double *dy)
{
  int j;
  double *ymtmp;
  
  /*printf("dpolin2 ");*/
  ymtmp = dvector(1,m);
  for (j=1; j<=m; j++) {                              /* Loop over rows. */
    dpolint(x2a, ya[j], n, x2, &ymtmp[j], dy);  /* Interpolate answer into
						   temporary storage. */
    /*printf("%f ", ymtmp[j]);*/
  }
  /*printf("\n");*/
  dpolint(x1a, ymtmp, m, x1, y, dy);      /* Do the final interpolation. */
  free_dvector(ymtmp,1,m);
  
  return;
}
