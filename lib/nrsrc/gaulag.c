#include <math.h>

#include "nrsag.h"

#define EPS 3.0e-14    /* Increase EPS if you don't have this precision */
#define MAXIT 10

/* 
 * Given alf, the parameter alpha of the Laguerre polynomials,
 * this routine returns arrays x[1..n] and w[1..n] containing the
 * abscissas and weights of the n-point Gauss-Laguerre quadrature
 * formula. The smallest abscissa is returned in x[1], the largest 
 * in x[n].
 */
void dgaulag(double x[], double w[], int n, double alf)
{
  double dgammln(double xx);
  void nrerror(char error_text[]);
  int i, its, j;
  double ai;
  double p1, p2, p3, pp, z, z1;  /* High precision is a good idea
				    for this routine. */
  
  for (i=1; i<=n; i++) {         /* Loop over the desired roots. */
    if (i == 1) {                /* Initial guess for the smallest root. */
      z = (1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
    } else if (i == 2) {         /* Initial guess for the second root. */
      z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
    } else {                     /* Initial guess for the other roots. */
      ai = i-2;
      z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
	    (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
    }
    for (its=1; its<=MAXIT; its++) { /* Refinement by Newton's method. */
      p1 = 1.0;
      p2 = 0.0;
      for (j=1; j<=n; j++) {     /* Loop up the recurrence relation to get
				    the Laguerre polynomial evaluated at z. */
	p3 = p2;
	p2 = p1;
	p1 = ((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
      }
      /* p1 is now the desired Laguerre polynomial. We next compute
	 pp, its derivative, by a standard relation involving also
	 p2, the polynomial of one lower order. */
      pp = (n*p1-(n+alf)*p2)/z;
      z1 = z;
      z = z1-p1/pp;              /* Newton's formula. */
      if (fabs(z-z1) <= EPS) break;
    }
    if (its > MAXIT) nrerror("too many iterations in gaulag");
    x[i] = z;                    /* Store the root and the weight. */
    w[i] = -exp(dgammln(alf+n)-dgammln((float)n))/(pp*n*p2);
  }
}
