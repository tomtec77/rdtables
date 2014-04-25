/**
 * @file nrsag.h
 * @brief Custom version of nr.h for use with SAG
 *
 * In this adaptation of the Numerical Recipes code, many programs
 * have been converted to double precision.
 */
double dbessi0(double);
double dbessi1(double);
double dbessk0(double);
double dbessk1(double);

double dqromb(double (*func)(double), double, double);
void dpolint(double [], double [], int, double, double *, double *);
double dtrapzd(double (*func)(double), double, double, int);

void dgaulag(double [], double [], int, double);
double dgammln(double);
void dlocate(double [], unsigned long, double, unsigned long *);
void dpolin2(double [], double [], double **, int, int, double, double, double *, double *);
