/*
  Header file for myutils library (definitions and prototypes)

  T. E. Tecce 2010
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/utsname.h>


unsigned int SeedX, SeedY, SeedZ, SeedC;

double dmin(double, double);
double dmax(double, double);

unsigned int KISS();
double rndnum();
double gaussdev();
void rnd_setup(unsigned int);
unsigned int devrand();

void errormsg(char *, int);

void checkwrite(int, int);

unsigned int indexrm(int, int, int, int);

unsigned int validate_continuation(void);

int mygethostname(char *, int);
