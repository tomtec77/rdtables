#include <stdio.h>
#include "nrsag.h"

/*
  Given an array xx[1..n] and given a value x, returns a value j 
  such that x is between xx[j] and xx[j+1]. xx must be monotonic, 
  either increasing or decreasing.
  j = 0 or j = n is returned to indicate that x is out of range.
*/
void dlocate(double xx[], unsigned long n, double x, unsigned long *j)
{
  unsigned long ju, jm, jl;
  int ascnd;
  /*int i;*/

  jl = 0;                       // Initialise lower
  ju = n+1;                     // and upper limits.
  ascnd = (xx[n] >= xx[1]);
  while (ju-jl > 1) {           // If we are not yet done,
    jm = (ju+jl) >> 1;          // compute a midpoint,
    if (x >= xx[jm] == ascnd)
      jl = jm;                  // and replace either the lower limit
    else
      ju = jm;                  // or the upper limit, as appropriate.
  }                             // Repeat until the test condition is satisfied.
  if (x == xx[1]) *j = 1;       // Then set the output
  else if (x == xx[n]) *j = n-1;
  else *j = jl;
}                               // and return.
