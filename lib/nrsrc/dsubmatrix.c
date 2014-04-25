/*
  dsubmatrix.c

  Double precision version of submatrix.c
*/
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "nrutil.h"

#define NR_END 1
#define FREE_ARG char*

double **dsubmatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* Point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	double **m;

	/* Allocate array of pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in dsubmatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_dsubmatrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by dsubmatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}
