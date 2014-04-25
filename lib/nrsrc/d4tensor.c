/**
 * @file d4tensor.c
 * @brief Allocate a double 4tensor with range 
 * t[nrl..nrh][ncl..nch][ndl..ndh][nel..neh] 
 */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "nrutil.h"

#define NR_END 1
#define FREE_ARG char*

double ****d4tensor(long nrl, long nrh, 
		    long ncl, long nch, 
		    long ndl, long ndh,
		    long nel, long neh)
{
  long a,b,c, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1,net=neh-nel+1;

  double ****t;

  /* Allocate pointers to pointers to rows */
  t=(double ****) malloc((size_t)((nrow+NR_END)*sizeof(double***)));
  if (!t) nrerror("allocation failure 1 in d4tensor()");
  t += NR_END;
  t -= nrl;

  /* Allocate pointers to rows and set pointers to them */
  t[nrl]=(double ***) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double**)));
  if (!t[nrl]) nrerror("allocation failure 2 in d4tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* Allocate rows and set pointers to them */
  t[nrl][ncl]=(double **) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double*)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in d4tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
 
  /* Allocate ... */
  t[nrl][ncl][ndl]=(double *) malloc((size_t)((nrow*ncol*ndep*net+NR_END)*sizeof(double)));
  if (!t[nrl][ncl][ndl]) nrerror("allocation failure 4 in d4tensor()");
  t[nrl][ncl][ndl] += NR_END;
  t[nrl][ncl][ndl] -= nel;

  /* I */
  for(a=ndl+1;a<=ndh;a++)
    t[nrl][ncl][a]=t[nrl][ncl][a-1]+net;

  for(a=ncl+1;a<=nch;a++)
    {
      /* II */
      t[nrl][a]=t[nrl][a-1]+ndep;
      
      /* III */
      t[nrl][a][ndl]=t[nrl][a-1][ndl]+ndep*net;

      /* IV */
      for(b=ndl+1;b<=ndh;b++)
	t[nrl][a][b]=t[nrl][a][b-1]+net; 
    }

  for(a=nrl+1;a<=nrh;a++)
    {
      /* V */
      t[a]=t[a-1]+ncol;
      
      /* VI */
      t[a][ncl]=t[a-1][ncl]+ncol*ndep; //ndep**net;
      
      /* VII */
      t[a][ncl][ndl]=t[a-1][ncl][ndl]+ncol*ndep*net;

      /* VIII */
      for(b=ndl+1;b<=ndh;b++)
	t[a][ncl][b]=t[a][ncl][b-1]+net;

      for(b=ncl+1;b<=nch;b++)
	{
	  /* IX */
	  t[a][b]=t[a][b-1]+ndep;
	  
	  /* X */
	  t[a][b][ndl]=t[a][b-1][ndl]+ndep*net;

	  /* XI */
	  for(c=ndl+1;c<=ndh;c++)
	    t[a][b][c]=t[a][b][c-1]+net;
	}
    }
  return (t);
}


/* 
 * Free a double d4tensor allocated by d4tensor() 
 */
void free_d4tensor(double ****t, 
		   long nrl, long nrh, 
		   long ncl, long nch, 
		   long ndl, long ndh, 
		   long nel, long neh)
{
  free((FREE_ARG) (t[nrl][ncl][ndl]+nel-NR_END));
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}
