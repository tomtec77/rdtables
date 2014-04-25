/*
  indexrm.c

  Returns the row-major offset of element (i,j) of a matrix of dimensions
  sizex X sizey. Element indices are assumed to run from 0 to size-1. i is
  the row index, j the column index.
*/
#include <stdio.h>
#include <stdlib.h>

#include "myutils.h"

unsigned int indexrm(int i, int j, int sizex, int sizey)
{
  if (i>sizex || j>sizey) {
    fprintf(stderr, "Error (indexrm): index value too large - Exit\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  return (unsigned int)(j + (sizey*i));
} 
