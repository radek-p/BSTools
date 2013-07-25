
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

/* /////////////////////////////////////////////////////////////////////////// */
boolean pkn_SPsubMSortByRows ( unsigned int nrows, unsigned int ncols,
                               unsigned int nnz, index3 *ai, unsigned int *permut )
{
  int  i;

  if ( nnz < 1 )
    return false;
  else if ( nnz > 1 ) {
    for ( i = 0; i < nnz; i++ )
      permut[i] = i;
    if ( pkv_SortKernel ( sizeof(int), ID_UNSIGNED, sizeof(index3),
                          0, nnz,  &ai->j, permut ) != SORT_OK )
      return false;
    if ( pkv_SortKernel ( sizeof(int), ID_UNSIGNED, sizeof(index3),
                          0, nnz, &ai->i, permut ) != SORT_OK )
      return false;
  }
  else
    permut[0] = 0;
  return true;
} /*pkn_SPsubMSortByRows*/

boolean pkn_SPsubMSortByCols ( unsigned int nrows, unsigned int ncols,
                               unsigned int nnz, index3 *ai, unsigned int *permut )
{
  int  i;

  if ( nnz < 1 )
    return false;
  else if ( nnz > 1 ) {
    for ( i = 0; i < nnz; i++ )
      permut[i] = i;
    if ( pkv_SortKernel ( sizeof(int), ID_UNSIGNED, sizeof(index3),
                          0, nnz,  &ai->i, permut ) != SORT_OK )
      return false;
    if ( pkv_SortKernel ( sizeof(int), ID_UNSIGNED, sizeof(index3),
                          0, nnz, &ai->j, permut ) != SORT_OK )
      return false;
  }
  else
    permut[0] = 0;
  return true;
} /*pkn_SPsubMSortByCols*/

boolean pkn_SPsubMFindRows ( unsigned int nrows, unsigned int ncols, unsigned int nnz,
                             index3 *ai, unsigned int *permut, boolean ro,
                             int *rows )
{
  int i, k, l, c;

  if ( !ro ) {
    if ( !pkn_SPsubMSortByRows ( nrows, ncols, nnz, ai, permut ) )
      return false;
  }
  memset ( rows, 0, (nrows+1)*sizeof(int) );
  for ( k = c = 0, i = -1;  k < nnz;  k++ ) {
    l = permut[k];
    if ( ai[l].i != i ) {
      rows[i+1] = c;
      c = 1;
      i = ai[l].i;
    }
    else
      c ++;
  }
  rows[i+1] = c;
  for ( i = 1; i <= nrows; i++ )
    rows[i] += rows[i-1];
  return true;
} /*pkn_SPsubMFindRows*/

boolean pkn_SPsubMFindCols ( unsigned int nrows, unsigned int ncols, unsigned int nnz,
                             index3 *ai, unsigned int *permut, boolean co,
                             int *cols )
{
  int j, k, l, c;

  if ( !co ) {
    if ( !pkn_SPsubMSortByCols ( nrows, ncols, nnz, ai, permut ) )
      return false;
  }
  memset ( cols, 0, (ncols+1)*sizeof(int) );
  for ( k = c = 0, j = -1;  k < nnz;  k++ ) {
    l = permut[k];
    if ( ai[l].j != j ) {
      cols[j+1] = c;
      c = 1;
      j = ai[l].j;
    }
    else
      c ++;
  }
  cols[j+1] = c;
  for ( j = 1; j <= ncols; j++ )
    cols[j] += cols[j-1];
  return true;
} /*pkn_SPsubMFindCols*/

