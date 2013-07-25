
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Mateusz Markowski                                */
/* and modified by Przemyslaw Kiciak                                         */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "g1blendingf.h"

static const float coeff[] =
  { 13.0,   53.0,   48.0,   53.0, 13.0,
    53.0, -272.0, -282.0, -272.0, 53.0,
    48.0, -282.0, 1548.0, -282.0, 48.0};

boolean g1bl_SetupClosedBiharmAMatrixf ( int lastknotu, int lastknotv,
                               int *n, int **prof, float **Amat, float ***arow )
{
  int   bldim, blnum, _n, asize;
  int   *_prof;
  float *_Amat, **_arow;
  int   i, j, k, j0, j1;

        /* compute the size and number of blocks */
  _prof = NULL;
  _Amat = NULL;
  _arow = NULL;
  bldim = lastknotv-6;
  blnum = lastknotu-4;
  if ( bldim < 1 || blnum < 1 )
    goto failure;
  _n = bldim*blnum;

  PKV_MALLOC ( _prof, _n*sizeof(int) )
  if ( !_prof )
    goto failure;

  for ( i = 0; i < blnum; i++ )
    for ( j = 0; j < bldim; j++ ) {
      _prof[i*bldim+j] = max ( 0, j-2 );
      if ( i > 2 && i < blnum-2 )
        _prof[i*bldim+j] += (i-2)*bldim;
    }
        /* compute the size of the coefficients array */
  asize = pkn_NRBArraySize ( _n, _prof );
        /* setup the blocks */
  PKV_MALLOC ( _Amat, asize*sizeof(float) )
  PKV_MALLOC ( _arow, _n*sizeof(float*) )
  if ( !_Amat || !_arow )
    goto failure;
  memset ( _Amat, 0, asize*sizeof(float) );
  pkn_NRBFindRowsf ( _n, _prof, _Amat, _arow );
  for ( i = k = 0;  i < blnum;  i++ ) {
    for ( j = 0;  j < bldim;  j++, k++ ) {
      j0 = max ( 0, j-2 );
      j1 = min ( 3, bldim-j );
      memcpy ( &_arow[k][k-(j-j0)], &coeff[10+2-(j-j0)],
               (j-j0+1)*sizeof(float) );
      if ( i > 0 )
        memcpy ( &_arow[k][k-(j-j0)-bldim], &coeff[5+2-(j-j0)],
                 (j-j0+j1)*sizeof(float) );
      if ( i > 1 )
        memcpy ( &_arow[k][k-(j-j0)-2*bldim], &coeff[2-(j-j0)],
                 (j-j0+j1)*sizeof(float) );
    }
  }
  for ( i = blnum-2, k = i*bldim;  i < blnum;  i++ ) {
    for ( j = 0;  j < bldim;  j++, k++ ) {
      j0 = max ( 0, j-2 );
      j1 = min ( 3, bldim-j );
      memcpy ( &_arow[k][(i-blnum+2)*bldim+j0], &coeff[2-(j-j0)],
               (j-j0+j1)*sizeof(float) );
      if ( i > blnum-2 )
        memcpy ( &_arow[k][(i-blnum+1)*bldim+j0], &coeff[5+2-(j-j0)],
                 (j-j0+j1)*sizeof(float) );
    }
  }
  *n = _n;
  *prof = _prof;
  *Amat = _Amat;
  *arow = _arow;
  return true;

failure:
  if ( !_prof ) PKV_FREE ( _prof )
  if ( !_Amat ) PKV_FREE ( _Amat )
  if ( !_arow ) PKV_FREE ( _arow )
  *prof = NULL;
  *Amat = NULL;
  *arow = NULL;
  return false;
} /*g1bl_SetupClosedBiharmAMatrixf*/

boolean g1bl_SetupClosedBiharmRHSf ( int lastknotu, int lastknotv,
                               int spdimen, int pitch, const float *cpoints,
                               float *rhs )
{
  int bldim, blnum, n;
  int i, j, k, l;

  bldim = lastknotv-6;
  blnum = lastknotu-4;
  if ( bldim < 1 || blnum < 1 )
    return false;
  n = bldim*blnum;

  memset ( rhs, 0, spdimen*n*sizeof(float) );
        /* for each column of the control net */
  for ( i = 0; i < blnum; i++ ) {
    for ( j = 0; j < 5; j++ ) {
      k = 2 - abs( j-2 );
      l = i - (j-2);
      if ( l < 0 )           l += blnum;
      else if ( l >= blnum ) l -= blnum;
      pkn_MultMatrixSubf ( 1, 2, 0, &coeff[5*k], spdimen, spdimen,
                 &cpoints[i*pitch], 0, &rhs[l*bldim*spdimen] );
      if ( bldim > 1 )
        pkn_MultMatrixSubf ( 1, 1, 0, &coeff[5*k], spdimen, spdimen,
                 &cpoints[i*pitch+spdimen], 0, &rhs[(l*bldim+1)*spdimen] );
      pkn_MultMatrixSubf ( 1, 2, 0, &coeff[5*k+3], spdimen, spdimen,
                 &cpoints[i*pitch+(lastknotv-4)*spdimen], 0,
                 &rhs[((l+1)*bldim-1)*spdimen]);
      if ( bldim > 1 )
        pkn_MultMatrixSubf ( 1, 1, 0, &coeff[5*k+4], spdimen, spdimen,
                 &cpoints[i*pitch+(lastknotv-4)*spdimen], 0,
                 &rhs[((l+1)*bldim-2)*spdimen] );
    }
  }
  return true;
} /*g1bl_SetupClosedBiharmRHSf*/

