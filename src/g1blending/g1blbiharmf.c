
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Mateusz Markowski                                */

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

boolean g1bl_SetupBiharmAMatrixf ( int lastknotu, int lastknotv,
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
  blnum = lastknotu-6;
  if ( bldim < 1 || blnum < 1 )
    goto failure;
  _n = bldim*blnum;

  PKV_MALLOC ( _prof, _n*sizeof(int) )
  if ( !_prof )
    goto failure;

  for ( i = 0; i < blnum; i++ )
    for ( j = 0; j < bldim; j++ ) {
      _prof[i*bldim+j] = max ( 0, j-2 );
      if ( i > 2 )
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
} /*g1bl_SetupBiharmAMatrixf*/

static void RHSInnerBlock ( int spdimen, int pitch, const float *cpoints,
                            int c1, int bldim, int i, int j,
                            float *rhs )
{
  pkn_MultMatrixSubf ( 1, 2, 0, &coeff[c1], spdimen, spdimen,
               &cpoints[(i+2)*pitch], 0, &rhs[j*bldim*spdimen] );
  if ( bldim > 1 ) {
    pkn_MultMatrixSubf ( 1, 1, 0, &coeff[c1], spdimen, spdimen,
                 &cpoints[(i+2)*pitch+spdimen], 0, &rhs[(j*bldim+1)*spdimen] );
    pkn_MultMatrixSubf ( 1, 1, 0, &coeff[c1+4], spdimen, spdimen,
                 &cpoints[(i+2)*pitch+(bldim+2)*spdimen], 0,
                 &rhs[((j+1)*bldim-2)*spdimen] );
  }
  pkn_MultMatrixSubf ( 1, 2, 0, &coeff[c1+3], spdimen, spdimen,
               &cpoints[(i+2)*pitch+(bldim+2)*spdimen], 0,
               &rhs[((j+1)*bldim-1)*spdimen] );
} /*RHSInnerBlock*/

boolean g1bl_SetupBiharmRHSf ( int lastknotu, int lastknotv,
                               int spdimen, int pitch, const float *cpoints,
                               float *rhs )
{
  int bldim, blnum, n;
  int i, j;

  bldim = lastknotv-6;
  blnum = lastknotu-6;
  if ( bldim < 1 || blnum < 1 )
    return false;
  n = bldim*blnum;

  memset ( rhs, 0, spdimen*n*sizeof(float) );
        /* the first two columns of the control net */
  for ( j = 0; j < bldim; j++ ) {
    pkn_MultMatrixSubf ( 1, 5, 0, &coeff[0], spdimen, spdimen,
                 &cpoints[j*spdimen], 0, &rhs[j*spdimen] );
    pkn_MultMatrixSubf ( 1, 5, 0, &coeff[5], spdimen, spdimen,
                 &cpoints[j*spdimen+pitch], 0, &rhs[j*spdimen] );
    if ( j+bldim < n )
      pkn_MultMatrixSubf ( 1, 5, 0, &coeff[0], spdimen, spdimen,
                   &cpoints[j*spdimen+pitch], 0, &rhs[(j+bldim)*spdimen] );
  }
        /* the first two and last two points of the middle columns */
  for ( i = 0; i < blnum; i++ ) {
    RHSInnerBlock ( spdimen, pitch, cpoints, 10, bldim, i, i, rhs );
    if ( i > 0 )
      RHSInnerBlock ( spdimen, pitch, cpoints, 5, bldim, i-1, i, rhs );
    if ( i > 1 )
      RHSInnerBlock ( spdimen, pitch, cpoints, 0, bldim, i-2, i, rhs );
    if ( i < lastknotu-7 )
      RHSInnerBlock ( spdimen, pitch, cpoints, 5, bldim, i+1, i, rhs );
    if ( i < lastknotu-8 )
      RHSInnerBlock ( spdimen, pitch, cpoints, 0, bldim, i+2, i, rhs );
  }
        /* the last two columns of the control net */
  for ( j = 0; j < bldim; j++ ) {
    pkn_MultMatrixSubf ( 1, 5, 0, &coeff[5], spdimen, spdimen,
                         &cpoints[(blnum+2)*pitch+j*spdimen],
                         0, &rhs[((blnum-1)*bldim+j)*spdimen] );
    pkn_MultMatrixSubf ( 1, 5, 0, &coeff[0], spdimen, spdimen,
                         &cpoints[(blnum+3)*pitch+j*spdimen],
                         0, &rhs[((blnum-1)*bldim+j)*spdimen] );
    if ( blnum >= 2 )
      pkn_MultMatrixSubf ( 1, 5, 0, &coeff[0], spdimen, spdimen,
                           &cpoints[(blnum+2)*pitch+j*spdimen],
                           0, &rhs[((blnum-2)*bldim+j)*spdimen] );
  }

  return true;
} /*g1bl_SetupBiharmRHSf*/

