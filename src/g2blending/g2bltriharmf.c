
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "g2blendingf.h"

static const float coeff[] =
  { -22.0, -309.0,   -666.0,   -526.0,   -666.0, -309.0,  -22.0,
   -309.0,  720.0,   4941.0,   4416.0,   4941.0,  720.0, -309.0,
   -666.0, 4941.0, -15030.0, -16290.0, -15030.0, 4941.0, -666.0,
   -526.0, 4416.0, -16290.0,  75200.0, -16290.0, 4416.0, -526.0};

boolean g2bl_SetupTriharmAMatrixf ( int lastknotu, int lastknotv,
                               int *n, int **prof, float **Amat, float ***arow )
{
  int    bldim, blnum, _n, asize;
  int    *_prof;
  float *_Amat, **_arow;
  int    i, j, k, j0, j1;

        /* compute the size and number of blocks */
  _prof = NULL;
  _Amat = NULL;
  _arow = NULL;
  bldim = lastknotv-9;
  blnum = lastknotu-9;
  if ( bldim < 1 || blnum < 1 )
    goto failure;
  _n = bldim*blnum;

  PKV_MALLOC ( _prof, _n*sizeof(int) )
  if ( !_prof )
    goto failure;

  for ( i = 0; i < blnum; i++ )
    for ( j = 0; j < bldim; j++ ) {
      _prof[i*bldim+j] = max ( 0, j-3 );
      if ( i > 3 )
        _prof[i*bldim+j] += (i-3)*bldim;
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
      j0 = max ( 0, j-3 );
      j1 = min ( 4, bldim-j );
      memcpy ( &_arow[k][k-(j-j0)], &coeff[21+3-(j-j0)],
               (j-j0+1)*sizeof(float) );
      if ( i > 0 )
        memcpy ( &_arow[k][k-(j-j0)-bldim], &coeff[14+3-(j-j0)],
                 (j-j0+j1)*sizeof(float) );
      if ( i > 1 )
        memcpy ( &_arow[k][k-(j-j0)-2*bldim], &coeff[7+3-(j-j0)],
                 (j-j0+j1)*sizeof(float) );
      if ( i > 2 )
        memcpy ( &_arow[k][k-(j-j0)-3*bldim], &coeff[3-(j-j0)],
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
} /*g2bl_SetupTriharmAMatrixf*/

static void RHSInnerBlock ( int spdimen, int pitch, const float *cpoints,
                            int c1, int bldim, int i, int j, float *rhs )
{
  pkn_MultMatrixSubf ( 1, 3, 0, &coeff[c1], spdimen, spdimen,
               &cpoints[(i+3)*pitch], 0, &rhs[j*bldim*spdimen] );
  if ( bldim > 1 )
    pkn_MultMatrixSubf ( 1, 2, 0, &coeff[c1], spdimen, spdimen,
                 &cpoints[(i+3)*pitch+spdimen], 0, &rhs[(j*bldim+1)*spdimen] );
  if ( bldim > 2 ) {
    pkn_MultMatrixSubf ( 1, 1, 0, &coeff[c1], spdimen, spdimen,
                 &cpoints[(i+3)*pitch+2*spdimen], 0, &rhs[(j*bldim+2)*spdimen] );
    pkn_MultMatrixSubf ( 1, 1, 0, &coeff[c1+6], spdimen, spdimen,
                 &cpoints[(i+3)*pitch+(bldim+3)*spdimen], 0,
                 &rhs[((j+1)*bldim-3)*spdimen] );
  }
  if ( bldim > 1 )
    pkn_MultMatrixSubf ( 1, 2, 0, &coeff[c1+5], spdimen, spdimen,
                 &cpoints[(i+3)*pitch+(bldim+3)*spdimen], 0,
                 &rhs[((j+1)*bldim-2)*spdimen] );
  pkn_MultMatrixSubf ( 1, 3, 0, &coeff[c1+4], spdimen, spdimen,
               &cpoints[(i+3)*pitch+(bldim+3)*spdimen], 0,
               &rhs[((j+1)*bldim-1)*spdimen] );
} /*RHSInnerBlock*/

boolean g2bl_SetupTriharmRHSf ( int lastknotu, int lastknotv,
                           int spdimen, int pitch, const float *cpoints,
                           float *rhs )
{
  int bldim, blnum, n;
  int i, j;

  bldim = lastknotv-9;
  blnum = lastknotu-9;
  if ( bldim < 1 || blnum < 1 )
    return false;
  n = bldim*blnum;

  memset ( rhs, 0, spdimen*n*sizeof(float) );
        /* the first three columns of the control net */
  for ( j = 0; j < bldim; j++ ) {
    pkn_MultMatrixSubf ( 1, 7, 0, &coeff[0], spdimen, spdimen,
                 &cpoints[j*spdimen], 0, &rhs[j*spdimen] );
    pkn_MultMatrixSubf ( 1, 7, 0, &coeff[7], spdimen, spdimen,
                 &cpoints[j*spdimen+pitch], 0, &rhs[j*spdimen] );
    pkn_MultMatrixSubf ( 1, 7, 0, &coeff[14], spdimen, spdimen,
                 &cpoints[j*spdimen+2*pitch], 0, &rhs[j*spdimen] );
    if ( j+bldim < n ) {
      pkn_MultMatrixSubf ( 1, 7, 0, &coeff[0], spdimen, spdimen,
                   &cpoints[j*spdimen+pitch], 0, &rhs[(j+bldim)*spdimen] );
      pkn_MultMatrixSubf ( 1, 7, 0, &coeff[7], spdimen, spdimen,
                   &cpoints[j*spdimen+2*pitch], 0, &rhs[(j+bldim)*spdimen] );
    }
    if ( j+2*bldim < n )
      pkn_MultMatrixSubf ( 1, 7, 0, &coeff[0], spdimen, spdimen,
                   &cpoints[j*spdimen+2*pitch], 0, &rhs[(j+2*bldim)*spdimen] );
  }
        /* the first three and last three points of the middle columns */
  for ( i = 0; i < lastknotu-9; i++ ) {
    RHSInnerBlock ( spdimen, pitch, cpoints, 21, bldim, i, i, rhs );
    if ( i > 0 )
      RHSInnerBlock ( spdimen, pitch, cpoints, 14, bldim, i-1, i, rhs );
    if ( i > 1 )
      RHSInnerBlock ( spdimen, pitch, cpoints, 7, bldim, i-2, i, rhs );
    if ( i > 2 )
      RHSInnerBlock ( spdimen, pitch, cpoints, 0, bldim, i-3, i, rhs );
    if ( i < lastknotu-10 )
      RHSInnerBlock ( spdimen, pitch, cpoints, 14, bldim, i+1, i, rhs );
    if ( i < lastknotu-11 )
      RHSInnerBlock ( spdimen, pitch, cpoints, 7, bldim, i+2, i, rhs );
    if ( i < lastknotu-12 )
      RHSInnerBlock ( spdimen, pitch, cpoints, 0, bldim, i+3, i, rhs );
  }
        /* the last three columns of the control net */
  for ( j = 0; j < bldim; j++ ) {
    pkn_MultMatrixSubf ( 1, 7, 0, &coeff[14], spdimen, spdimen,
                         &cpoints[(blnum+3)*pitch+j*spdimen],
                         0, &rhs[((blnum-1)*bldim+j)*spdimen] );
    pkn_MultMatrixSubf ( 1, 7, 0, &coeff[7], spdimen, spdimen,
                         &cpoints[(blnum+4)*pitch+j*spdimen],
                         0, &rhs[((blnum-1)*bldim+j)*spdimen] );
    pkn_MultMatrixSubf ( 1, 7, 0, &coeff[0], spdimen, spdimen,
                         &cpoints[(blnum+5)*pitch+j*spdimen],
                         0, &rhs[((blnum-1)*bldim+j)*spdimen] );
    if ( blnum >= 2 ) {
      pkn_MultMatrixSubf ( 1, 7, 0, &coeff[7], spdimen, spdimen,
                           &cpoints[(blnum+3)*pitch+j*spdimen],
                           0, &rhs[((blnum-2)*bldim+j)*spdimen] );
      pkn_MultMatrixSubf ( 1, 7, 0, &coeff[0], spdimen, spdimen,
                           &cpoints[(blnum+4)*pitch+j*spdimen],
                           0, &rhs[((blnum-2)*bldim+j)*spdimen] );
    }
    if ( blnum >= 3 )
      pkn_MultMatrixSubf ( 1, 7, 0, &coeff[0], spdimen, spdimen,
                           &cpoints[(blnum+3)*pitch+j*spdimen],
                           0, &rhs[((blnum-3)*bldim+j)*spdimen] );
  }

  return true;
} /*g2bl_SetupTriharmRHSf*/

