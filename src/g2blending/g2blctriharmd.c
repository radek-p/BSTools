
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

#include "g2blendingd.h"
#include "msgpool.h"

static const double coeff[] =
  { -22.0, -309.0,   -666.0,   -526.0,   -666.0, -309.0,  -22.0,
   -309.0,  720.0,   4941.0,   4416.0,   4941.0,  720.0, -309.0,
   -666.0, 4941.0, -15030.0, -16290.0, -15030.0, 4941.0, -666.0,
   -526.0, 4416.0, -16290.0,  75200.0, -16290.0, 4416.0, -526.0};

boolean g2bl_SetupClosedTriharmAMatrixd ( int lastknotu, int lastknotv,
                           int *n, int **prof, double **Amat, double ***arow )
{
  int    bldim, blnum, _n, asize;
  int    *_prof;
  double *_Amat, **_arow;
  int    i, j, k, j0, j1;

        /* compute the size and number of blocks */
  _prof = NULL;
  _Amat = NULL;
  _arow = NULL;
  bldim = lastknotv-9;
  blnum = lastknotu-6;
  if ( bldim < 1 || blnum < 7 )
    goto failure;
  _n = bldim*blnum;

  PKV_MALLOC ( _prof, _n*sizeof(int) );
  if ( !_prof ) {
    pkv_SignalError ( LIB_G2BLENDING, 1, ERRMSG_1 );
    goto failure;
  }

  for ( i = 0; i < blnum; i++ )
    for ( j = 0; j < bldim; j++ ) {
      _prof[i*bldim+j] = max ( 0, j-3 );
      if ( i > 3 && i < blnum-3 )
        _prof[i*bldim+j] += (i-3)*bldim;
    }
        /* compute the size of the coefficients array */
  asize = pkn_NRBArraySize ( _n, _prof );
        /* setup the blocks */
  PKV_MALLOC ( _Amat, asize*sizeof(double) )
  PKV_MALLOC ( _arow, _n*sizeof(double*) );
  if ( !_Amat || !_arow ) {
    pkv_SignalError ( LIB_G2BLENDING, 2, ERRMSG_1 );
    goto failure;
  }
  memset ( _Amat, 0, asize*sizeof(double) );
  pkn_NRBFindRowsd ( _n, _prof, _Amat, _arow );
  for ( i = k = 0;  i < blnum;  i++ ) {   
    for ( j = 0;  j < bldim;  j++, k++ ) {
      j0 = max ( 0, j-3 );
      j1 = min ( 4, bldim-j );
      memcpy ( &_arow[k][k-(j-j0)], &coeff[21+3-(j-j0)],
               (j-j0+1)*sizeof(double) );
      if ( i > 0 )
        memcpy ( &_arow[k][k-(j-j0)-bldim], &coeff[14+3-(j-j0)],
                 (j-j0+j1)*sizeof(double) );
      if ( i > 1 )
        memcpy ( &_arow[k][k-(j-j0)-2*bldim], &coeff[7+3-(j-j0)],
                 (j-j0+j1)*sizeof(double) );
      if ( i > 2 )
        memcpy ( &_arow[k][k-(j-j0)-3*bldim], &coeff[3-(j-j0)],
                 (j-j0+j1)*sizeof(double) );
    }
  }
  for ( i = blnum-3, k = i*bldim;  i < blnum;  i++ ) {
    for ( j = 0;  j < bldim;  j++, k++ ) {
      j0 = max ( 0, j-3 );
      j1 = min ( 4, bldim-j );
      memcpy ( &_arow[k][(i-blnum+3)*bldim+j0], &coeff[3-(j-j0)],
               (j-j0+j1)*sizeof(double) );
      if ( i > blnum-3 )
        memcpy ( &_arow[k][(i-blnum+2)*bldim+j0], &coeff[7+3-(j-j0)],
                 (j-j0+j1)*sizeof(double) );
      if ( i > blnum-2 )
        memcpy ( &_arow[k][(i-blnum+1)*bldim+j0], &coeff[14+3-(j-j0)],
                 (j-j0+j1)*sizeof(double) );
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
} /*g2bl_SetupClosedTriharmAMatrixd*/

boolean g2bl_SetupClosedTriharmRHSd ( int lastknotu, int lastknotv,
                           int spdimen, int pitch, const double *cpoints,
                           double *rhs )
{
  int bldim, blnum, n;
  int i, j, k, l;

  bldim = lastknotv-9;
  blnum = lastknotu-6;
  if ( bldim < 1 || blnum < 7 ) {
    pkv_SignalError ( LIB_G2BLENDING, 3, ERRMSG_2 );
    return false;
  }

  n = bldim*blnum;
  memset ( rhs, 0, spdimen*n*sizeof(double) );
        /* for each column of the control net */
  for ( i = 0; i < blnum; i++ ) {
    for ( j = 0;  j < 7;  j++ ) {
      k = 3 - abs ( j-3 );
      l = i - (j-3);
      if ( l < 0 )           l += blnum;
      else if ( l >= blnum ) l -= blnum;
      pkn_MultMatrixSubd ( 1, 3, 0, &coeff[7*k], spdimen, spdimen,
                 &cpoints[i*pitch], 0, &rhs[l*bldim*spdimen] );
      if ( bldim > 1 )
        pkn_MultMatrixSubd ( 1, 2, 0, &coeff[7*k], spdimen, spdimen,
                   &cpoints[i*pitch+spdimen], 0, &rhs[(l*bldim+1)*spdimen] );
      if ( bldim > 2 )
        pkn_MultMatrixSubd ( 1, 1, 0, &coeff[7*k], spdimen, spdimen,
                   &cpoints[i*pitch+2*spdimen], 0, &rhs[(l*bldim+2)*spdimen] );
      pkn_MultMatrixSubd ( 1, 3, 0, &coeff[7*k+4], spdimen, spdimen,
                   &cpoints[i*pitch+(lastknotv-6)*spdimen], 0,
                   &rhs[((l+1)*bldim-1)*spdimen] );
      if ( bldim > 1 )
        pkn_MultMatrixSubd ( 1, 2, 0, &coeff[7*k+5], spdimen, spdimen,
                   &cpoints[i*pitch+(lastknotv-6)*spdimen], 0,
                   &rhs[((l+1)*bldim-2)*spdimen] );
      if ( bldim > 2 )
        pkn_MultMatrixSubd ( 1, 1, 0, &coeff[7*k+6], spdimen, spdimen,
                   &cpoints[i*pitch+(lastknotv-6)*spdimen], 0,
                   &rhs[((l+1)*bldim-3)*spdimen] );
    }
  }
  return true;
} /*g2bl_SetupClosedTriharmRHSd*/

