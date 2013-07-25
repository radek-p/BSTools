
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2011                            */
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
boolean pkn_PCGf ( int n, void *usrdata, float *b, float *x,
            boolean (*multAx)( int n, void *usrdata, const float *x, float *Ax ),
            boolean (*multQIx)( int n, void *usrdata, const float *x, float *Qix ),
            int maxit, float eps, float delta, int *itm )
{
  void  *sp;
  int   k;
  float *r, *z, *v, c, d, t, rn;

  sp = pkv_GetScratchMemTop ();
  r = pkv_GetScratchMemf ( n );
  z = pkv_GetScratchMemf ( n );
  v = pkv_GetScratchMemf ( n );
  k = -1;
  if ( !r || !z || !v )
    goto failure;

        /* the conjugate gradient method */
  if ( !multAx ( n, usrdata, x, z ) )
    goto failure;
  pkn_SubtractMatrixf ( 1, n, 0, b, 0, z, 0, r );
  rn = pkn_ScalarProductf ( n, r, r );
  if ( multQIx ) {    /* preconditioning */
    if ( !multQIx ( n, usrdata, r, z ) )
      goto failure;
    c = pkn_ScalarProductf ( n, z, r );
  }
  else {              /* no preconditioning */
    memcpy ( z, r, n*sizeof(float) );
    c = rn;
  }
  memcpy ( v, z, n*sizeof(float) );
  for ( k = 0; k < maxit; k++ ) {
    d = pkn_ScalarProductf ( n, v, v );
    if ( d < delta )
      goto success;
    if ( !multAx ( n, usrdata, v, z ) )
      goto failure;
    t = c/pkn_ScalarProductf ( n, v, z );
    pkn_AddMatrixMf ( 1, n, 0, x, 0, v, t, 0, x );
    pkn_AddMatrixMf ( 1, n, 0, r, 0, z, -t, 0, r );
    t = pkn_ScalarProductf ( n, r, r );
    if ( multQIx ) {    /* preconditioning */
      if ( !multQIx ( n, usrdata, r, z ) )
        goto failure;
      d = pkn_ScalarProductf ( n, z, r );
    }
    else {              /* no preconditioning */
      memcpy ( z, r, n*sizeof(float) );
      d = t;
    }
    if ( d < eps && t < eps )
      goto success;
    pkn_AddMatrixMf ( 1, n, 0, z, 0, v, d/c, 0, v );
    c = d;
    rn = t;
  }

success:
  if ( itm )
    *itm = k;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  if ( itm )
    *itm = k;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*pkn_PCGf*/

