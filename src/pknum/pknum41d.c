
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
boolean pkn_PCGd ( int n, void *usrdata, double *b, double *x,
            boolean (*multAx)( int n, void *usrdata, const double *x, double *Ax ),
            boolean (*multQIx)( int n, void *usrdata, const double *x, double *Qix ),
            int maxit, double eps, double delta, int *itm )
{
  void  *sp;
  int   k;
  double *r, *z, *v, c, d, t, rn;

  sp = pkv_GetScratchMemTop ();
  r = pkv_GetScratchMemd ( n );
  z = pkv_GetScratchMemd ( n );
  v = pkv_GetScratchMemd ( n );
  k = -1;
  if ( !r || !z || !v )
    goto failure;

        /* the conjugate gradient method */
  if ( !multAx ( n, usrdata, x, z ) )
    goto failure;
  pkn_SubtractMatrixd ( 1, n, 0, b, 0, z, 0, r );
  rn = pkn_ScalarProductd ( n, r, r );
  if ( multQIx ) {    /* preconditioning */
    if ( !multQIx ( n, usrdata, r, z ) )
      goto failure;
    c = pkn_ScalarProductd ( n, z, r );
  }
  else {              /* no preconditioning */
    memcpy ( z, r, n*sizeof(double) );
    c = rn;
  }
  memcpy ( v, z, n*sizeof(double) );
  for ( k = 0; k < maxit; k++ ) {
    d = pkn_ScalarProductd ( n, v, v );
    if ( d < delta )
      goto success;
    if ( !multAx ( n, usrdata, v, z ) )
      goto failure;
    t = c/pkn_ScalarProductd ( n, v, z );
    pkn_AddMatrixMd ( 1, n, 0, x, 0, v, t, 0, x );
    pkn_AddMatrixMd ( 1, n, 0, r, 0, z, -t, 0, r );
    t = pkn_ScalarProductd ( n, r, r );
    if ( multQIx ) {    /* preconditioning */
      if ( !multQIx ( n, usrdata, r, z ) )
        goto failure;
      d = pkn_ScalarProductd ( n, z, r );
    }
    else {              /* no preconditioning */
      memcpy ( z, r, n*sizeof(double) );
      d = t;
    }
    if ( d < eps && t < eps )
      goto success;
    pkn_AddMatrixMd ( 1, n, 0, z, 0, v, d/c, 0, v );
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
} /*pkn_PCGd*/

