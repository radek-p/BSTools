
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2007, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

boolean mbs_ClosedKnotsCorrectf ( int degree, int lastknot, float *knots,
                                  float T, int K, float tol )
{
  int   i, m;

  /* are they sorted? */
  for ( i = 0; i < lastknot-1; i++ )
    if ( knots[i] > knots[i+1] ) {
      /* printf ( "wrong knot order: %d\n", i ); */
      return false;
    }
  /* are they periodic? */
  for ( i = 1; i+K < lastknot; i++ )
    if ( fabs(knots[i+K]-knots[i]-T) > tol ) {
      /* printf ( "wrong knot periodicity: %d\n", i ); */
      return false;
    }
  /* are their multiplicities not too high? */
  for ( i = 1; i < lastknot; i += m ) {
    for ( m = 1;  i+m < lastknot && knots[i+m]-knots[i] < tol;  m++ )
      ;
    if ( m > degree )
      return false;
  }    
  return true;
} /*mbs_ClosedKnotsCorrectf*/

