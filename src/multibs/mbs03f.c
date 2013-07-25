
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* changes: 24.08.2004: P.Kiciak - fixed error in mbs_FindKnotIntervalf  */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
int mbs_KnotMultiplicityf ( int lastknot, const float *knots, float t )
{
  int k, l, r;

                       /* binary search */
  k = 0;  l = lastknot;
  while ( l-k > 1 ) {
    r = (k+l)/2;
    if ( t > knots[r] )
      k = r;
    else
      l = r;
  }
  while ( k < lastknot && t == knots[k+1] )
    k++;
  if ( t != knots[k] )
    k--;
  r = 0;
  while ( k >= 0 && t == knots[k] ) {
    k--;
    r++;
  }
  return r;
} /*mbs_KnotMultiplicityf*/

int mbs_FindKnotIntervalf ( int degree, int lastknot, const float *knots,
                            float t, int *mult )
{
  int k, l, r;

    /* if degree < 1 then the purpose of this procedure is just to find */
    /* the proper interval between the knots; it may be -1 or lastknot  */
    /* if t is less than or greater than all knots */

  if ( degree < 0 ) {
    if ( t < knots[0] ) {
      if ( mult ) *mult = 0;
      return -1;
    }
    else if ( t > knots[lastknot] ) {
      if ( mult ) *mult = 0;
      return lastknot;
    }
    else if ( t == knots[lastknot] ) {
      if ( mult ) {
        k = lastknot-1;  r = 1;
        while ( t == knots[k] ) { r++; k--; }
        *mult = r;
      }
      return lastknot;
    }
    degree = 0;
  }

    /* if degree >= 0 then the purpose of this procedure is to find   */
    /* the domain of the polynomial arc, whose element is t. The arcs */
    /* are numbered from degree to lastknot-degree-1; in that case    */
    /* the parameter mult is the multiplicity of the left end of the  */
    /* interval found, if t is equal to it, or else 0. */

  if ( t < knots[degree] ) {
    if ( mult ) *mult = 0;
    k = degree;  t = knots[k];
    while ( t == knots[k+1] ) k++;
    return k;
  }
  else if ( t > knots[lastknot-degree] ) {
    if ( mult ) *mult = 0;
    k = lastknot-degree-1;  t = knots[lastknot-degree];
    while ( t == knots[k] ) k--;
    return k;
  }
  else {
                       /* binary search */
    k = degree;  l = lastknot-degree;
    while ( l-k > 1 ) {
      r = (k+l)/2;
      if ( t > knots[r] )
        k = r;
      else
        l = r;
    }
    
    while ( k < lastknot-degree-1 && t == knots[k+1] )
      k++;
    if ( mult ) {
      r = 0;
      while ( k-r > 0 && t == knots[k-r] )
        r++;
      *mult = r;
    }
    return k;
  }
} /*mbs_FindKnotIntervalf*/

