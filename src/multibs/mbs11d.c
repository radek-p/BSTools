
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* /////////////////////////////////////////// */
/* degree elevation of Bezier curves */

boolean mbs_multiBCDegElevd ( int ncurves, int spdimen,
                              int inpitch, int indegree, const double *inctlpoints,
                              int deltadeg,
                              int outpitch, int *outdegree, double *outctlpoints )
{
  int   i, j, k, n, m;
  double *p, *q;

  if ( deltadeg < 0 ) {/* We can do the degree elevation, not reduction here. */
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    return false;
  }
    
/* If inctlpoints is NULL then we assume that data are already pointed by */
/* outctlpoints and we do not have to copy. But we may have to rearrange  */
/* the data if the final pitch is different than the initial pitch.       */
  if ( inctlpoints )
    pkv_Selectd ( ncurves, spdimen*(indegree+1), inpitch, outpitch,
                  inctlpoints, outctlpoints );
  else if ( outpitch != inpitch )
    pkv_Rearranged ( ncurves, spdimen*(indegree+1), inpitch, outpitch,
                     outctlpoints );

/* Now the proper degree elevation, by 1 with every loop execution. */
  for ( n = indegree; n < indegree+deltadeg; n++ ) {
    m = n+1;
    for ( k = 0, p = outctlpoints;  k < ncurves;  k++, p+= outpitch ) {
      q = p + n*spdimen;
      memcpy ( q+spdimen, q, spdimen*sizeof(double) );
      for ( i = n;  i > 0;  i--, q -= spdimen )
        for ( j = 0; j < spdimen; j++ )
          q[j] = ((double)i*q[j-spdimen]+(double)(m-i)*q[j])/(double)m;
    }
  }
  *outdegree = indegree+deltadeg;
  return true;
} /*mbs_multiBCDegElevd*/

/* /////////////////////////////////////////// */
/* degree elevation of B-spline curves */

/* It is assumed below that inknots[1]=...=inknots[indegree] and          */
/* inknots[(*inlastknot)-indegree]=...=inknots[(*inlastknot)-1], which is */
/* too restrictive and the procedure should be reimplemented.             */
/* It must be preceded by the appropriate reimplementation of the         */
/* procedure mbs_multiMaxKnotInsd, which assumes the same thing.          */

boolean mbs_multiBSDegElevd ( int ncurves, int spdimen,
                              int indegree, int inlastknot, const double *inknots,
                              int inpitch, const double *inctlpoints,
                              int deltadeg,
                              int *outdegree, int *outlastknot,
                              double *outknots, int outpitch, double *outctlpoints,
                              boolean freeend )
{
  void          *stp;
  int           auxpitch, ap, akns, skipl, skipr;
  int           i, j, k, l;
  int           ki, lkn, deg;
  double         *akn, *auxcp;
  double         t;

  if ( !deltadeg ) {  /* There is no degree elevation, but data may have to */
                      /* be rearranged or copied. */
    if ( inctlpoints )
      pkv_Selectd ( ncurves, spdimen*(inlastknot-indegree),
                    inpitch, outpitch, inctlpoints, outctlpoints );
    else
      pkv_Rearranged ( ncurves, spdimen*(inlastknot-indegree),
                       inpitch, outpitch, outctlpoints );
    *outdegree = indegree;
    *outlastknot = inlastknot;
    memcpy ( outknots, inknots, (inlastknot+1)*sizeof(double) );
    return true;
  }
  if ( deltadeg < 0 ) {/* We can do the degree elevation, not reduction here. */
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    return false;
  }

/* Create a working area. */
  deg = indegree+deltadeg;
  stp = pkv_GetScratchMemTop ();
  lkn = inlastknot;
  ap = (deg+1)*spdimen;
  akns = mbs_LastknotMaxInsd ( indegree, lkn, inknots, &ki );
  if ( ki <= 0 ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_5, ERRMSG_5 );
    goto failure;
  }
  auxpitch = max ( ki*ap, (akns-indegree)*spdimen );
  auxcp = pkv_GetScratchMemd ( ncurves*auxpitch + deg*spdimen );
  akn = pkv_GetScratchMemd ( (ki+3)*(deg+1) - 2 );

  if ( !auxcp || !akn ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    goto failure;
  }

/* Now the conversion to the piecewise Bezier form and then the actual */
/* degree elevation of the Bezier curve arcs, each time the loop is    */
/* executed, all arcs of one curve are processed.                      */
  if ( !inctlpoints )
    inctlpoints = outctlpoints;
  if ( !mbs_multiMaxKnotInsd ( ncurves, spdimen, indegree, inlastknot, inknots,
                               inpitch, inctlpoints, &lkn, akn, auxpitch, auxcp,
                               &skipl, &skipr ) )
    goto failure;
  auxcp += skipl*spdimen;
  lkn -= skipl+skipr;
  for ( k = 0; k < ncurves; k++ ) {
    if ( !mbs_multiBCDegElevd ( ki, spdimen, (indegree+1)*spdimen, indegree, NULL,
                                deltadeg, ap, outdegree, &auxcp[k*auxpitch] ) )
      goto failure;
  }
  
/* construct the final knot sequence - each knot has the multiplicity */
/* increased by deltadeg, in the array outknots, and the intermediate */
/* knot sequence in bkn - each knot of multiplicity equal to the      */
/* final degree plus one, as the part of the piecewise Bezier form.   */
  mbs_SetKnotPatternd ( inlastknot-2*indegree, &inknots[indegree],
                        deg+1, &lkn, akn );
  akn[0] = min ( inknots[0], akn[1] );
  akn[lkn] = max ( inknots[inlastknot], akn[lkn-1] );

  outknots[0] = akn[0];
  for ( j = 1; j <= deg; j++ )
    outknots[j] = inknots[indegree];
  i = indegree+1;
  while ( inknots[i] == inknots[indegree] )
    i++;
  while ( i < inlastknot-indegree &&
          inknots[i] < inknots[inlastknot-indegree] ) {
    t = inknots[i];
    if ( t != inknots[i-1] ) {
      for ( k = 0; k < deltadeg; k++, j++ )
        outknots[j] = t;
    }
    outknots[j] = inknots[i];
    i++;
    j++;
  }
  for ( k = 1; k <= deg; k++, j++ )
    outknots[j] = inknots[inlastknot-indegree];
  outknots[j] = akn[lkn];
  *outlastknot = j;

/* remove knots by calling the procedure, which solves */
/* the least-squares problem with the Oslo matrix. */

  mbs_multiOsloRemoveKnotsLSQd ( ncurves, spdimen, deg, lkn, akn, auxpitch,
                                 auxcp, j, outknots, outpitch, outctlpoints );

/* if these are free end points curves, modify appropriately */
/* their first and last knots */

  if ( freeend ) {
    pkv_SetScratchMemTop ( stp );
    akn = pkv_GetScratchMemd ( (*outdegree)+1 );

        /* generate the first outdegree+1 knots of the final knot sequence */
    for ( k = indegree, t = inknots[k];  t == inknots[k+1];  k++ );
    for ( i = *outdegree; i >= 0; i--, k-- ) {
      akn[i] = t;
      if ( k == 0 || inknots[k-1] != t ) {
        for ( j = 0; j < deltadeg; j++ )
          akn[--i] = t;
        if ( k > 0 )
          t = inknots[k-1];
        else
          break;
      }
    }
    mbs_multiBSChangeLeftKnotsd ( ncurves, spdimen, *outdegree, outknots,
                                  outpitch, outctlpoints, akn );

        /* generate the last outdegree+1 knots of the final knot sequence */
    for ( k = inlastknot-indegree, t = inknots[k];  t == inknots[k-1];  k-- );
    for ( i = l = (*outlastknot)-(*outdegree); i <= *outlastknot; i++, k++ ) {
      akn[i-l] = t;
      if ( k == inlastknot || inknots[k+1] != t ) {
        for ( j = 0; j < deltadeg; j++ )
          akn[(++i)-l] = t;
        if ( k < inlastknot )
          t = inknots[k+1];
        else
          break;
      }
    }
    mbs_multiBSChangeRightKnotsd ( ncurves, spdimen,
                                   *outdegree, *outlastknot, outknots,
                                   outpitch, (double*)outctlpoints, akn );
  }

  pkv_SetScratchMemTop ( stp );
  return true;

failure:
  pkv_SetScratchMemTop ( stp );
  return false;
} /*mbs_multiBSDegElevd*/

