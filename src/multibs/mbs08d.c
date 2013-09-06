
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
/* The procedure below removes knots whose multiplicity exceeds  */
/* degree+1, which is done by moving data in the arrays.         */
void mbs_multiRemoveSuperfluousKnotsd ( int ncurves, int spdimen, int degree,
                                        int *lastknot, double *knots,
                                        int inpitch, int outpitch,
                                        double *ctlpoints )
{
  int   i, j, k, r, l, lkn, nlkn;
  double u0, uN, ta, tb;

  nlkn = lkn = *lastknot;
  u0 = knots[0];  uN = knots[lkn];
  knots[0] = knots[1];  knots[lkn] = knots[lkn-1];

  i = j = 0;
  while ( j+degree+1 <= lkn ) {
    ta = knots[j];
    if ( knots[j+degree+1] == ta ) {
      for ( r = 1; j+degree+r < lkn && knots[j+degree+r+1] == ta; r++ ) ;
      nlkn -= r;
      j += r;
      for ( k = j+degree+1, l = degree+1;  k <= lkn-degree;  k++, l++ ) {
        tb = knots[k];
        if ( knots[k+degree+1] == tb )
          break;
      }
      if ( k >= lkn-degree ) {
        memmove ( &knots[i], &knots[j], (lkn-j+1)*sizeof(double) );
        l = lkn-j-degree;
        if ( l > 0 )
          pkv_Moved ( ncurves, l*spdimen, inpitch, (i-j)*spdimen,
                      &ctlpoints[j*spdimen] );
        break;
      }
      memmove ( &knots[i], &knots[j], l*sizeof(double) );
      pkv_Moved ( ncurves, l*spdimen, inpitch, (i-j)*spdimen,
                  &ctlpoints[j*spdimen] );
      i += l;
      j = k;
    }
    else {
      i++;
      j++;
    }
  }

  *lastknot = nlkn;
  knots[0] = u0;  knots[nlkn] = uN;
  if ( outpitch != inpitch )
    pkv_Rearranged ( ncurves, (*lastknot-degree)*spdimen,
                     inpitch, outpitch, ctlpoints );
} /*mbs_multiRemoveSuperfluousKnotsd*/

/* *** the old and naive version *** */
/*
void mbs_multiRemoveSuperfluousKnotsd ( int ncurves, int spdimen, int degree,
                                        int *lastknot, double *knots,
                                        int inpitch, int outpitch,
                                        double *ctlpoints )
{
  int i;

  while ( knots[*lastknot-degree-1] == knots[*lastknot-1] )
    (*lastknot) --;

  for ( i = (*lastknot)-degree-2; i > 0; i-- ) {
    if ( knots[i] == knots[i+degree+1] )
      mbs_multiKnotRemoved ( degree, lastknot, knots, ncurves,
               spdimen, inpitch, inpitch, ctlpoints, i );
  }
  if ( knots[degree+1] == knots[1] ) {
    (*lastknot) --;
    memmove ( &knots[1], &knots[2], (*lastknot-1)*sizeof(double) );
    pkv_Moved ( ncurves, (*lastknot-degree)*spdimen, inpitch,
                -spdimen, &ctlpoints[spdimen] );
  }
  if ( outpitch != inpitch )
    pkv_Rearranged ( ncurves, (*lastknot-degree)*spdimen,
                     inpitch, outpitch, ctlpoints );
}*/ /*mbs_multiRemoveSuperfluousKnotsd*/


/* /////////////////////////////////////////// */
/* computing the number of knots after the maximal knot insertion */
int mbs_LastknotMaxInsd ( int degree, int lastknot, const double *knots,
                          int *numknotintervals )
{
  int nki, skip, i;

  nki = mbs_NumKnotIntervalsd ( degree, lastknot, knots );
  if ( numknotintervals )
    *numknotintervals = nki;
  skip = 0;
  for ( i = 1; knots[i] != knots[degree]; i++ )
    skip++;
  for ( i = lastknot-1; knots[i] != knots[lastknot-degree]; i-- )
    skip++;
  return (degree+1)*(nki+1) + skip - 1;
} /*mbs_LastknotMaxInsd*/

/* /////////////////////////////////////////// */
/* maximal knot insertion - conversion to piecewise Bezier form */

void mbs_multiMaxKnotInsd ( int ncurves, int spdimen, int degree,
                            int inlastknot, const double *inknots,
                            int inpitch, const double *inctlpoints,
                            int *outlastknot, double *outknots,
                            int outpitch, double *outctlpoints,
                            int *skipl, int *skipr )
{
  int    i, j, k, l, r, pitch, lkn;
  double *auxkn, *auxcp;
  void   *stp;

/* In order to leave the input data intact, we have to make a working */
/* copy of the data. Unfortunately it is possible that the output     */
/* representation takes less space than input data, so we cannot just */
/* copy input to output arrays.                                       */
  stp = pkv_GetScratchMemTop ();
  lkn = inlastknot;
  auxkn = pkv_GetScratchMemd ( lkn+1 );
  pitch = spdimen*(lkn-degree);
  auxcp = pkv_GetScratchMemd ( pitch*ncurves );
  if ( !auxkn || !auxcp ) {
    PKV_SIGNALERROR ( LIB_MULTIBS, ERRCODE_2, ERRMSG_2 );
    exit ( 1 );
  }
  memcpy ( auxkn, inknots, (lkn+1)*sizeof(double) );
  pkv_Selectd ( ncurves, pitch, inpitch, pitch, inctlpoints, auxcp );
  mbs_multiRemoveSuperfluousKnotsd ( ncurves, spdimen, degree, &lkn,
                                     auxkn, pitch, pitch, auxcp );

/* Find the final knot sequence. */
  for ( j = 0; auxkn[j] < auxkn[degree]; j++ )
    outknots[j] = auxkn[j];
  r = degree;  if ( auxkn[0] == auxkn[degree] ) r++;
  for ( i = 0; i < r; i++, j++ )
    outknots[j] = auxkn[degree];
  while ( auxkn[i] == auxkn[degree] )
    i++;
  k = lkn-degree-1;
  while ( auxkn[k] == auxkn[lkn-degree] )
    k--;
  if ( k >= i ) {
    mbs_SetKnotPatternd ( k-i, &auxkn[i], degree+1, &l, &outknots[j] );
    j += l+1;
  }
  for ( i = 0; i < degree; i++, j++ )
    outknots[j] = auxkn[lkn-degree];
  i = lkn-degree+1;
  while ( i < lkn && auxkn[i] == auxkn[lkn-degree] ) i++;
  for ( ; i <= lkn; i++, j++ )
    outknots[j] = auxkn[i];
  j--;

  if ( skipl ) {
    for ( i = 1, r = 0; outknots[i] < outknots[degree]; i++, r++ ) ;
    *skipl = r;
  }
  if ( skipr ) {
    for ( i = j-1, r = 0; outknots[i] > outknots[j-degree]; i--, r++ ) ;
    *skipr = r;
  }

  *outlastknot = j;
  if ( j == lkn )  /* no knot is to be inserted - just copy data */
    pkv_Selectd ( ncurves, spdimen*(lkn-degree), pitch, outpitch,
                  auxcp, outctlpoints );

  else
/* The last thing - multiple knot insertion, with use of the Oslo algorithm. */
/* We do not call mbs_OsloKnotsCorrectf, as we have just created the correct */
/* knot sequences. */
    mbs_multiOsloInsertKnotsd ( ncurves, spdimen, degree, lkn, auxkn, pitch,
                                auxcp, j, outknots, outpitch, outctlpoints );

/* Cleanup. */
  pkv_SetScratchMemTop ( stp );
} /*mbs_multiMaxKnotInsd*/

