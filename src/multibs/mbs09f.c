
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"


/* /////////////////////////////////////////// */
/* spline algebra - processing knot sequences */

int mbs_NumKnotIntervalsf ( int degree, int lastknot, const float *knots )
{
  int i, k;

  k = 0;
  for ( i = degree+1; i <= lastknot-degree; i++ )
    if ( knots[i] > knots[i-1] )
      k++;
  return k;
} /*mbs_NumKnotIntervalsf*/

int mbs_BSProdRepSizef ( int degree1, int lastknot1, const float *knots1,
                         int degree2, int lastknot2, const float *knots2 )
{
  int i, j, k, nm, mult_u, mult_v, mult_z;
  float za, zb, ckn;

  nm = degree1+degree2;

  za = knots1[degree1] > knots2[degree2] ? knots1[degree1] : knots2[degree2];
  zb = knots1[lastknot1-degree1] > knots2[lastknot2-degree2] ?
         knots1[lastknot1-degree1] : knots2[lastknot2-degree2];
  if ( za >= zb )
    return 0;

  k = 2*nm+1;
  i = degree1;  while ( knots1[i] <= za ) i++;
  j = degree2;  while ( knots2[j] <= za ) j++;
  while ( knots1[i] < zb || knots2[j] < zb ) {
    ckn = knots1[i] < knots2[j] ? knots1[i] : knots2[j];
    mult_u = 0;  while ( knots1[i] == ckn ) { mult_u++;  i++; }
    mult_v = 0;  while ( knots2[j] == ckn ) { mult_v++;  j++; }
    if ( mult_u == 0 )
      mult_z = degree1+mult_v;
    else if ( mult_v == 0 )
      mult_z = degree2+mult_u;
    else
      mult_z = degree2+mult_u > degree1+mult_v ?
               degree2+mult_u : degree1+mult_v;
    k += mult_z;
  }
  return k;
} /*mbs_BSProdRepSizef*/

int mbs_NumMaxKnotsf ( int degree, int lastknot, const float *knots )
{
  int i, k;

  for ( i = degree+1, k = 0; i <= lastknot-degree-1; i++ )
    if ( knots[i] != knots[i-1] )
      k++;
  return (k+2)*(degree+1);
} /*mbs_NumMaxKnotsf*/

void mbs_SetBSProdKnotsf ( int degree1, int lastknot1, const float *knots1,
                           int degree2, int lastknot2, const float *knots2,
                           int *degree, int *lastknot, float *knots )
{
  int i, j, k, mult_u, mult_v, mult_z, nm;
  float za, zb, ckn;

  nm = degree1+degree2;

  za = knots1[degree1] > knots2[degree2] ? knots1[degree1] : knots2[degree2];
  zb = knots1[lastknot1-degree1] > knots2[lastknot2-degree2] ?
         knots1[lastknot1-degree1] : knots2[lastknot2-degree2];
  if ( za >= zb ) {
    *lastknot = 0;
    return;
  }
  knots[0] = knots1[0] > knots2[0] ? knots1[0] : knots2[0];
  k = 1;
  mult_z = nm;
  while ( mult_z > 0 ) {
    knots[k] = za;
    k++;
    mult_z--;
  }
  i = degree1;  while ( knots1[i] <= za ) i++;
  j = degree2;  while ( knots2[j] <= za ) j++;
  while ( knots1[i] < zb || knots2[j] < zb ) {
    ckn = knots1[i] < knots2[j] ? knots1[i] : knots2[j];
    mult_u = 0;  while ( knots1[i] == ckn ) { mult_u++;  i++; }
    mult_v = 0;  while ( knots2[j] == ckn ) { mult_v++;  j++; }
    if ( mult_u == 0 )
      mult_z = degree1+mult_v;
    else if ( mult_v == 0 )
      mult_z =  degree2+mult_u;
    else {
      mult_z = degree2+mult_u > degree1+mult_v ?
               degree2+mult_u : degree1+mult_v;
      mult_z = nm+1 < mult_z ? nm+1 : mult_z;
    }
    while ( mult_z > 0 ) {
      knots[k] = ckn;
      k++;
      mult_z--;
    }
  }
  mult_z = nm;
  while ( mult_z > 0 ) {
    knots[k] = zb;
    k++;
    mult_z--;
  }
  knots[k] = knots1[lastknot1] < knots2[lastknot2] ?
             knots1[lastknot1] : knots2[lastknot2];
  *lastknot = k;
  *degree = nm;
} /*mbs_SetBSProdKnotsf*/

/* The procedure below creates a regular knot pattern, i.e. a sequence of */
/* knots of the same multiplicity, specified by the parameter multipl.    */
/* The knots are taken from the sequence of length lastinknot+1, given in */
/* the array inknots. The procedure is intended for the preparation of a  */
/* piecewise Bezier representation of a curve or a patch.                 */
/* To allocate the array for the output knot sequence you can compute its */
/* length by calling the mbs_NumKnotIntervalsf procedure. Add 1 to the    */
/* result and multiply by degree+1. This works if                         */
/* inknots[1]=...=inknots[degree] and                                     */
/*  inknots[lastinknot-degree]=...=inknot[lastinknot-1].                  */
void mbs_SetKnotPatternf ( int lastinknot, const float *inknots,
                           int multipl,
                           int *lastoutknot, float *outknots )
{
  int   i, j, k;
  float u;

  u = inknots[lastinknot];
  if ( u != inknots[0] ) {
    for ( i = k = 0; i <= lastinknot; i++ ) {
      if ( u != inknots[i] ) {
        u = inknots[i];
        for ( j = 0; j < multipl; j++, k++ )
          outknots[k] = u;
      }
    }
  }
  else {
    for ( k = 0; k < multipl; k++ )
      outknots[k] = u;
  }
  *lastoutknot = k-1;
} /*mbs_SetKnotPatternf*/



/* /////////////////////////////////////////// */
/* conversion between scaled and unscaled Bernstein bases */
void mbs_multiBezScalef ( int degree, int narcs, int ncurves, int spdimen,
                          int pitch, float *ctlpoints )
{
  int   i, j, k, l;
  float *a, *b, *c;

  if ( degree < 2 )
    return;
  else if ( degree <= 29 ) {
    int binom;

    binom = degree;
    for ( i = 1, a = &ctlpoints[spdimen];  i < degree;  i++, a += spdimen ) {
      for ( j = 0, b = a;  j < ncurves;  j++, b += pitch )
        for ( k = 0, c = b;  k < narcs;  k++, c += (degree+1)*spdimen )
          for ( l = 0;  l < spdimen;  l++ )
            c[l] *= (double)binom;
      binom = (binom*(degree-i))/(i+1);
    }
  }
  else if ( degree <= 61 ) {
          /* for high degrees 64 bit integers are needed to avoid */
          /* overflow in computing binomial coefficients */
          /* however, in single precision it may have not much sense */
#if __WORDSIZE == 64
    long int binom; 
#else
    long long int binom;
#endif

    binom = degree;
    for ( i = 1, a = &ctlpoints[spdimen];  i < degree;  i++, a += spdimen ) {
      for ( j = 0, b = a;  j < ncurves;  j++, b += pitch )
        for ( k = 0, c = b;  k < narcs;  k++, c += (degree+1)*spdimen )
          for ( l = 0;  l < spdimen;  l++ )
            c[l] *= (double)binom;
      binom = (binom*(degree-i))/(i+1);
    }
  }
} /*mbs_multiBezScalef*/

void mbs_multiBezUnscalef ( int degree, int narcs, int ncurves, int spdimen,
                            int pitch, float *ctlpoints )
{
  int    i, j, k, l;
  float  *a, *b, *c;
  double binv;

  if ( degree < 2 )
    return;
  else if ( degree <= 29 ) {
    int binom;

    binom = degree;
    for ( i = 1, a = &ctlpoints[spdimen];  i < degree;  i++, a += spdimen ) {
      binv = 1.0/(double)binom;
      for ( j = 0, b = a;  j < ncurves;  j++, b += pitch )
        for ( k = 0, c = b;  k < narcs;  k++, c += (degree+1)*spdimen )
          for ( l = 0;  l < spdimen;  l++ )
            c[l] *= binv;
      binom = (binom*(degree-i))/(i+1);
    }
  }
  else if ( degree <= 61 ) {
          /* for high degrees 64 bit integers are needed to avoid */
          /* overflow in computing binomial coefficients */
          /* however, in single precision it may have not much sense */
#if __WORDSIZE == 64
    long int binom; 
#else
    long long int binom;
#endif

    binom = degree;
    for ( i = 1, a = &ctlpoints[spdimen];  i < degree;  i++, a += spdimen ) {
      binv = 1.0/(double)binom;
      for ( j = 0, b = a;  j < ncurves;  j++, b += pitch )
        for ( k = 0, c = b;  k < narcs;  k++, c += (degree+1)*spdimen )
          for ( l = 0;  l < spdimen;  l++ )
            c[l] *= binv;
      binom = (binom*(degree-i))/(i+1);
    }
  }
} /*mbs_multiBezUnscalef*/

