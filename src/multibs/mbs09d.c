
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

int mbs_NumKnotIntervalsd ( int degree, int lastknot, const double *knots )
{
  int i, k;

  k = 0;
  for ( i = degree+1; i <= lastknot-degree; i++ )
    if ( knots[i] > knots[i-1] )
      k++;
  return k;
} /*mbs_NumKnotIntervalsd*/

int mbs_BSProdRepSized ( int degree1, int lastknot1, const double *knots1,
                         int degree2, int lastknot2, const double *knots2 )
{
  int i, j, k, nm, mult_u, mult_v, mult_z;
  double za, zb, ckn;

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
} /*mbs_BSProdRepSized*/

int mbs_NumMaxKnotsd ( int degree, int lastknot, const double *knots )
{
  int i, k;

  for ( i = degree+1, k = 0; i <= lastknot-degree-1; i++ )
    if ( knots[i] != knots[i-1] )
      k++;
  return (k+2)*(degree+1);
} /*mbs_NumMaxKnotsd*/

void mbs_SetBSProdKnotsd ( int degree1, int lastknot1, const double *knots1,
                           int degree2, int lastknot2, const double *knots2,
                           int *degree, int *lastknot, double *knots )
{
  int i, j, k, mult_u, mult_v, mult_z, nm;
  double za, zb, ckn;

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
} /*mbs_SetBSProdKnotsd*/

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
void mbs_SetKnotPatternd ( int lastinknot, const double *inknots,
                           int multipl,
                           int *lastoutknot, double *outknots )
{
  int   i, j, k;
  double u;

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
} /*mbs_SetKnotPatternd*/



/* /////////////////////////////////////////// */
/* conversion between scaled and unscaled Bernstein bases */
void mbs_multiBezScaled ( int degree, int narcs, int ncurves, int spdimen,
                          int pitch, double *ctlpoints )
{
  int    i, j, k, l;
  double *a, *b, *c;

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
} /*mbs_multiBezScaled*/

void mbs_multiBezUnscaled ( int degree, int narcs, int ncurves, int spdimen,
                            int pitch, double *ctlpoints )
{
  int    i, j, k, l;
  double *a, *b, *c, binv;

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
} /*mbs_multiBezUnscaled*/

