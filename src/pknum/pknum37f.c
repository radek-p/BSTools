
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "pkvaria.h"

#undef CONST_
#define CONST_

#include "pknum.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_QuadGaussLegendre4f ( float a, float b, int n,
                                  float *qknots, float *qcoeff )
{
  int    i;
  float c, d, h, k1, k2;

  if ( n < 2 || (n & 0x01) )
    return false;
  d = b-a;
  if ( qknots ) {
    h = 1.0*d/(float)n;
    k1 = 1.0-1.0/SQRT3;
    k2 = 1.0+1.0/SQRT3;
    for ( i = 0; i < n; i += 2 ) {
      qknots[i] = a + h*((float)i+k1);
      qknots[i+1] = a + h*((float)i+k2);
    }
  }
  if ( qcoeff ) {
    c = d/(float)n;
    for ( i = 0; i < n; i++ )
      qcoeff[i] = c;
  }
  return true;
} /*pkn_QuadGaussLegendre4f*/

boolean pkn_QuadGaussLegendre6f ( float a, float b, int n,
                                  float *qknots, float *qcoeff )
{
  int    i, j;
  float kn[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
  float qc[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
  float c, d, h;

  if ( n < 3 || (n % 3) )
    return false;
  d = b-a;
  if ( qknots ) {
    h = d/(float)(n/3);
    for ( i = 0; i < n; i += 3 )
      for ( j = 0; j < 3; j++ )
        qknots[i+j] = a + h*((float)(i/3)+0.5+0.5*kn[j]);
  }
  if ( qcoeff ) {
    c = 1.5*d/(float)n;
    for ( i = 0; i < n; i += 3 )
      for ( j = 0; j < 3; j++ )
        qcoeff[i+j] = c*qc[j];
  }
  return true;
} /*pkn_QuadGaussLegendre6f*/

boolean pkn_QuadGaussLegendre8f ( float a, float b, int n,
                                  float *qknots, float *qcoeff )
{
  int    i, j;
  float kn[4] = {-0.8611363115940526, -0.3399810435848563,
                   0.3399810435848563,  0.8611363115940526};
  float qc[4] = {0.3478548451374539, 0.6521451548625461,
                  0.6521451548625461, 0.3478548451374539};
  float c, d, h;

  if ( n < 4 || (n % 4) )
    return false;
  d = b-a;
  if ( qknots ) {
    h = d/(float)(n/4);
    for ( i = 0; i < n; i += 4 )
      for ( j = 0; j < 4; j++ )
        qknots[i+j] = a + h*((float)(i/4)+0.5+0.5*kn[j]);
  }
  if ( qcoeff ) {
    c = 2.0*d/(float)n;
    for ( i = 0; i < n; i += 4 )
      for ( j = 0; j < 4; j++ )
        qcoeff[i+j] = c*qc[j];
  }
  return true;
} /*pkn_QuadGaussLegendre8f*/

boolean pkn_QuadGaussLegendre10f ( float a, float b, int n,
                                   float *qknots, float *qcoeff )
{
  int    i, j;
  float kn[5] = {-0.906179845938664, -0.5384693101056831, 0.0,
                   0.5384693101056831, 0.906179845938664};
  float qc[5] = {0.2369268850561891, 0.4786286704993665, 128.0/225.0,
                  0.4786286704993665, 0.2369268850561891};
  float c, d, h;

  if ( n < 5 || (n % 5) )
    return false;
  d = b-a;
  if ( qknots ) {
    h = d/(float)(n/5);
    for ( i = 0; i < n; i += 5 )
      for ( j = 0; j < 5; j++ )
        qknots[i+j] = a + h*((float)(i/5)+0.5+0.5*kn[j]);
  }
  if ( qcoeff ) {
    c = 2.5*d/(float)n;
    for ( i = 0; i < n; i += 5 )
      for ( j = 0; j < 5; j++ )
        qcoeff[i+j] = c*qc[j];
  }
  return true;
} /*pkn_QuadGaussLegendre10f*/

boolean pkn_QuadGaussLegendre12f ( float a, float b, int n,
                                   float *qknots, float *qcoeff )
{
  int    i, j;
  float kn[6] = {-0.932469514203152, -0.6612093864662645, -0.2386191860831969,
                   0.2386191860831969, 0.6612093864662645,  0.932469514203152};
  float qc[6] = {0.1713244923791703, 0.3607615730481386, 0.467913934572691,
                  0.467913934572691, 0.3607615730481386, 0.1713244923791703};
  float c, d, h;

  if ( n < 6 || (n % 6) )
    return false;
  d = b-a;
  if ( qknots ) {
    h = d/(float)(n/6);
    for ( i = 0; i < n; i += 6 )
      for ( j = 0; j < 6; j++ )
        qknots[i+j] = a + h*((float)(i/6)+0.5+0.5*kn[j]);
  }
  if ( qcoeff ) {
    c = 3.0*d/(float)n;
    for ( i = 0; i < n; i += 6 )
      for ( j = 0; j < 6; j++ )
        qcoeff[i+j] = c*qc[j];
  }
  return true;
} /*pkn_QuadGaussLegendre12f*/

boolean pkn_QuadGaussLegendre14f ( float a, float b, int n,
                                   float *qknots, float *qcoeff )
{
  int    i, j;
  float kn[7] = {-0.9491079123427585, -0.7415311855993944, -0.4058451513773972,
              0.0, 0.4058451513773972,  0.7415311855993944,  0.9491079123427585};
  float qc[7] = {0.1294849661688697,0.2797053914892767,0.3818300505051189,
     512.0/1225.0,0.3818300505051189,0.2797053914892767,0.1294849661688697};
  float c, d, h;

  if ( n < 7 || (n % 7) )
    return false;
  d = b-a;
  if ( qknots ) {
    h = d/(float)(n/7);
    for ( i = 0; i < n; i += 7 )
      for ( j = 0; j < 7; j++ )
        qknots[i+j] = a + h*((float)(i/7)+0.5+0.5*kn[j]);
  }
  if ( qcoeff ) {
    c = 3.5*d/(float)n;
    for ( i = 0; i < n; i += 7 )
      for ( j = 0; j < 7; j++ )
        qcoeff[i+j] = c*qc[j];
  }
  return true;
} /*pkn_QuadGaussLegendre14f*/

boolean pkn_QuadGaussLegendre16f ( float a, float b, int n,
                                   float *qknots, float *qcoeff )
{
  int    i, j;
  float kn[8] = {-0.9602898564975362, -0.7966664774136267, -0.525532409916329,
                  -0.1834346424956498,  0.1834346424956498,  0.525532409916329,
                   0.7966664774136267,  0.9602898564975362};
  float qc[8] = {0.1012285362903763,0.2223810344533745,0.3137066458778873,
                  0.362683783378362,0.362683783378362,0.3137066458778873,
                  0.2223810344533745,0.1012285362903763};
  float c, d, h;

  if ( n < 8 || (n % 8) )
    return false;
  d = b-a;
  if ( qknots ) {
    h = d/(float)(n/8);
    for ( i = 0; i < n; i += 8 )
      for ( j = 0; j < 8; j++ )
        qknots[i+j] = a + h*((float)(i/8)+0.5+0.5*kn[j]);
  }
  if ( qcoeff ) {
    c = 4.0*d/(float)n;
    for ( i = 0; i < n; i += 8 )
      for ( j = 0; j < 8; j++ )
        qcoeff[i+j] = c*qc[j];
  }
  return true;
} /*pkn_QuadGaussLegendre16f*/

boolean pkn_QuadGaussLegendre18f ( float a, float b, int n,
                                   float *qknots, float *qcoeff )
{
  int    i, j;
  float kn[9] = {-0.9681602395076261, -0.8360311073266358, -0.6133714327005904,
                  -0.3242534234038089, 0.0, 0.3242534234038089,
                  0.6133714327005904, 0.8360311073266358, 0.9681602395076261};
  float qc[9] = {0.08127438836157441, 0.1806481606948574, 0.2606106964029355,
                  0.3123470770400028, 32768.0/99225.0,  0.3123470770400028,
                  0.2606106964029355, 0.1806481606948574, 0.08127438836157441};
  float c, d, h;

  if ( n < 9 || (n % 9) )
    return false;
  d = b-a;
  if ( qknots ) {
    h = d/(float)(n/9);
    for ( i = 0; i < n; i += 9 )
      for ( j = 0; j < 9; j++ )
        qknots[i+j] = a + h*((float)(i/9)+0.5+0.5*kn[j]);
  }
  if ( qcoeff ) {
    c = 4.5*d/(float)n;
    for ( i = 0; i < n; i += 9 )
      for ( j = 0; j < 9; j++ )
        qcoeff[i+j] = c*qc[j];
  }
  return true;
} /*pkn_QuadGaussLegendre18f*/

boolean pkn_QuadGaussLegendre20f ( float a, float b, int n,
                                   float *qknots, float *qcoeff )
{
  int    i, j;
  float kn[10] = {-0.9739065285171717, -0.8650633666889845, -0.6794095682990244,
                   -0.4333953941292472, -0.1488743389816312, 0.1488743389816312,
                   0.4333953941292472, 0.6794095682990244, 0.8650633666889845,
                   0.9739065285171717};
  float qc[10] = {0.06667134430868814, 0.1494513491505806, 0.219086362515982,
                   0.2692667193099964, 0.2955242247147529, 0.2955242247147529,
                   0.2692667193099964, 0.219086362515982, 0.1494513491505806,
                   0.06667134430868814};
  float c, d, h;

  if ( n < 10 || (n % 10) )
    return false;
  d = b-a;
  if ( qknots ) {
    h = d/(float)(n/10);
    for ( i = 0; i < n; i += 10 )
      for ( j = 0; j < 10; j++ )
        qknots[i+j] = a + h*((float)(i/10)+0.5+0.5*kn[j]);
  }
  if ( qcoeff ) {
    c = 5.0*d/(float)n;
    for ( i = 0; i < n; i += 10 )
      for ( j = 0; j < 10; j++ )
        qcoeff[i+j] = c*qc[j];
  }
  return true;
} /*pkn_QuadGaussLegendre20f*/

