
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>
#include <stdio.h>  /* for testing */

#include "pkvaria.h"
#include "pkgeom.h"
#include "multibs.h"

/* workspace length must be at least (3*(degreev+1)+36)*spdimen doubles */
boolean _mbs_BCHornerDer2Pd ( int degreeu, int degreev, int spdimen,   
                              const double *ctlpoints,   
                              double u, double v,   
                              double *p, double *du, double *dv,
                              double *duu, double *duv, double *dvv,
                              double *workspace )
{
  double *q, *r;
  int    n, m, i, k, pitch;

  pitch = (degreev+1)*spdimen;
  r = workspace;
  q = &r[36*spdimen];
  if ( degreeu <= 2 ) {
    n = degreeu+1;
    memcpy ( q, ctlpoints, n*pitch*sizeof(double) );
  }
  else {
    n = 3;
    if ( !mbs_multiBCHornerd ( degreeu-2, 3, pitch, pitch, ctlpoints, u, q ) )
      return false;
  }
  if ( degreev <= 2 ) {
    m = degreev+1;
    pkv_Selectd ( n, m*spdimen, pitch, 3*spdimen, q, r );
  }
  else {
    m = 3;
    for ( i = 0; i < n; i++ )
      if ( !mbs_multiBCHornerd ( degreev-2, 3, spdimen, spdimen,
                                 &q[i*pitch], v, &r[i*3*spdimen] ) )
        return false;
  }
                    /* now the last steps of the de Casteljau algorithm */
  switch ( n ) {
case 3:
    pkn_MatrixLinCombd ( 2, m*spdimen, 3*spdimen, r, 1.0-u,
          3*spdimen, &r[3*spdimen], u, 3*spdimen, &r[9*spdimen] );
    pkn_MatrixLinCombd ( 1, m*spdimen, 0, &r[9*spdimen], 1.0-u,
          0, &r[12*spdimen], u, 0, &r[15*spdimen] );
    k = 6;
    break;
case 2:
    memcpy ( &r[9*spdimen], r, 6*spdimen*sizeof(double) );
    pkn_MatrixLinCombd ( 1, m*spdimen, 0, &r[9*spdimen], 1.0-u,
          0, &r[12*spdimen], u, 0, &r[15*spdimen] );
    k = 3;
    break;
case 1:
    memcpy ( &r[15*spdimen], r, 3*spdimen*sizeof(double) );
    k = 1;
    break;
default:
    k = 0;  /* to suppress a warning */
  }
  switch ( m ) {
case 3:
    pkn_MatrixLinCombd ( k, spdimen, 3*spdimen,
          &r[3*(6-k)*spdimen], 1.0-v,
          3*spdimen, &r[(3*(6-k)+1)*spdimen], v,
          2*spdimen, &r[(2*(6-k)+18)*spdimen] );
    pkn_MatrixLinCombd ( k, spdimen, 3*spdimen,
          &r[(3*(6-k)+1)*spdimen], 1.0-v,
          3*spdimen, &r[(3*(6-k)+2)*spdimen], v,
          2*spdimen, &r[(2*(6-k)+19)*spdimen] );
    pkn_MatrixLinCombd ( k, spdimen, 2*spdimen,
          &r[(2*(6-k)+18)*spdimen], 1.0-v,
          2*spdimen, &r[(2*(6-k)+19)*spdimen], v,
          spdimen, &r[(6-k+30)*spdimen] );
    break;
case 2:
    pkv_Selectd ( k, 2*spdimen, 3*spdimen, 2*spdimen,
                  &r[(3*(6-k))*spdimen], &r[(2*(6-k)+18)*spdimen] );
    pkn_MatrixLinCombd ( k, spdimen, 2*spdimen,
          &r[(2*(6-k)+18)*spdimen], 1.0-v,
          2*spdimen, &r[(2*(6-k)+19)*spdimen], v,
          spdimen, &r[(6-k+30)*spdimen] );
    break;
case 1:
    pkv_Selectd ( k, spdimen, 3*spdimen, spdimen,
                  &r[(3*(6-k))*spdimen], &r[(6-k+30)*spdimen] );
    break;
  }
                    /* compute the patch point and derivatives */
  memcpy ( p, &r[35*spdimen], spdimen*sizeof(double) );
  if ( degreeu > 0 ) {
    pkn_MatrixMDifferenced ( 1, spdimen, 0, &r[34*spdimen],
                             0, &r[33*spdimen], (double)degreeu, 0, du );
    if ( degreeu > 1 ) {
      pkn_AddMatrixd ( 1, spdimen, 0, &r[30*spdimen], 0, &r[32*spdimen], 0, duu );
      pkn_AddMatrixMd ( 1, spdimen, 0, duu, 0, &r[31*spdimen], -2.0, 0, duu );
      pkn_MultMatrixNumd ( 1, spdimen, 0, duu, (double)(degreeu*(degreeu-1)),
                           0, duu );
    }
    else
      memset ( duu, 0, spdimen*sizeof(double) );
    if ( degreev > 0 ) {
      pkn_AddMatrixd ( 1, spdimen, 0, &r[24*spdimen], 0, &r[27*spdimen], 0, duv );
      pkn_SubtractMatrixd ( 1, spdimen, 0, duv, 0, &r[25*spdimen], 0, duv );
      pkn_SubtractMatrixd ( 1, spdimen, 0, duv, 0, &r[26*spdimen], 0, duv );
      pkn_MultMatrixNumd ( 1, spdimen, 0, duv, (double)(degreeu*degreev), 0, duv );
    }
    else
      memset ( duv, 0, spdimen*sizeof(double) );
  }
  else {
    memset ( du, 0, spdimen*sizeof(double) );
    memset ( duu, 0, spdimen*sizeof(double) );
    memset ( duv, 0, spdimen*sizeof(double) );
  }
  if ( degreev > 0 ) {
    pkn_MatrixMDifferenced ( 1, spdimen, 0, &r[29*spdimen],
                             0, &r[28*spdimen], (double)degreev, 0, dv );
    if ( degreev > 1 ) {
      pkn_AddMatrixd ( 1, spdimen, 0, &r[15*spdimen], 0, &r[17*spdimen], 0, dvv );
      pkn_AddMatrixMd ( 1, spdimen, 0, dvv, 0, &r[16*spdimen], -2.0, 0, dvv );
      pkn_MultMatrixNumd ( 1, spdimen, 0, dvv, (double)(degreev*(degreev-1)),
                           0, dvv );
    }
    else
      memset ( dvv, 0, spdimen*sizeof(double) );
  }
  else {
    memset ( dv, 0, spdimen*sizeof(double) );
    memset ( dvv, 0, spdimen*sizeof(double) );
  }
  return true;
} /*_mbs_BCHornerDer2Pd*/

