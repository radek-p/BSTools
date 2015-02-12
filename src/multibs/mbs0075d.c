
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

static int FindMaxInt ( int a, int b, int c, int d, int e, int f, int g )
{
  if ( b > a ) a = b;  if ( c > a ) a = c;  if ( d > a ) a = d;
  if ( e > a ) a = e;  if ( f > a ) a = f;
  return g > a ? g : a;
} /*FindMaxInt*/

boolean _mbs_BezC2CoonsToBezd ( int spdimen,
                                int degc00, const double *c00,
                                int degc01, const double *c01,
                                int degc02, const double *c02,
                                int degc10, const double *c10,
                                int degc11, const double *c11,
                                int degc12, const double *c12,
                                int degd00, const double *d00,
                                int degd01, const double *d01,
                                int degd02, const double *d02,
                                int degd10, const double *d10,
                                int degd11, const double *d11,
                                int degd12, const double *d12,
                                int *n, int *m, double *p,
                                double *workspace )
{
  int    degu, degv, d, du, dv;
  int    pitch;
  double *pc, *bc, *aux;
  double *p1, *p2, *p3; 
  int    i;

  *n = degu = FindMaxInt ( degc00, degc01, degc02, degc10, degc11, degc12, 5 );
  *m = degv = FindMaxInt ( degd00, degd01, degd02, degd10, degd11, degd12, 5 );
  pc = workspace;
  p1 = &pc[36*spdimen];
  p2 = &p1[(degu+1)*(degv+1)*spdimen];
  p3 = &p2[(degu+1)*(degv+1)*spdimen];
  bc = &p3[(degu+1)*(degv+1)*spdimen];
  aux = &bc[6*(max(degu,degv)+1)*spdimen];
        /* construct the patch p1 */
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degc00, c00, degu-degc00,
                              0, &d, &bc[0] ) )
    return false;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degc01, c01, degu-degc01,
                              0, &d, &bc[(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degc02, c02, degu-degc02,
                              0, &d, &bc[2*(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degc10, c10, degu-degc10,
                              0, &d, &bc[3*(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degc11, c11, degu-degc11,
                              0, &d, &bc[4*(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degc12, c12, degu-degc12,
                              0, &d, &bc[5*(degu+1)*spdimen] ) )
    return false;
  pitch = (degu+1)*spdimen;
  if ( !mbs_multiInterp2knHermiteBezd ( 1, pitch, 5, 3, pitch, bc,
           3, pitch, &bc[3*(degu+1)*spdimen], pitch, p1 ) )
    return false;
  pkv_TransposeMatrixc ( 6, degu+1, spdimen*sizeof(double),
                         (degu+1)*spdimen*sizeof(double), (char*)p1,
                         6*spdimen*sizeof(double), (char*)p );
  if ( !mbs_BCDegElevPd ( spdimen, degu, 5, p, 0, degv-5, &du, &dv, p ) )
    return false;

        /* construct the patch p2 */
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degd00, d00, degv-degd00,
                              0, &d, &bc[0] ) )
    return false;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degd01, d01, degv-degd01,
                              0, &d, &bc[(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degd02, d02, degv-degd02,
                              0, &d, &bc[2*(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degd10, d10, degv-degd10,
                              0, &d, &bc[3*(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degd11, d11, degv-degd11,
                              0, &d, &bc[4*(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevd ( 1, spdimen, 0, degd12, d12, degv-degd12,
                              0, &d, &bc[5*(degv+1)*spdimen] ) )
    return false;
  pitch = (degv+1)*spdimen;
  if ( !mbs_multiInterp2knHermiteBezd ( 1, pitch, 5, 3, pitch, bc,
           3, pitch, &bc[3*(degv+1)*spdimen], pitch, p2 ) )
    return false;
  if ( !mbs_BCDegElevPd ( spdimen, 5, degv, p2, degu-5, 0, &du, &dv, p2 ) )
    return false;
  pkn_AddMatrixd ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p2, 0, p );

        /* construct the patch p3 */
  if ( !_mbs_BezC2CoonsFindCornersd ( spdimen,
                       degc00, c00, degc01, c01, degc02, c02,
                       degc10, c10, degc11, c11, degc12, c12, pc, aux ) )
    return false;

  pkv_Selectd ( 3, 6*spdimen, 6*2*spdimen, 6*spdimen, pc, p3 );
  pkv_Selectd ( 3, 6*spdimen, 6*2*spdimen, 6*spdimen, &pc[6*spdimen], &p3[6*3*spdimen] );
  for ( i = 0; i < 6; i++ ) {
    pkv_Selectd ( 3, spdimen, 2*spdimen, spdimen, &p3[6*i*spdimen], &bc[6*i*spdimen] );
    pkv_Selectd ( 3, spdimen, 2*spdimen, spdimen, &p3[(6*i+1)*spdimen], &bc[(6*i+3)*spdimen] );
  }
  if ( !mbs_multiInterp2knHermiteBezd ( 1, 6*spdimen, 5, 3, 0, bc,
                                        3, 0, &bc[6*3*spdimen], 0, p3 ) )
    return false;
  if ( !mbs_multiInterp2knHermiteBezd ( 6, spdimen, 5, 3, 6*spdimen, p3,
                                  3, 6*spdimen, &p3[3*spdimen], 6*spdimen, p2 ) )
    return false;
  if ( !mbs_BCDegElevPd ( spdimen, 5, 5, p2, degu-5, degv-5, &du, &dv, p3 ) )
    return false;
  pkn_SubtractMatrixd ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p3, 0, p );
  return true;
} /*_mbs_BezC2CoonsToBezd*/

