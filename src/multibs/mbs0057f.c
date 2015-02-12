
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

static int FindMaxInt ( int a, int b, int c, int d, int e )
{
  if ( b > a ) a = b;  if ( c > a ) a = c;  if ( d > a ) a = d;
  return e > a ? e : a;
} /*FindMaxInt*/

boolean _mbs_BezC1CoonsToBezf ( int spdimen,
                                int degc00, const float *c00,
                                int degc01, const float *c01,
                                int degc10, const float *c10,
                                int degc11, const float *c11,
                                int degd00, const float *d00,
                                int degd01, const float *d01,
                                int degd10, const float *d10,
                                int degd11, const float *d11,
                                int *n, int *m, float *p,
                                float *workspace )
{
  int   degu, degv, d, du, dv;
  int   pitch;
  float *pc, *bc;
  float *p1, *p2, *p3; 
  int   i;

  *n = degu = FindMaxInt ( degc00, degc01, degc10, degc11, 3 );
  *m = degv = FindMaxInt ( degd00, degd01, degd10, degd11, 3 );
  pc = &workspace[2*spdimen];
  p1 = &pc[16*spdimen];
  p2 = &p1[(degu+1)*(degv+1)*spdimen];
  p3 = &p2[(degu+1)*(degv+1)*spdimen];
  bc = &p3[(degu+1)*(degv+1)*spdimen];

        /* construct the patch p1 */
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc00, c00, degu-degc00,
                              0, &d, &bc[0] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc01, c01, degu-degc01,
                              0, &d, &bc[(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc10, c10, degu-degc10,
                              0, &d, &bc[2*(degu+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degc11, c11, degu-degc11,
                              0, &d, &bc[3*(degu+1)*spdimen] ) )
    return false;
  pitch = (degu+1)*spdimen;
  if ( !mbs_multiInterp2knHermiteBezf ( 1, pitch, 3, 2, pitch, bc,
                       2, pitch, &bc[2*(degu+1)*spdimen], pitch, p1 ) )
    return false;
  pkv_TransposeMatrixc ( 4, degu+1, spdimen*sizeof(float),
                         (degu+1)*spdimen*sizeof(float), (char*)p1,
                         4*spdimen*sizeof(float), (char*)p );
  if ( !mbs_BCDegElevPf ( spdimen, degu, 3, p, 0, degv-3, &du, &dv, p ) )
    return false;

        /* construct the patch p2 */
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd00, d00, degv-degd00,
                              0, &d, &bc[0] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd01, d01, degv-degd01,
                              0, &d, &bc[(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd10, d10, degv-degd10,
                              0, &d, &bc[2*(degv+1)*spdimen] ) )
    return false;
  if ( !mbs_multiBCDegElevf ( 1, spdimen, 0, degd11, d11, degv-degd11,
                              0, &d, &bc[3*(degv+1)*spdimen] ) )
    return false;
  pitch = (degv+1)*spdimen;
  if ( !mbs_multiInterp2knHermiteBezf ( 1, pitch, 3, 2, pitch, bc,
                       2, pitch, &bc[2*(degv+1)*spdimen], pitch, p2 ) )
    return false;
  if ( !mbs_BCDegElevPf ( spdimen, 3, degv, p2, degu-3, 0, &du, &dv, p2 ) )
    return false;
  pkn_AddMatrixf ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p2, 0, p );

        /* construct the patch p3 */
  if ( !_mbs_BezC1CoonsFindCornersf ( spdimen, degc00, c00, degc01, c01,
                       degc10, c10, degc11, c11, pc, workspace ) )
    return false;

  pkv_Selectf ( 2, 4*spdimen, 4*2*spdimen, 4*spdimen, pc, p3 );
  pkv_Selectf ( 2, 4*spdimen, 4*2*spdimen, 4*spdimen, &pc[4*spdimen], &p3[4*2*spdimen] );
  for ( i = 0; i < 4; i++ ) {
    pkv_Selectf ( 2, spdimen, 2*spdimen, spdimen, &p3[4*i*spdimen], &bc[4*i*spdimen] );
    pkv_Selectf ( 2, spdimen, 2*spdimen, spdimen, &p3[(4*i+1)*spdimen], &bc[(4*i+2)*spdimen] );
  }
  if ( !mbs_multiInterp2knHermiteBezf ( 1, 4*spdimen, 3, 2, 0, bc,
                                        2, 0, &bc[4*2*spdimen], 0, p3 ) )
    return false;
  if ( !mbs_multiInterp2knHermiteBezf ( 4, spdimen, 3, 2, 4*spdimen, p3,
                                  2, 4*spdimen, &p3[2*spdimen], 4*spdimen, p2 ) )
    return false;
  if ( !mbs_BCDegElevPf ( spdimen, 3, 3, p2, degu-3, degv-3, &du, &dv, p3 ) )
    return false;
  pkn_SubtractMatrixf ( 1, spdimen*(degu+1)*(degv+1), 0, p, 0, p3, 0, p );
  return true;
} /*_mbs_BezC1CoonsToBezf*/

