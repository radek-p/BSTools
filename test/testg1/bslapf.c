
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "bslapf.h"

float ComputeBezLaplacianf ( int n, int m, const point3f *cp, float u, float v )
{
  vector3f p, pu, pv, puu, puv, pvv;
  float    fx, fy, fxx, fxy, fyy;
  vector2f gx, gy, gxx, gxy, gyy;
  float    A21[6], A22[9], trd[5], lap;

  mbs_BCHornerDer2P3f ( n, m, cp, u, v, &p, &pu, &pv, &puu, &puv, &pvv );

  pkn_f2iDerivatives2f ( pu.x, pu.y, pv.x, pv.y, puu.x, puu.y, puv.x, puv.y,
                         pvv.x, pvv.y, (float*)&gx, (float*)&gy,
                         (float*)&gxx, (float*)&gxy, (float*)&gyy );
  pkn_Setup2DerA21Matrixf ( gxx.x, gxx.y, gxy.x, gxy.y, gyy.x, gyy.y, A21 );
  pkn_Setup2DerA22Matrixf ( gx.x, gx.y, gy.x, gy.y, A22 );
  trd[0]  = A21[0]+A21[4];   trd[1]  = A21[1]+A21[5];
  trd[2]  = A22[0]+A22[6];   trd[3]  = A22[1]+A22[7];   trd[4] = A22[2]+A22[8];
  lap = trd[0]*pu.z + trd[1]*pv.z + trd[2]*puu.z + trd[3]*puv.z + trd[4]*pvv.z;

  pkn_Comp2iDerivatives2f ( pu.x, pu.y, pv.x, pv.y, puu.x, puu.y,
                            puv.x, puv.y, pvv.x, pvv.y,
                            1, &pu.z, &pv.z, &puu.z, &puv.z, &pvv.z,
                            &fx, &fy, &fxx, &fxy, &fyy );
  return fxx+fyy;
} /*ComputeBezLaplacianf*/

void ComputeBezLapGradf ( int n, int m, const point3f *cp, float u, float v,
                          vector2f *lapgr )
{
  vector3f p, pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv;
  float    fx, fy, fxx, fxy, fyy, fxxx, fxxy, fxyy, fyyy;

  mbs_BCHornerDer3P3f ( n, m, cp, u, v, &p, &pu, &pv, &puu, &puv, &pvv,
                        &puuu, &puuv, &puvv, &pvvv );
  pkn_Comp2iDerivatives3f ( pu.x, pu.y, pv.x, pv.y, puu.x, puu.y,
                    puv.x, puv.y, pvv.x, pvv.y, puuu.x, puuu.y,
                    puuv.x, puuv.y, puvv.x, puvv.y, pvvv.x, pvvv.y,
                    1, &pu.z, &pv.z, &puu.z, &puv.z, &pvv.z,
                    &puuu.z, &puuv.z, &puvv.z, &pvvv.z,
                    &fx, &fy, &fxx, &fxy, &fyy, &fxxx, &fxxy, &fxyy, &fyyy );
  SetVector2f ( lapgr, fxxx+fxyy, fxxy+fyyy );
} /*ComputeBezLapGradf*/

