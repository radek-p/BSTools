
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

  mbs_BCHornerDer2P3f ( n, m, cp, u, v, &p, &pu, &pv, &puu, &puv, &pvv );
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

