
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

#include "bslapd.h"

double ComputeBezLaplaciand ( int n, int m, const point3d *cp, double u, double v )
{
  vector3d p, pu, pv, puu, puv, pvv;
  double   fx, fy, fxx, fxy, fyy;

  mbs_BCHornerDer2P3d ( n, m, cp, u, v, &p, &pu, &pv, &puu, &puv, &pvv );
  pkn_Comp2iDerivatives2d ( pu.x, pu.y, pv.x, pv.y, puu.x, puu.y,
                            puv.x, puv.y, pvv.x, pvv.y,
                            1, &pu.z, &pv.z, &puu.z, &puv.z, &pvv.z,
                            &fx, &fy, &fxx, &fxy, &fyy );
  return fxx+fyy;
} /*ComputeBezLaplaciand*/

void ComputeBezLapGradd ( int n, int m, const point3d *cp, double u, double v,
                          vector2d *lapgr )
{
  vector3d p, pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv;
  double   fx, fy, fxx, fxy, fyy, fxxx, fxxy, fxyy, fyyy;

  mbs_BCHornerDer3P3d ( n, m, cp, u, v, &p, &pu, &pv, &puu, &puv, &pvv,
                        &puuu, &puuv, &puvv, &pvvv );
  pkn_Comp2iDerivatives3d ( pu.x, pu.y, pv.x, pv.y, puu.x, puu.y,
                    puv.x, puv.y, pvv.x, pvv.y, puuu.x, puuu.y,
                    puuv.x, puuv.y, puvv.x, puvv.y, pvvv.x, pvvv.y,
                    1, &pu.z, &pv.z, &puu.z, &puv.z, &pvv.z,
                    &puuu.z, &puuv.z, &puvv.z, &pvvv.z,
                    &fx, &fy, &fxx, &fxy, &fyy, &fxxx, &fxxy, &fxyy, &fyyy );
  SetVector2d ( lapgr, fxxx+fxyy, fxxy+fyyy );
} /*ComputeBezLapGradd*/

