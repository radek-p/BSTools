
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2008, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"

boolean mbs_deBoorDer2Pd ( int degreeu, int lastknotu, const double *knotsu,
                           int degreev, int lastknotv, const double *knotsv,
                           int spdimen, int pitch, const double *ctlpoints,
                           double u, double v,
                           double *ppoint, double *uder, double *vder,   
                           double *uuder, double *uvder, double *vvder )
{
  void   *sp;
  int    k, l;
  double *q;

  sp = pkv_GetScratchMemTop ();
  if ( !(q = pkv_GetScratchMemd (3*(degreeu+1)*spdimen)) )
    goto failure;

  k = mbs_FindKnotIntervald ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervald ( degreev, lastknotv, knotsv, v, NULL );
  mbs_multideBoorDer2d ( degreev, 2*degreev+1,
        &knotsv[l-degreev], degreeu+1, spdimen, pitch,
        &ctlpoints[(k-degreeu)*pitch+(l-degreev)*spdimen],
        v, q, &q[(degreeu+1)*spdimen], &q[2*(degreeu+1)*spdimen] );
  mbs_multideBoorDer2d ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                         0, q, u, ppoint, uder, uuder );
  mbs_multideBoorDerd ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                        0, &q[(degreeu+1)*spdimen], u, vder, uvder );
  mbs_multideBoord ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                     0, &q[2*(degreeu+1)*spdimen], u, vvder );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_deBoorDer2Pd*/

boolean mbs_deBoorDer3Pd ( int degreeu, int lastknotu, const double *knotsu,
                           int degreev, int lastknotv, const double *knotsv,
                           int spdimen, int pitch, const double *ctlpoints,
                           double u, double v,
                           double *ppoint, double *uder, double *vder,    
                           double *uuder, double *uvder, double *vvder,   
                           double *uuuder, double *uuvder, double *uvvder,
                           double *vvvder )
{
  void   *sp;
  int    k, l;
  double *q;

  sp = pkv_GetScratchMemTop ();
  if ( !(q = pkv_GetScratchMemd (4*(degreeu+1)*spdimen) ) )
    goto failure;

  k = mbs_FindKnotIntervald ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervald ( degreev, lastknotv, knotsv, v, NULL );
  mbs_multideBoorDer3d ( degreev, 2*degreev+1,
        &knotsv[l-degreev], degreeu+1, spdimen, pitch,
        &ctlpoints[(k-degreeu)*pitch+(l-degreev)*spdimen],
        v, q, &q[(degreeu+1)*spdimen], &q[2*(degreeu+1)*spdimen],
        &q[3*(degreeu+1)*spdimen] );
  mbs_multideBoorDer3d ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                         0, q, u, ppoint, uder, uuder, uuuder );
  mbs_multideBoorDer2d ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                         0, &q[(degreeu+1)*spdimen], u, vder, uvder, uuvder );
  mbs_multideBoorDerd ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                        0, &q[2*(degreeu+1)*spdimen], u, vvder, uvvder );
  mbs_multideBoord ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                     0, &q[3*(degreeu+1)*spdimen], u, vvvder );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_deBoorDer3Pd*/

