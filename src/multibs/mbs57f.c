
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

boolean mbs_deBoorDer2Pf ( int degreeu, int lastknotu, const float *knotsu,
                           int degreev, int lastknotv, const float *knotsv,
                           int spdimen, int pitch, const float *ctlpoints,
                           float u, float v,
                           float *ppoint, float *uder, float *vder,   
                           float *uuder, float *uvder, float *vvder )
{
  void  *sp;
  int   k, l;
  float *q;

  sp = pkv_GetScratchMemTop ();
  if ( !(q = pkv_GetScratchMemf (3*(degreeu+1)*spdimen)) )
    goto failure;

  k = mbs_FindKnotIntervalf ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervalf ( degreev, lastknotv, knotsv, v, NULL );
  mbs_multideBoorDer2f ( degreev, 2*degreev+1,
        &knotsv[l-degreev], degreeu+1, spdimen, pitch,
        &ctlpoints[(k-degreeu)*pitch+(l-degreev)*spdimen],
        v, q, &q[(degreeu+1)*spdimen], &q[2*(degreeu+1)*spdimen] );
  mbs_multideBoorDer2f ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                         0, q, u, ppoint, uder, uuder );
  mbs_multideBoorDerf ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                        0, &q[(degreeu+1)*spdimen], u, vder, uvder );
  mbs_multideBoorf ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                     0, &q[2*(degreeu+1)*spdimen], u, vvder );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_deBoorDer2Pf*/

boolean mbs_deBoorDer3Pf ( int degreeu, int lastknotu, const float *knotsu,
                           int degreev, int lastknotv, const float *knotsv,
                           int spdimen, int pitch, const float *ctlpoints,
                           float u, float v,
                           float *ppoint, float *uder, float *vder,    
                           float *uuder, float *uvder, float *vvder,   
                           float *uuuder, float *uuvder, float *uvvder,
                           float *vvvder )
{
  void  *sp;
  int   k, l;
  float *q;

  sp = pkv_GetScratchMemTop ();
  if ( !(q = pkv_GetScratchMemf (4*(degreeu+1)*spdimen) ) )
    goto failure;

  k = mbs_FindKnotIntervalf ( degreeu, lastknotu, knotsu, u, NULL );
  l = mbs_FindKnotIntervalf ( degreev, lastknotv, knotsv, v, NULL );
  mbs_multideBoorDer3f ( degreev, 2*degreev+1,
        &knotsv[l-degreev], degreeu+1, spdimen, pitch,
        &ctlpoints[(k-degreeu)*pitch+(l-degreev)*spdimen],
        v, q, &q[(degreeu+1)*spdimen], &q[2*(degreeu+1)*spdimen],
        &q[3*(degreeu+1)*spdimen] );
  mbs_multideBoorDer3f ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                         0, q, u, ppoint, uder, uuder, uuuder );
  mbs_multideBoorDer2f ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                         0, &q[(degreeu+1)*spdimen], u, vder, uvder, uuvder );
  mbs_multideBoorDerf ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                        0, &q[2*(degreeu+1)*spdimen], u, vvder, uvvder );
  mbs_multideBoorf ( degreeu, 2*degreeu+1, &knotsu[k-degreeu], 1, spdimen,
                     0, &q[3*(degreeu+1)*spdimen], u, vvvder );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mbs_deBoorDer3Pf*/

