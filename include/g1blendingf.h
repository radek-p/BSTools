
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Mateusz Markowski                                */
/* and modified by Przemyslaw Kiciak                                         */

#ifndef G1BLENDINGF_H
#define G1BLENDINGF_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif
#ifndef PKNUM_H
#include "pknum.h"
#endif
#ifndef PKGEOM_H
#include "pkgeom.h"
#endif
#ifndef MULTIBS_H
#include "multibs.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif


boolean g1bl_SetupBiharmAMatrixf ( int lastknotu, int lastknotv,
                               int *n, int **prof, float **Amat, float ***arow );
boolean g1bl_SetupBiharmRHSf ( int lastknotu, int lastknotv,
                               int spdimen, int pitch, const float *cpoints,
                               float *rhs );

boolean g1bl_SetupClosedBiharmAMatrixf ( int lastknotu, int lastknotv,
                               int *n, int **prof, float **Amat, float ***arow );
boolean g1bl_SetupClosedBiharmRHSf ( int lastknotu, int lastknotv,
                               int spdimen, int pitch, const float *cpoints,
                               float *rhs );

#ifdef __cplusplus
}
#endif

#endif /*G1BLENDINGF_H*/

