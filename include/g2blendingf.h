
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef G2BLENDINGF_H
#define G2BLENDINGF_H

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


/* procedures implementing the construction of triharmonic patches */
boolean g2bl_SetupTriharmAMatrixf ( int lastknotu, int lastknotv,
                           int *n, int **prof, float **Amat, float ***arow );
boolean g2bl_SetupTriharmRHSf ( int lastknotu, int lastknotv,
                           int spdimen, int pitch, const float *cpoints,
                           float *rhs );

boolean g2bl_SetupClosedTriharmAMatrixf ( int lastknotu, int lastknotv,
                           int *n, int **prof, float **Amat, float ***arow );
boolean g2bl_SetupClosedTriharmRHSf ( int lastknotu, int lastknotv,
                           int spdimen, int pitch, const float *cpoints,
                           float *rhs );

/* procedures implementing the construction of minimal blending patches */
/* of a nonlinear functional measuring surface shape badness have been  */
/* removed; floating point overflow and underflow make the computations */
/* in single precision impossible */

#ifdef __cplusplus
}
#endif

#endif /*G2BLENDINGF_H*/

