
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2009                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libmultibs library of C procedures -              */
/* processing B-spline and Bezier curves and surfaces                    */ 

#ifndef MULTIBS_H
#define MULTIBS_H

#ifndef PKNUM_H
#include "pknum.h"
#endif

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#ifdef __cplusplus   
extern "C" {
#endif


/* boundary conditions identifiers for cubic splines of interpolation */

#define BS3_BC_FIRST_DER   0
#define BS3_BC_FIRST_DER0  1
#define BS3_BC_SECOND_DER  2
#define BS3_BC_SECOND_DER0 3
#define BS3_BC_THIRD_DER   4
#define BS3_BC_THIRD_DER0  5
#define BS3_BC_BESSEL      6
#define BS3_BC_NOT_A_KNOT  7

void mbs_BezP3NormalDeg ( int degreeu, int degreev, int *ndegu, int *ndegv );
void mbs_BezP3RNormalDeg ( int degreeu, int degreev, int *ndegu, int *ndegv );

#include "multibsf.h"
#include "multibsd.h"

#ifdef __cplusplus
}
#endif

#endif /*MULTIBS_H*/

