
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the FillG1Holef procedure, constructing biquintic */
/* Bezier patches filling a polygonal hole in a surface made of      */
/* bicubic Bezier patches                                            */

#ifndef CONST_  /* a dirty trick to suppress many compiler warning messages */
#define CONST_ const
#endif

#ifndef G1HOLEF_H
#define G1HOLEF_H

#ifdef __cplusplus   
extern "C" {
#endif

/* Output routines may be hooked here in order to snoop into       */
/* the computation; in fact, the G1OutCentralPointf is allowed to  */
/* modify the central point. If the routine pointers are NULL then */
/* no output is done during the computations.                      */

extern void (*G1OutCentralPointf)( point3f *p );
extern void (*G1OutAuxCurvesf)( int ncurves, int degree, CONST_ point3f *accp,
                                float t );
extern void (*G1OutStarCurvesf)( int ncurves, int degree, CONST_ point3f *sccp );
extern void (*G1OutAuxPatchesf)( int npatches, int degu, int degv,
                                 CONST_ point3f *apcp );

boolean FillG1Holef ( int hole_k, point3f*(*GetBezp)(int i, int j),
                      float beta1, float beta2,
                      point3f *hpcp );

#ifdef __cplusplus
}
#endif

#endif

