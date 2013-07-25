
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the FillG1Holed procedure, constructing biquintic */
/* Bezier patches filling a polygonal hole in a surface made of      */
/* bicubic Bezier patches                                            */

#ifndef CONST_  /* a dirty trick to suppress many compiler warning messages */
#define CONST_ const
#endif

#ifndef G1HOLED_H
#define G1HOLED_H

#ifdef __cplusplus   
extern "C" {
#endif

/* Output routines may be hooked here in order to snoop into       */
/* the computation; in fact, the G1OutCentralPointd is allowed to  */
/* modify the central point. If the routine pointers are NULL then */
/* no output is done during the computations.                      */

extern void (*G1OutCentralPointd)( point3d *p );
extern void (*G1OutAuxCurvesd)( int ncurves, int degree, CONST_ point3d *accp,
                                double t );
extern void (*G1OutStarCurvesd)( int ncurves, int degree, CONST_ point3d *sccp );
extern void (*G1OutAuxPatchesd)( int npatches, int degu, int degv,
                                 CONST_ point3d *apcp );

boolean FillG1Holed ( int hole_k, point3d*(*GetBezp)(int i, int j),
                      double beta1, double beta2,
                      point3d *hpcp );

#ifdef __cplusplus
}
#endif

#endif

