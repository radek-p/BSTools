
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef G2BLENDINGD_H
#define G2BLENDINGD_H

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
boolean g2bl_SetupTriharmAMatrixd ( int lastknotu, int lastknotv,
                           int *n, int **prof, double **Amat, double ***arow );
boolean g2bl_SetupTriharmRHSd ( int lastknotu, int lastknotv,
                           int spdimen, int pitch, const double *cpoints,
                           double *rhs );

boolean g2bl_SetupClosedTriharmAMatrixd ( int lastknotu, int lastknotv,
                           int *n, int **prof, double **Amat, double ***arow );
boolean g2bl_SetupClosedTriharmRHSd ( int lastknotu, int lastknotv,
                           int spdimen, int pitch, const double *cpoints,
                           double *rhs );

/* procedures implementing the construction of minimal blending patches */
/* of a nonlinear functional measuring surface shape badness */

int g2bl_NiSize ( int nkn );
int g2bl_NijSize ( int nkn );
int g2bl_MijSize ( int nkn );

int _g2bl_SetupHessian1Profile ( int lastknotu, int lastknotv, int *prof );

double g2bl_UFuncd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv, int pitch, point3d *cp,
                     char *dirty,
                     double tC, double *ftab );
void g2bl_UFuncGradd ( int nkn, const double *qcoeff, double *Nitab,
                       int lastknotu, int lastknotv,
                       int pitch, point3d *cp, char *dirty,
                       double tC, double *ftab, double *gtab,
                       double *func, double *grad );
void g2bl_UFuncGradHessiand ( int nkn, const double *qcoeff, double *Nitab,
                              double *Nijtab, double *Mijtab,
                              int lastknotu, int lastknotv,
                              int pitch, point3d *cp, char *dirty,
                              double tC, double *ftab, double *gtab, double *htab,
                              double *func, double *grad,
                              int hsize, const int *prof, double **hrows );
void g2bl_ClosedUFuncGradd ( int nkn, const double *qcoeff, double *Nitab,
                       int lastknotu, int lastknotv,
                       int pitch, point3d *cp, char *dirty,
                       double tC, double *ftab, double *gtab,
                       double *func, double *grad );
void g2bl_ClosedUFuncGradHessiand ( int nkn, const double *qcoeff, double *Nitab,
                              double *Nijtab, double *Mijtab,
                              int lastknotu, int lastknotv,
                              int pitch, point3d *cp, char *dirty,
                              double tC, double *ftab, double *gtab, double *htab,
                              double *func, double *grad,
                              int hsize, const int *prof, double **hrows );

double g2bl_SurfNetDiameterSqd ( int lastknotu, int lastknotv,
                                 int pitch, const point3d *cp );
double g2bl_ClosedSurfNetDiameterSqd ( int lastknotu, int lastknotv,
                                       int pitch, const point3d *cp );


boolean g2bl_InitBlSurfaceOptLMTd ( int lastknotu, int lastknotv, int pitch,
                                    point3d *cp,
                                    double C, double dO, double dM,
                                    int nkn1, int nkn2,
                                    void **data );
boolean g2bl_IterBlSurfaceOptLMTd ( void *data, boolean *finished );
void g2bl_OptLMTDeallocated ( void **data );
boolean g2bl_FindBlSurfaceLMTd ( int lastknotu, int lastknotv, int pitch,
                                 point3d *cp,
                                 double C, double dO, double dM,
                                 int maxit, int nkn1, int nkn2 );


boolean g2bl_InitBlSurfaceConstrOptLMTd ( int lastknotu, int lastknotv, int pitch,
                                          point3d *cp,
                                          int nconstr, double *constrmat,
                                          double *constrrhs,
                                          double C, double dO, double dM,
                                          int nkn1, int nkn2,
                                          void **data );
boolean g2bl_IterBlSurfaceConstrOptLMTd ( void *data, boolean *finished );
void g2bl_ConstrOptLMTDeallocated ( void **data );
boolean g2bl_FindBlSurfaceConstrLMTd ( int lastknotu, int lastknotv, int pitch,
                                       point3d *cp,
                                       int nconstr, double *constrmat,
                                       double *constrrhs,
                                       double C, double dO, double dM,
                                       int maxit, int nkn1, int nkn2 );


boolean g2bl_ClosedInitBlSurfaceOptLMTd ( int lastknotu, int lastknotv, int pitch,
                                          point3d *cp,
                                          double C, double dO, double dM,
                                          int nkn1, int nkn2,
                                          void **data );
boolean g2bl_ClosedIterBlSurfaceOptLMTd ( void *data, boolean *finished );
void g2bl_ClosedOptLMTDeallocated ( void **data );
boolean g2bl_ClosedFindBlSurfaceLMTd ( int lastknotu, int lastknotv, int pitch,
                                       point3d *cp,
                                       double C, double dO, double dM,
                                       int maxit, int nkn1, int nkn2 );


boolean g2bl_ClosedInitBlSurfaceConstrOptLMTd (
                   int lastknotu, int lastknotv, int pitch, point3d *cp,
                   int nconstr, double *constrmat, double *constrrhs,
                   double C, double dO, double dM, int nkn1, int nkn2,
                   void **data );
boolean g2bl_ClosedIterBlSurfaceConstrOptLMTd ( void *data, boolean *finished );
void g2bl_ClosedConstrOptLMTDeallocated ( void **data );
boolean g2bl_ClosedFindBlSurfaceConstrLMTd (
                   int lastknotu, int lastknotv, int pitch, point3d *cp,
                   int nconstr, double *constrmat, double *constrrhs,
                   double C, double dO, double dM,
                   int maxit, int nkn1, int nkn2 );


boolean g2bl_SetupULConstraintsd ( int lastknotu, int lastknotv, int spdimen,
                                   int ppitch, double *cp,
                                   int nucurv, double *ucknots,
                                   int cpitch, double *uccp,
                                   int *nconstr, double *cmat, double *crhs );

boolean g2bl_SetupUNLConstraintsd ( int lastknotu, int lastknotv,
                                    int ppitch, point3d *cp,
                                    int nucurv, double *ucknots,
                                    int cpitch, point3d *uccp,
                                    int *nconstr, double *cmat, double *crhs );


boolean g2bl_SetupClosedULConstraintsd ( int lastknotu, int lastknotv, int spdimen,
                                         int ppitch, double *cp,
                                         int nucurv, double *ucknots,
                                         int cpitch, double *uccp,
                                         int *nconstr, double *cmat, double *crhs );

boolean g2bl_SetupClosedUNLConstraintsd ( int lastknotu, int lastknotv,
                                          int ppitch, point3d *cp,
                                          int nucurv, double *ucknots,
                                          int cpitch, point3d *uccp,
                                          int *nconstr, double *cmat, double *crhs );


boolean g2bl_FuncTSQFd ( int nkn,
                         int lastknotu, int lastknotv, int pitch, point3d *cp,
                         double tC,
                         double *fT, double *fS, double *fQ, double *fF );

#ifdef __cplusplus
}
#endif

#endif /*G2BLENDINGD_H*/

