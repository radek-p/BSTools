
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

#ifndef G1BLENDINGD_H
#define G1BLENDINGD_H

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

int g1bl_NiSize ( int nkn );
int g1bl_NijSize ( int nkn );
int g1bl_MijSize ( int nkn ); 

boolean g1bl_SetupBiharmAMatrixd ( int lastknotu, int lastknotv,
                               int *n, int **prof, double **Amat, double ***arow );
boolean g1bl_SetupBiharmRHSd ( int lastknotu, int lastknotv,
                               int spdimen, int pitch, const double *cpoints,
                               double *rhs );

boolean g1bl_SetupClosedBiharmAMatrixd ( int lastknotu, int lastknotv,
                               int *n, int **prof, double **Amat, double ***arow );
boolean g1bl_SetupClosedBiharmRHSd ( int lastknotu, int lastknotv,
                               int spdimen, int pitch, const double *cpoints,
                               double *rhs );

void g1bl_TabNid ( int nkn, double *bf, double *dbf, double *ddbf,
                   double *Nitab );
void g1bl_TabNijd ( int nkn, double *bf, double *dbf, double *ddbf,
                    double *Nijtab );
double g1bl_UFuncd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv, int pitch, point3d *cp,
                     char *dirty,
                     double tC, double *ftab );
		     

double g1bl_QFuncd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv, int pitch, point3d *cp,
                     char *dirty,
                     double tC, double *ftab );
		     
double g1bl_biharmFuncd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv, int pitch, point3d *cp,
                     char *dirty,
                     double tC, double *ftab );
		     
void g1bl_UFuncGradd ( int nkn, const double *qcoeff, double *Nitab,
                       int lastknotu, int lastknotv,
                       int pitch, point3d *cp, char *dirty,
                       double tC, double *ftab, double *gtab,
                       double *func, double *grad );
void g1bl_UFuncGradHessiand ( int nkn, const double *qcoeff, double *Nitab,
                              double *Nijtab, double *Mijtab,
                              int lastknotu, int lastknotv,
                              int pitch, point3d *cp, char *dirty,
                              double tC, double *ftab, double *gtab, double *htab,
                              double *func, double *grad,
                              int hsize, const int *prof, double **hrows );
double g1bl_SurfNetDiameterSqd ( int lastknotu, int lastknotv,
                                 int pitch, const point3d *cp );

boolean g1bl_InitBlSurfaceOptLMTd ( int lastknotu, int lastknotv, int pitch,
                                    point3d *cp,
                                    double C, double dO, double dM,
                                    int nkn1, int nkn2,
                                    void **data );
boolean g1bl_IterBlSurfaceOptLMTd ( void *data, boolean *finished );
void g1bl_OptLMTDeallocated ( void **data );
boolean g1bl_FindBlSurfaceLMTd ( int lastknotu, int lastknotv, int pitch,
                                 point3d *cp,
                                 double C, double dO, double dM,
                                 int maxit, int nkn1, int nkn2 );


boolean g1bl_InitBlSurfaceConstrOptLMTd ( int lastknotu, int lastknotv, int pitch,
                                          point3d *cp,
                                          int nconstr, double *constrmat,
                                          double *constrrhs,
                                          double C, double dO, double dM,
                                          int nkn1, int nkn2,
                                          void **data );
boolean g1bl_IterBlSurfaceConstrOptLMTd ( void *data, boolean *finished );
void g1bl_ConstrOptLMTDeallocated ( void **data );
boolean g1bl_FindBlSurfaceConstrLMTd ( int lastknotu, int lastknotv, int pitch,
                                       point3d *cp,
                                       int nconstr, double *constrmat,
                                       double *constrrhs,
                                       double C, double dO, double dM,
                                       int maxit, int nkn1, int nkn2 );


boolean g1bl_ClosedInitBlSurfaceOptLMTd ( int lastknotu, int lastknotv, int pitch,
                                          point3d *cp,
                                          double C, double dO, double dM,
                                          int nkn1, int nkn2,
                                          void **data );
boolean g1bl_ClosedIterBlSurfaceOptLMTd ( void *data, boolean *finished );
void g1bl_ClosedOptLMTDeallocated ( void **data );
boolean g1bl_ClosedFindBlSurfaceLMTd ( int lastknotu, int lastknotv, int pitch,
                                       point3d *cp,
                                       double C, double dO, double dM,
                                       int maxit, int nkn1, int nkn2 );


boolean g1bl_ClosedInitBlSurfaceConstrOptLMTd (
                   int lastknotu, int lastknotv, int pitch, point3d *cp,
                   int nconstr, double *constrmat, double *constrrhs,
                   double C, double dO, double dM, int nkn1, int nkn2,
                   void **data );
boolean g1bl_ClosedIterBlSurfaceConstrOptLMTd ( void *data, boolean *finished );
void g1bl_ClosedConstrOptLMTDeallocated ( void **data );
boolean g1bl_ClosedFindBlSurfaceConstrLMTd (
                   int lastknotu, int lastknotv, int pitch, point3d *cp,
                   int nconstr, double *constrmat, double *constrrhs,
                   double C, double dO, double dM,
                   int maxit, int nkn1, int nkn2 );

boolean g1bl_SetupULConstraintsd ( int lastknotu, int lastknotv, int spdimen,
                                   int ppitch, double *cp,
                                   int nucurv, double *ucknots,
                                   int cpitch, double *uccp,
                                   int *nconstr, double *cmat, double *crhs );

boolean g1bl_SetupUNLConstraintsd ( int lastknotu, int lastknotv,
                                    int ppitch, point3d *cp,
                                    int nucurv, double *ucknots,
                                    int cpitch, point3d *uccp,
                                    int *nconstr, double *cmat, double *crhs );


boolean g1bl_SetupClosedULConstraintsd ( int lastknotu, int lastknotv, int spdimen,
                                         int ppitch, double *cp,
                                         int nucurv, double *ucknots,
                                         int cpitch, double *uccp,
                                         int *nconstr, double *cmat, double *crhs );

boolean g1bl_SetupClosedUNLConstraintsd ( int lastknotu, int lastknotv,
                                          int ppitch, point3d *cp,
                                          int nucurv, double *ucknots,
                                          int cpitch, point3d *uccp,
                                          int *nconstr, double *cmat, double *crhs );


boolean g1bl_FuncTSQFd ( int nkn,
                         int lastknotu, int lastknotv, int pitch, point3d *cp,
                         double tC,
                         double *fT, double *fS, double *fQ, double *fF );

#ifdef __cplusplus
}
#endif

#endif /*G1BLENDINGD_H*/

