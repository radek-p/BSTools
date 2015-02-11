
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2015                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libmultibs library of C procedures -              */
/* processing B-spline and Bezier curves and surfaces                    */ 

#ifndef CONST_  /* a dirty trick to suppress many compiler warning messages */
#define CONST_ const
#endif

#ifndef MULTIBSD_H
#define MULTIBSD_H

#ifndef PKNUM_H
#include "pknum.h"
#endif

#ifndef PKGEOM_H
#include "pkgeom.h"
#endif

#ifdef __cplusplus   
extern "C" {
#endif

#ifndef MULTIBS_H
/* boundary conditions identifiers for cubic splines of interpolation */

#define BS3_BC_FIRST_DER   0
#define BS3_BC_FIRST_DER0  1
#define BS3_BC_SECOND_DER  2
#define BS3_BC_SECOND_DER0 3
#define BS3_BC_THIRD_DER   4
#define BS3_BC_THIRD_DER0  5
#define BS3_BC_BESSEL      6
#define BS3_BC_NOT_A_KNOT  7
#endif


int mbs_KnotMultiplicityd ( int lastknot, const double *knots, double t );

int mbs_FindKnotIntervald ( int degree, int lastknot, const double *knots,
                            double t, int *mult );

double mbs_GrevilleAbscissad ( int degree, double *knots, int i );

void mbs_TransformAffKnotsd ( int degree, int lastknot, const double *inknots,
                              double a, double b, double *outknots );

void mbs_multiReverseBSCurved ( int degree, int lastknot, double *knots,
                                int ncurves, int spdimen, int pitch,
                                double *ctlpoints );
boolean mbs_ClosedKnotsCorrectd ( int degree, int lastknot, double *knots,
                                  double T, int K, double tol );

int mbs_SetKnotd ( int lastknot, double *knots,
                   int knotnum, int mult, double t );
int mbs_SetKnotClosedd ( int degree, int lastknot, double *knots, double T,
                         int knotnum, int mult, double t );
                                                  
void _mbs_multideBoorKerneld ( int degree, const double *knots,
                               int ncurves, int spdimen,
                               int pitch, const double *ctlpoints,
                               double t, int k, int r, int lj,
                               int dpitch, double *d );

int mbs_multideBoord ( int degree, int lastknot,
                       const double *knots,
                       int ncurves, int spdimen, int pitch,
                       const double *ctlpoints,
                       double t,
                       double *cpoints );

#define mbs_deBoorC1d(degree,lastknot,knots,coeff,t,value) \
  mbs_multideBoord ( degree, lastknot, knots, 1, 1, 0, coeff, t, value )
#define mbs_deBoorC2d(degree,lastknot,knots,ctlpoints,t,cpoint) \
  mbs_multideBoord ( degree, lastknot, knots, 1, 2, 0, (double*)ctlpoints, t, \
    (double*)cpoint )
#define mbs_deBoorC3d(degree,lastknot,knots,ctlpoints,t,cpoint) \
  mbs_multideBoord ( degree, lastknot, knots, 1, 3, 0, (double*)ctlpoints, t, \
    (double*)cpoint )
#define mbs_deBoorC4d(degree,lastknot,knots,ctlpoints,t,cpoint) \
  mbs_multideBoord ( degree, lastknot, knots, 1, 4, 0, (double*)ctlpoints, t, \
    (double*)cpoint )


void mbs_deBoorC2Rd ( int degree, int lastknot,
                      const double *knots, const point3d *ctlpoints,
                      double t, point2d *cpoint );

void mbs_deBoorC3Rd ( int degree, int lastknot,
                      const double *knots, const point4d *ctlpoints,
                      double t, point3d *cpoint );

boolean mbs_deBoorPd ( int degreeu, int lastknotu, const double *knotsu,
                       int degreev, int lastknotv, const double *knotsv,
                       int spdimen, int pitch,
                       const double *ctlpoints,
                       double u, double v, double *ppoint );

boolean mbs_deBoorP3d ( int degreeu, int lastknotu, const double *knotsu,
                        int degreev, int lastknotv, const double *knotsv,
                        int pitch,
                        const point3d *ctlpoints,
                        double u, double v, point3d *ppoint );

boolean mbs_deBoorP3Rd ( int degreeu, int lastknotu, const double *knotsu,
                         int degreev, int lastknotv, const double *knotsv,
                         int pitch,
                         const point4d *ctlpoints,
                         double u, double v, point3d *ppoint );

boolean mbs_deBoorP4d ( int degreeu, int lastknotu, const double *knotsu,
                        int degreev, int lastknotv, const double *knotsv,
                        int pitch,
                        const point4d *ctlpoints,
                        double u, double v, point4d *ppoint );


int mbs_multideBoorDerd ( int degree, int lastknot,
                          const double *knots,
                          int ncurves, int spdimen, int pitch,
                          const double *ctlpoints,
                          double t,
                          double *cpoints, double *dervect );

#define mbs_deBoorDerC1d(degree,lastknot,knots,ctlpoints,t,cpoint,cder) \
  mbs_multideBoorDerd ( degree, lastknot, knots, 1, 1, 0, (double*)ctlpoints, \
    t, cpoint, cder )
#define mbs_deBoorDerC2d(degree,lastknot,knots,ctlpoints,t,cpoint,cder) \
  mbs_multideBoorDerd ( degree, lastknot, knots, 1, 2, 0, (double*)ctlpoints, \
    t, (double*)cpoint, (double*)cder )
#define mbs_deBoorDerC3d(degree,lastknot,knots,ctlpoints,t,cpoint,cder) \
  mbs_multideBoorDerd ( degree, lastknot, knots, 1, 3, 0, (double*)ctlpoints, \
    t, (double*)cpoint, (double*)cder )
#define mbs_deBoorDerC4d(degree,lastknot,knots,ctlpoints,t,cpoint,cder) \
  mbs_multideBoorDerd ( degree, lastknot, knots, 1, 4, 0, (double*)ctlpoints, \
    t, (double*)cpoint, (double*)cder )


int mbs_multideBoorDer2d ( int degree, int lastknot, const double *knots,
                           int ncurves, int spdimen,
                           int pitch, const double *ctlpoints,
                           double t, double *p, double *d1, double *d2 );

#define mbs_deBoorDer2C1d(degree,lastknot,knots,coeff,t,p,d1,d2) \
  mbs_multideBoorDer2d(degree,lastknot,knots,1,1,0,coeff,t,p,d1,d2)
#define mbs_deBoorDer2C2d(degree,lastknot,knots,ctlpoints,t,p,d1,d2) \
  mbs_multideBoorDer2d(degree,lastknot,knots,1,2,0,(double*)ctlpoints,t, \
    (double*)p,(double*)d1,(double*)d2)
#define mbs_deBoorDer2C3d(degree,lastknot,knots,ctlpoints,t,p,d1,d2) \
  mbs_multideBoorDer2d(degree,lastknot,knots,1,3,0,(double*)ctlpoints,t, \
    (double*)p,(double*)d1,(double*)d2)
#define mbs_deBoorDer2C4d(degree,lastknot,knots,ctlpoints,t,p,d1,d2) \
  mbs_multideBoorDer2d(degree,lastknot,knots,1,4,0,(double*)ctlpoints,t, \
    (double*)p,(double*)d1,(double*)d2)


int mbs_multideBoorDer3d ( int degree, int lastknot, const double *knots,
                           int ncurves, int spdimen,
                           int pitch, const double *ctlpoints, double t,
                           double *p, double *d1, double *d2, double *d3 );

#define mbs_deBoorDer3C1d(degree,lastknot,knots,coeff,t,p,d1,d2,d3) \
  mbs_multideBoorDer3d(degree,lastknot,knots,1,1,0,coeff,t,p,d1,d2,d3)
#define mbs_deBoorDer3C2d(degree,lastknot,knots,ctlpoints,t,p,d1,d2,d3) \
  mbs_multideBoorDer3d(degree,lastknot,knots,1,2,0,(double*)ctlpoints,t, \
    (double*)p,(double*)d1,(double*)d2,(double*)d3)
#define mbs_deBoorDer3C3d(degree,lastknot,knots,ctlpoints,t,p,d1,d2,d3) \
  mbs_multideBoorDer3d(degree,lastknot,knots,1,3,0,(double*)ctlpoints,t, \
    (double*)p,(double*)d1,(double*)d2,(double*)d3)
#define mbs_deBoorDer3C4d(degree,lastknot,knots,ctlpoints,t,p,d1,d2,d3) \
  mbs_multideBoorDer3d(degree,lastknot,knots,1,4,0,(double*)ctlpoints,t, \
    (double*)p,(double*)d1,(double*)d2,(double*)d3)


boolean mbs_deBoorDerPd ( int degreeu, int lastknotu, const double *knotsu,    
                          int degreev, int lastknotv, const double *knotsv,
                          int spdimen, int pitch, const double *ctlpoints,
                          double u, double v,
                          double *ppoint, double *uder, double *vder );

boolean mbs_deBoorDer2Pd ( int degreeu, int lastknotu, const double *knotsu,
                           int degreev, int lastknotv, const double *knotsv,
                           int spdimen, int pitch, const double *ctlpoints, 
                           double u, double v,
                           double *ppoint, double *uder, double *vder,
                           double *uuder, double *uvder, double *vvder );

boolean mbs_deBoorDer3Pd ( int degreeu, int lastknotu, const double *knotsu,
                           int degreev, int lastknotv, const double *knotsv,
                           int spdimen, int pitch, const double *ctlpoints, 
                           double u, double v,
                           double *ppoint, double *uder, double *vder,
                           double *uuder, double *uvder, double *vvder,
                           double *uuuder, double *uuvder, double *uvvder,
                           double *vvvder );

 
int mbs_multiKnotInsd ( int degree, int *lastknot,
                        double *knots,
                        int ncurves, int spdimen, int inpitch, int outpitch,
                        double *ctlpoints, double t );

#define mbs_KnotInsC1d(degree,lastknot,knots,coeff,t) \
  mbs_multiKnotInsd (degree,lastknot,knots,1,1,0,0,coeff,t)
#define mbs_KnotInsC2d(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsd (degree,lastknot,knots,1,2,0,0,(double*)ctlpoints,t)
#define mbs_KnotInsC3d(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsd (degree,lastknot,knots,1,3,0,0,(double*)ctlpoints,t)
#define mbs_KnotInsC4d(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsd (degree,lastknot,knots,1,4,0,0,(double*)ctlpoints,t)


int mbs_multiKnotInsClosedd ( int degree, int *lastknot, double *knots,
                              int ncurves, int spdimen, int inpitch, int outpitch,
                              double *ctlpoints, double t );

#define mbs_KnotInsClosedC1d(degree,lastknot,knots,coeff,t) \
  mbs_multiKnotInsClosedd(degree,lastknot,knots,1,1,0,0,coeff,t)
#define mbs_KnotInsClosedC2d(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsClosedd(degree,lastknot,knots,1,2,0,0,(double*)ctlpoints,t)
#define mbs_KnotInsClosedC3d(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsClosedd(degree,lastknot,knots,1,3,0,0,(double*)ctlpoints,t)
#define mbs_KnotInsClosedC4d(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsClosedd(degree,lastknot,knots,1,4,0,0,(double*)ctlpoints,t)


int mbs_multiKnotRemoved ( int degree, int *lastknot,
                           double *knots,
                           int ncurves, int spdimen, int inpitch, int outpitch,
                           double *ctlpoints, int knotnum );

#define mbs_KnotRemoveC1d(degree,lastknot,knots,coeff,knotnum) \
  mbs_multiKnotRemoved(degree,lastknot,knots,1,1,0,0,coeff,knotnum)
#define mbs_KnotRemoveC2d(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemoved(degree,lastknot,knots,1,2,0,0, \
    (double*)ctlpoints,knotnum)
#define mbs_KnotRemoveC3d(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemoved(degree,lastknot,knots,1,3,0,0, \
    (double*)ctlpoints,knotnum)
#define mbs_KnotRemoveC4d(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemoved(degree,lastknot,knots,1,4,0,0, \
    (double*)ctlpoints,knotnum)


int mbs_multiKnotRemoveClosedd ( int degree, int *lastknot,
                                 double *knots,
                                 int ncurves, int spdimen, int inpitch, int outpitch,
                                 double *ctlpoints, int knotnum );

#define mbs_KnotRemoveClosedC1d(degree,lastknot,knots,coeff,knotnum) \
  mbs_multiKnotRemoveClosedd (degree,lastknot,knots,1,1,0,0,coeff,knotnum)
#define mbs_KnotRemoveClosedC2d(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemoveClosedd (degree,lastknot,knots,1,2,0,0, \
    (double*)ctlpoints,knotnum)
#define mbs_KnotRemoveClosedC3d(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemoveClosedd (degree,lastknot,knots,1,3,0,0, \
    (double*)ctlpoints,knotnum)
#define mbs_KnotRemoveClosedC4d(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemoveClosedd (degree,lastknot,knots,1,4,0,0, \
    (double*)ctlpoints,knotnum)


void mbs_multiRemoveSuperfluousKnotsd ( int ncurves, int spdimen, int degree,
                                        int *lastknot, double *knots,
                                        int inpitch, int outpitch,
                                        double *ctlpoints );


boolean mbs_multiMaxKnotInsd ( int ncurves, int spdimen, int degree,
                               int inlastknot, const double *inknots,
                               int inpitch, const double *inctlpoints,
                               int *outlastknot, double *outknots,
                               int outpitch, double *outctlpoints,
                               int *skipl, int *skipr );

#define mbs_MaxKnotInsC1d(degree,inlastknot,inknots,incoeff, \
                          outlastknot,outknots,outcoeff,skipl,skipr) \
  mbs_multiMaxKnotInsd(1,1,degree,inlastknot,inknots,0,incoeff, \
    outlastknot,outknots,0,outcoeff,skipl,skipr)
#define mbs_MaxKnotInsC2d(degree,inlastknot,inknots,inctlpoints, \
                          outlastknot,outknots,outctlpoints,skipl,skipr) \
  mbs_multiMaxKnotInsd(1,2,degree,inlastknot,inknots,0,(double*)inctlpoints, \
    outlastknot,outknots,0,(double*)outctlpoints,skipl,skipr)
#define mbs_MaxKnotInsC3d(degree,inlastknot,inknots,inctlpoints, \
                          outlastknot,outknots,outctlpoints,skipl,skipr) \
  mbs_multiMaxKnotInsd(1,3,degree,inlastknot,inknots,0,(double*)inctlpoints, \
    outlastknot,outknots,0,(double*)outctlpoints,skipl,skipr)
#define mbs_MaxKnotInsC4d(degree,inlastknot,inknots,inctlpoints, \
                          outlastknot,outknots,outctlpoints,skipl,skipr) \
  mbs_multiMaxKnotInsd(1,4,degree,inlastknot,inknots,0,(double*)inctlpoints, \
    outlastknot,outknots,0,(double*)outctlpoints,skipl,skipr)


boolean mbs_multiBSCurvesToBezd ( int spdimen, int ncurves,
                                  int degree, int lastinknot, const double *inknots,
                                  int inpitch, const double *inctlp,
                                  int *kpcs, int *lastoutknot, double *outknots,
                                  int outpitch, double *outctlp );

#define mbs_BSToBezC1d(degree,lastinknot,inknots,incoeff,kpcs, \
    lastoutknot,outknots,outcoeff) \
  mbs_multiBSCurvesToBezd(1,1,degree,lastinknot,inknots,0,incoeff, \
    kpcs,lastoutknot,outknots,0,outcoeff)
#define mbs_BSToBezC2d(degree,lastinknot,inknots,inctlp,kpcs, \
    lastoutknot,outknots,outctlp) \
  mbs_multiBSCurvesToBezd(2,1,degree,lastinknot,inknots,0,(double*)inctlp, \
    kpcs,lastoutknot,outknots,0,(double*)outctlp)
#define mbs_BSToBezC3d(degree,lastinknot,inknots,inctlp,kpcs, \
    lastoutknot,outknots,outctlp) \
  mbs_multiBSCurvesToBezd(3,1,degree,lastinknot,inknots,0,(double*)inctlp, \
    kpcs,lastoutknot,outknots,0,(double*)outctlp)
#define mbs_BSToBezC4d(degree,lastinknot,inknots,inctlp,kpcs, \
    lastoutknot,outknots,outctlp) \
  mbs_multiBSCurvesToBezd(4,1,degree,lastinknot,inknots,0,(double*)inctlp, \
    kpcs,lastoutknot,outknots,0,(double*)outctlp)


boolean mbs_BSPatchToBezd ( int spdimen,
                            int degreeu, int lastuknot, const double *uknots,
                            int degreev, int lastvknot, const double *vknots,
                            int inpitch, const double *inctlp,
                            int *kupcs, int *lastoutuknot, double *outuknots,
                            int *kvpcs, int *lastoutvknot, double *outvknots,
                            int outpitch, double *outctlp );


int mbs_NumKnotIntervalsd ( int degree, int lastknot, const double *knots );

int mbs_LastknotMaxInsd ( int degree, int lastknot, const double *knots, 
                          int *numknotintervals );

int mbs_BSProdRepSized ( int degree1, int lastknot1, const double *knots1,
                          int degree2, int lastknot2, const double *knots2 );

int mbs_NumMaxKnotsd ( int degree, int lastknot, const double *knots );

void mbs_SetBSProdKnotsd ( int degree1, int lastknot1, const double *knots1,
                           int degree2, int lastknot2, const double *knots2,
                           int *degree, int *lastknot, double *knots );

void mbs_SetKnotPatternd ( int lastinknot, const double *inknots,
                           int multipl,
                           int *lastoutknot, double *outknots );

void mbs_multiBezScaled ( int degree, int narcs, int ncurves, int spdimen,  
                          int pitch, double *ctlpoints );

void mbs_multiBezUnscaled ( int degree, int narcs, int ncurves, int spdimen,  
                            int pitch, double *ctlpoints );

boolean mbs_multiMultBezCd ( int nscf, int degscf, int scfpitch,
                             const double *scfcoef,
                             int spdimen,
                             int nvecf, int degvecf, int vecfpitch,
                             const double *vecfcp,
                             int *degprod, int prodpitch, double *prodcp );

boolean mbs_multiMultBSCd ( int nscf, int degscf,
                            int scflastknot, const double *scfknots,
                            int scfpitch, const double *scfcoef,
                            int spdimen,
                            int nvecf, int degvecf,
                            int vecflastknot, const double *vecfknots,
                            int vecfpitch, const double *vecfcp,
                            int *degprod, int *prodlastknot, double *prodknots,
                            int prodpitch, double *prodcp );


boolean mbs_multiBCDegElevd ( int ncurves, int spdimen,
                              int inpitch, int indegree, const double *inctlpoints,
                              int deltadeg,
                              int outpitch, int *outdegree, double *outctlpoints );

#define mbs_BCDegElevC1d(indegree,incoeff,deltadeg,outdegree,outcoeff) \
  mbs_multiBCDegElevd ( 1, 1, 0, indegree, incoeff, deltadeg, \
    0, outdegree, outcoeff )
#define mbs_BCDegElevC2d(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegElevd ( 1, 2, 0, indegree, (double*)inctlpoints, deltadeg, \
    0, outdegree, (double*)outctlpoints )
#define mbs_BCDegElevC3d(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegElevd ( 1, 3, 0, indegree, (double*)inctlpoints, deltadeg, \
    0, outdegree, (double*)outctlpoints)
#define mbs_BCDegElevC4d(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegElevd ( 1, 4, 0, indegree, (double*)inctlpoints, deltadeg, \
    0, outdegree, (double*)outctlpoints )


boolean mbs_BCDegElevPd ( int spdimen,
                          int indegreeu, int indegreev, const double *inctlp,
                          int deltadegu, int deltadegv,
                          int *outdegreeu, int *outdegreev,
                          double *outctlp );

#define mbs_BCDegElevP1d(indegreeu,indegreev,incoeff,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outcoeff) \
  mbs_BCDegElevPd ( 1, indegreeu, indegreev, incoeff, deltadegu, deltadegv, \
    outdegreeu, outdegreev, outcoeff )
#define mbs_BCDegElevP2d(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegElevPd ( 2, indegreeu, indegreev, (double*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (double*)outctlp )
#define mbs_BCDegElevP3d(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegElevPd ( 3, indegreeu, indegreev, (double*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (double*)outctlp )
#define mbs_BCDegElevP4d(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegElevPd ( 4, indegreeu, indegreev, (double*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (double*)outctlp )


void mbs_multiBCDegRedd ( int ncurves, int spdimen,
                          int inpitch, int indegree, const double *inctlpoints,
                          int deltadeg,
                          int outpitch, int *outdegree, double *outctlpoints );

#define mbs_BCDegRedC1d(indegree,incoeff,deltadeg,outdegree,outcoeff) \
  mbs_multiBCDegRedd ( 1, 1, 0, indegree, incoeff, deltadeg, \
    0, outdegree, outcoeff )
#define mbs_BCDegRedC2d(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegRedd ( 1, 2, 0, indegree, (double*)inctlpoints, deltadeg, \
    0, outdegree, (double*)outctlpoints )
#define mbs_BCDegRedC3d(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegRedd ( 1, 3, 0, indegree, (double*)inctlpoints, deltadeg, \
    0, outdegree, (double*)outctlpoints)
#define mbs_BCDegRedC4d(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegRedd ( 1, 4, 0, indegree, (double*)inctlpoints, deltadeg, \
    0, outdegree, (double*)outctlpoints )


boolean mbs_BCDegRedPd ( int spdimen,
                         int indegreeu, int indegreev, const double *inctlp,
                         int deltadegu, int deltadegv,
                         int *outdegreeu, int *outdegreev,
                         double *outctlp );

#define mbs_BCDegRedP1d(indegreeu,indegreev,incoeff,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outcoeff) \
  mbs_BCDegRedPd ( 1, indegreeu, indegreev, incoeff, deltadegu, deltadegv, \
    outdegreeu, outdegreev, outcoeff )
#define mbs_BCDegRedP2d(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegRedPd ( 2, indegreeu, indegreev, (double*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (double*)outctlp )
#define mbs_BCDegRedP3d(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegRedPd ( 3, indegreeu, indegreev, (double*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (double*)outctlp )
#define mbs_BCDegRedP4d(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegRedPd ( 4, indegreeu, indegreev, (double*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (double*)outctlp )


boolean mbs_multiBSDegElevd ( int ncurves, int spdimen,
                              int indegree, int inlastknot, const double *inknots,
                              int inpitch, const double *inctlpoints,
                              int deltadeg,
                              int *outdegree, int *outlastknot,
                              double *outknots, int outpitch, double *outctlpoints,
                              boolean freeend );

#define mbs_BSDegElevC1d(indegree,inlastknot,inknots,incoeff, \
    deltadeg,outdegree,outlastknot,outknots,outcoeff,freeend) \
  mbs_multiBSDegElevd(1,1,indegree,inlastknot,inknots,0,incoeff, \
    deltadeg,outdegree,outlastknot,outknots,0,outcoeff,freeend)
#define mbs_BSDegElevC2d(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints,freeend) \
  mbs_multiBSDegElevd(1,2,indegree,inlastknot,inknots,0,(double*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(double*)outctlpoints,freeend)
#define mbs_BSDegElevC3d(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints,freeend) \
  mbs_multiBSDegElevd(1,3,indegree,inlastknot,inknots,0,(double*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(double*)outctlpoints,freeend)
#define mbs_BSDegElevC4d(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints,freeend) \
  mbs_multiBSDegElevd(1,4,indegree,inlastknot,inknots,0,(double*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(double*)outctlpoints,freeend)

boolean mbs_multiBSDegElevClosedd ( int ncurves, int spdimen,
                         int indegree, int inlastknot, const double *inknots,
                         int inpitch, const double *inctlpoints,
                         int deltadeg,
                         int *outdegree, int *outlastknot,
                         double *outknots, int outpitch, double *outctlpoints );

#define mbs_BSDegElevClosedC1d(indegree,inlastknot,inknots,incoeff, \
    deltadeg,outdegree,outlastknot,outknots,outcoeff) \
  mbs_multiBSDegElevClosedd(1,1,indegree,inlastknot,inknots,0,incoeff, \
    deltadeg,outdegree,outlastknot,outknots,0,outcoeff)
#define mbs_BSDegElevClosedC2d(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiBSDegElevClosedd(1,2,indegree,inlastknot,inknots,0,(double*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(double*)outctlpoints)
#define mbs_BSDegElevClosedC3d(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiBSDegElevClosedd(1,3,indegree,inlastknot,inknots,0,(double*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(double*)outctlpoints)
#define mbs_BSDegElevClosedC4d(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiBSDegElevClosedd(1,4,indegree,inlastknot,inknots,0,(double*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(double*)outctlpoints)


boolean mbs_multiBSDegRedd ( int ncurves, int spdimen,
                             int indegree, int inlastknot, const double *inknots,
                             int inpitch, CONST_ double *inctlpoints,
                             int deltadeg,
                             int *outdegree, int *outlastknot, double *outknots,
                             int outpitch, double *outctlpoints );

#define mbs_BSDegRedC1d(indegree,inlastknot,inknots,incoeff,deltadeg, \
    outdegree,outlastknot,outknots,outcoeff) \
  mbs_multiBSDegRedd(1,1,indegree,inlastknot,inknots,0,incoeff,deltadeg, \
    outdegree,outlastknot,outknots,0,outcoeff)
#define mbs_BSDegRedC2d(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedd(1,2,indegree,inlastknot,inknots,0,(double*)incpoints,deltadeg, \
    outdegree,outlastknot,outknots,0,(double*)outcpoints)
#define mbs_BSDegRedC3d(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedd(1,3,indegree,inlastknot,inknots,0,(double*)incpoints,deltadeg, \
    outdegree,outlastknot,outknots,0,(double*)outcpoints)
#define mbs_BSDegRedC4d(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedd(1,4,indegree,inlastknot,inknots,0,(double*)incpoints,deltadeg, \
    outdegree,outlastknot,outknots,0,(double*)outcpoints)


boolean mbs_multiBSDegRedClosedd ( int ncurves, int spdimen,
                             int indegree, int inlastknot, const double *inknots,
                             int inpitch, CONST_ double *inctlpoints,
                             int deltadeg,
                             int *outdegree, int *outlastknot, double *outknots,
                             int outpitch, double *outctlpoints );

#define mbs_BSDegRedClosedC1d(indegree,inlastknot,inknots,incoeff,deltadeg, \
    outdegree,outlastknot,outknots,outcoeff) \
  mbs_multiBSDegRedClosedd(1,1,indegree,inlastknot,inknots,0,incoeff,deltadeg, \
    outdegree,outlastknot,outknots,0,outcoeff)
#define mbs_BSDegRedClosedC2d(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedClosedd(1,2,indegree,inlastknot,inknots,0,(double*)incpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(double*)outcpoints)
#define mbs_BSDegRedClosedC3d(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedClosedd(1,3,indegree,inlastknot,inknots,0,(double*)incpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(double*)outcpoints)
#define mbs_BSDegRedClosedC4d(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedClosedd(1,4,indegree,inlastknot,inknots,0,(double*)incpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(double*)outcpoints)


boolean mbs_multiBCHornerd ( int degree, int ncurves, int spdimen, int pitch,
                             const double *ctlpoints, double t, double *cpoints );

#define mbs_BCHornerC1d(degree,coeff,t,value) \
  mbs_multiBCHornerd ( degree, 1, 1, 0, coeff, t, value )
#define mbs_BCHornerC2d(degree,ctlpoints,t,cpoint) \
  mbs_multiBCHornerd ( degree, 1, 2, 0, (double*)ctlpoints, t, (double*)cpoint )
#define mbs_BCHornerC3d(degree,ctlpoints,t,cpoint) \
  mbs_multiBCHornerd ( degree, 1, 3, 0, (double*)ctlpoints, t, (double*)cpoint )
#define mbs_BCHornerC4d(degree,ctlpoints,t,cpoint) \
  mbs_multiBCHornerd ( degree, 1, 4, 0, (double*)ctlpoints, t, (double*)cpoint )

boolean mbs_BCHornerC2Rd ( int degree, const point3d *ctlpoints, double t,
                           point2d *cpoint );

boolean mbs_BCHornerC3Rd ( int degree, const point4d *ctlpoints, double t,
                           point3d *cpoint );

boolean _mbs_BCHornerPd ( int degreeu, int degreev, int spdimen,
                          const double *ctlpoints,
                          double u, double v, double *ppoint,
                          double *workspace );

#define _mbs_BCHornerP1d(degreeu,degreev,coeff,u,v,ppoint,ws) \
  _mbs_BCHornerPd ( degreeu, degreev, 1, coeff, u, v, ppoint, ws )
#define _mbs_BCHornerP2d(degreeu,degreev,ctlpoints,u,v,ppoint,ws) \
  _mbs_BCHornerPd ( degreeu, degreev, 2, (double*)ctlpoints, \
    u, v, (double*)ppoint, ws )
#define _mbs_BCHornerP3d(degreeu,degreev,ctlpoints,u,v,ppoint,ws) \
  _mbs_BCHornerPd ( degreeu, degreev, 3, (double*)ctlpoints, \
    u, v, (double*)ppoint, ws )
#define _mbs_BCHornerP4d(degreeu,degreev,ctlpoints,u,v,ppoint,ws) \
  _mbs_BCHornerPd ( degreeu, degreev, 4, (double*)ctlpoints, \
    u, v, (double*)ppoint, ws)

boolean mbs_BCHornerPd ( int degreeu, int degreev, int spdimen,
                         const double *ctlpoints,
                         double u, double v, double *ppoint );

#define mbs_BCHornerP1d(degreeu,degreev,coeff,u,v,ppoint) \
  mbs_BCHornerPd ( degreeu, degreev, 1, coeff, u, v, ppoint )
#define mbs_BCHornerP2d(degreeu,degreev,ctlpoints,u,v,ppoint ) \
  mbs_BCHornerPd ( degreeu, degreev, 2, (double*)ctlpoints, \
    u, v, (double*)ppoint )
#define mbs_BCHornerP3d(degreeu,degreev,ctlpoints,u,v,ppoint ) \
  mbs_BCHornerPd ( degreeu, degreev, 3, (double*)ctlpoints, \
    u, v, (double*)ppoint )
#define mbs_BCHornerP4d(degreeu,degreev,ctlpoints,u,v,ppoint ) \
  mbs_BCHornerPd ( degreeu, degreev, 4, (double*)ctlpoints, \
    u, v, (double*)ppoint )

boolean _mbs_BCHornerP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                            double u, double v, point3d *p, double *workspace );

boolean mbs_BCHornerP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                           double u, double v, point3d *p );


boolean _mbs_multiBCHornerDerd ( int degree, int ncurves, int spdimen, int pitch,
                                 const double *ctlpoints,
                                 double t, double *p, double *d,
                                 double *workspace );

#define _mbs_BCHornerDerC1d(degree,coeff,t,p,d,workspace) \
  _mbs_multiBCHornerDerd ( degree, 1, 1, 0, coeff, t, p, d, workspace )
#define _mbs_BCHornerDerC2d(degree,ctlpoints,t,p,d,workspace) \
  _mbs_multiBCHornerDerd ( degree, 1, 2, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d, workspace )
#define _mbs_BCHornerDerC3d(degree,ctlpoints,t,p,d,workspace) \
  _mbs_multiBCHornerDerd ( degree, 1, 3, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d, workspace )
#define _mbs_BCHornerDerC4d(degree,ctlpoints,t,p,d,workspace) \
  _mbs_multiBCHornerDerd ( degree, 1, 4, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d, workspace )

boolean mbs_multiBCHornerDerd ( int degree, int ncurves, int spdimen, int pitch,
                                const double *ctlpoints,
                                double t, double *p, double *d );

#define mbs_BCHornerDerC1d(degree,coeff,t,p,d) \
  mbs_multiBCHornerDerd ( degree, 1, 1, 0, coeff, t, p, d )
#define mbs_BCHornerDerC2d(degree,ctlpoints,t,p,d) \
  mbs_multiBCHornerDerd ( degree, 1, 2, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d )
#define mbs_BCHornerDerC3d(degree,ctlpoints,t,p,d) \
  mbs_multiBCHornerDerd ( degree, 1, 3, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d )
#define mbs_BCHornerDerC4d(degree,ctlpoints,t,p,d) \
  mbs_multiBCHornerDerd ( degree, 1, 4, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d )

boolean mbs_BCHornerDerC2Rd ( int degree, const point3d *ctlpoints, double t,
                              point2d *p, vector2d *d );
boolean mbs_BCHornerDerC3Rd ( int degree, const point4d *ctlpoints, double t,
                              point3d *p, vector3d *d );


boolean _mbs_BCHornerDerPd ( int degreeu, int degreev, int spdimen,
                             const double *ctlpoints,
                             double u, double v,  
                             double *p, double *du, double *dv,
                             double *workspace );

#define _mbs_BCHornerDerP1d(degreeu,degreev,coeff,u,v,p,du,dv,workspace) \
  _mbs_BCHornerDerPd ( degreeu, degreev, 1, coeff, u, v, p, du, dv, workspace )
#define _mbs_BCHornerDerP2d(degreeu,degreev,ctlpoints,u,v,p,du,dv,workspace) \
  _mbs_BCHornerDerPd ( degreeu, degreev, 2, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv, workspace )
#define _mbs_BCHornerDerP3d(degreeu,degreev,ctlpoints,u,v,p,du,dv,workspace) \
  _mbs_BCHornerDerPd ( degreeu, degreev, 3, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv, workspace )
#define _mbs_BCHornerDerP4d(degreeu,degreev,ctlpoints,u,v,p,du,dv,workspace) \
  _mbs_BCHornerDerPd ( degreeu, degreev, 4, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv, workspace )

boolean mbs_BCHornerDerPd ( int degreeu, int degreev, int spdimen,
                            const double *ctlpoints,
                            double u, double v,  
                            double *p, double *du, double *dv );

#define mbs_BCHornerDerP1d(degreeu,degreev,coeff,u,v,p,du,dv) \
  mbs_BCHornerDerPd ( degreeu, degreev, 1, coeff, u, v, p, du, dv )
#define mbs_BCHornerDerP2d(degreeu,degreev,ctlpoints,u,v,p,du,dv) \
  mbs_BCHornerDerPd ( degreeu, degreev, 2, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv )
#define mbs_BCHornerDerP3d(degreeu,degreev,ctlpoints,u,v,p,du,dv) \
  mbs_BCHornerDerPd ( degreeu, degreev, 3, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv )
#define mbs_BCHornerDerP4d(degreeu,degreev,ctlpoints,u,v,p,du,dv) \
  mbs_BCHornerDerPd ( degreeu, degreev, 4, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv )

boolean _mbs_BCHornerDerP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                               double u, double v,
                               point3d *p, vector3d *du, vector3d *dv,
                               double *workspace );

boolean mbs_BCHornerDerP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                              double u, double v,
                              point3d *p, vector3d *du, vector3d *dv );

boolean _mbs_BCHornerNvP3d ( int degreeu, int degreev, const point3d *ctlpoints,
                             double u, double v,
                             point3d *p, vector3d *nv, double *workspace );

boolean mbs_BCHornerNvP3d ( int degreeu, int degreev, const point3d *ctlpoints,
                            double u, double v,
                            point3d *p, vector3d *nv );

boolean _mbs_BCHornerNvP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                              double u, double v,
                              point3d *p, vector3d *nv, double *workspace );

boolean mbs_BCHornerNvP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                             double u, double v,
                             point3d *p, vector3d *nv );


boolean _mbs_multiBCHornerDer2d ( int degree, int ncurves, int spdimen, int pitch,
                               const double *ctlpoints, 
                               double t, double *p, double *d1, double *d2,
                               double *workspace );

#define _mbs_BCHornerDer2C1d(degree,coeff,t,p,d1,d2,workspace) \
  _mbs_multiBCHornerDer2d ( degree, 1, 1, 0, coeff, t, p, d1, d2, workspace )
#define _mbs_BCHornerDer2C2d(degree,ctlpoints,t,p,d1,d2,workspace) \
  _mbs_multiBCHornerDer2d ( degree, 1, 2, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2, workspace )
#define _mbs_BCHornerDer2C3d(degree,ctlpoints,t,p,d1,d2,workspace) \
  _mbs_multiBCHornerDer2d ( degree, 1, 3, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2, workspace )
#define _mbs_BCHornerDer2C4d(degree,ctlpoints,t,p,d1,d2,workspace) \
  _mbs_multiBCHornerDer2d ( degree, 1, 4, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2, workspace )

boolean mbs_multiBCHornerDer2d ( int degree, int ncurves, int spdimen, int pitch,
                              const double *ctlpoints, 
                              double t, double *p, double *d1, double *d2 );

#define mbs_BCHornerDer2C1d(degree,coeff,t,p,d1,d2) \
  mbs_multiBCHornerDer2d ( degree, 1, 1, 0, coeff, t, p, d1, d2 )
#define mbs_BCHornerDer2C2d(degree,ctlpoints,t,p,d1,d2) \
  mbs_multiBCHornerDer2d ( degree, 1, 2, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2 )
#define mbs_BCHornerDer2C3d(degree,ctlpoints,t,p,d1,d2) \
  mbs_multiBCHornerDer2d ( degree, 1, 3, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2 )
#define mbs_BCHornerDer2C4d(degree,ctlpoints,t,p,d1,d2) \
  mbs_multiBCHornerDer2d ( degree, 1, 4, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2 )

boolean _mbs_BCHornerDer2C2Rd ( int degree, const point3d *ctlpoints, double t,
                                point2d *p, vector2d *d1, vector2d *d2,
                                double *workspace );
boolean mbs_BCHornerDer2C2Rd ( int degree, const point3d *ctlpoints, double t,
                               point2d *p, vector2d *d1, vector2d *d2 );
boolean _mbs_BCHornerDer2C3Rd ( int degree, const point4d *ctlpoints, double t,
                                point3d *p, vector3d *d1, vector3d *d2,
                                double *workspace );
boolean mbs_BCHornerDer2C3Rd ( int degree, const point4d *ctlpoints, double t,
                               point3d *p, vector3d *d1, vector3d *d2 );

boolean _mbs_BCHornerDer2Pd ( int degreeu, int degreev, int spdimen,
                             const double *ctlpoints,
                             double u, double v,
                             double *p, double *du, double *dv,
                             double *duu, double *duv, double *dvv,
                             double *workspace );

#define _mbs_BCHornerDer2P1d(degreeu,degreev,coeff,u,v,p,du,dv,duu,duv,dvv,workspace) \
  _mbs_BCHornerDer2Pd ( degreeu, degreev, 1, coeff, u, v, \
    p, du, dv, duu, duv, dvv, workspace )
#define _mbs_BCHornerDer2P2d(degreeu,degreev,ctlpoints,u,v,p,du,dv,duu,duv,dvv,workspace) \
  _mbs_BCHornerDer2Pd ( degreeu, degreev, 2, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv, (double*)duu, (double*)duv, (double*)dvv, workspace )
#define _mbs_BCHornerDer2P3d(degreeu,degreev,ctlpoints,u,v,p,du,dv,duu,duv,dvv,workspace) \
  _mbs_BCHornerDer2Pd ( degreeu, degreev, 3, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv, (double*)duu, (double*)duv, (double*)dvv, workspace )
#define _mbs_BCHornerDer2P4d(degreeu,degreev,ctlpoints,u,v,p,du,dv,duu,duv,dvv,workspace) \
  _mbs_BCHornerDer2Pd ( degreeu, degreev, 4, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv, (double*)duu, (double*)duv, (double*)dvv, workspace )

boolean mbs_BCHornerDer2Pd ( int degreeu, int degreev, int spdimen,
                             const double *ctlpoints,
                             double u, double v,
                             double *p, double *du, double *dv,
                             double *duu, double *duv, double *dvv );

#define mbs_BCHornerDer2P1d(degreeu,degreev,coeff,u,v,p,du,dv,duu,duv,dvv) \
  mbs_BCHornerDer2Pd ( degreeu, degreev, 1, coeff, u, v, \
    p, du, dv, duu, duv, dvv )
#define mbs_BCHornerDer2P2d(degreeu,degreev,ctlpoints,u,v,p,du,dv,duu,duv,dvv) \
  mbs_BCHornerDer2Pd ( degreeu, degreev, 2, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv, (double*)duu, (double*)duv, (double*)dvv )
#define mbs_BCHornerDer2P3d(degreeu,degreev,ctlpoints,u,v,p,du,dv,duu,duv,dvv) \
  mbs_BCHornerDer2Pd ( degreeu, degreev, 3, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv, (double*)duu, (double*)duv, (double*)dvv )
#define mbs_BCHornerDer2P4d(degreeu,degreev,ctlpoints,u,v,p,du,dv,duu,duv,dvv) \
  mbs_BCHornerDer2Pd ( degreeu, degreev, 4, (double*)ctlpoints, u, v, \
    (double*)p, (double*)du, (double*)dv, (double*)duu, (double*)duv, (double*)dvv )

boolean _mbs_BCHornerDer2P3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                                double u, double v,
                                point3d *p, vector3d *du, vector3d *dv,
                                vector3d *duu, vector3d *duv, vector3d *dvv,
                                double *workspace );
boolean mbs_BCHornerDer2P3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                               double u, double v,
                               point3d *p, vector3d *du, vector3d *dv,
                               vector3d *duu, vector3d *duv, vector3d *dvv );


boolean _mbs_multiBCHornerDer3d ( int degree, int ncurves, int spdimen, int pitch,
                                  const double *ctlpoints, double t,
                                  double *p, double *d1, double *d2, double *d3,
                                  double *workspace );

#define _mbs_BCHornerDer3C1d(degree,coeff,t,p,d1,d2,d3,workspace) \
  _mbs_multiBCHornerDer3d ( degree, 1, 1, 0, coeff, t, p, d1, d2, d3, workspace )
#define _mbs_BCHornerDer3C2d(degree,ctlpoints,t,p,d1,d2,d3,workspace) \
  _mbs_multiBCHornerDer3d ( degree, 1, 2, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2, (double*)d3, workspace )
#define _mbs_BCHornerDer3C3d(degree,ctlpoints,t,p,d1,d2,d3,workspace) \
  _mbs_multiBCHornerDer3d ( degree, 1, 3, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2, (double*)d3, workspace )
#define _mbs_BCHornerDer3C4d(degree,ctlpoints,t,p,d1,d2,d3,workspace) \
  _mbs_multiBCHornerDer3d ( degree, 1, 4, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2, (double*)d3, workspace )

boolean mbs_multiBCHornerDer3d ( int degree, int ncurves, int spdimen, int pitch,
                                 const double *ctlpoints, double t,
                                 double *p, double *d1, double *d2, double *d3 );

#define mbs_BCHornerDer3C1d(degree,coeff,t,p,d1,d2,d3) \
  mbs_multiBCHornerDer3d ( degree, 1, 1, 0, coeff, t, p, d1, d2, d3 )
#define mbs_BCHornerDer3C2d(degree,ctlpoints,t,p,d1,d2,d3) \
  mbs_multiBCHornerDer3d ( degree, 1, 2, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2, (double*)d3 )
#define mbs_BCHornerDer3C3d(degree,ctlpoints,t,p,d1,d2,d3) \
  mbs_multiBCHornerDer3d ( degree, 1, 3, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2, (double*)d3 )
#define mbs_BCHornerDer3C4d(degree,ctlpoints,t,p,d1,d2,d3) \
  mbs_multiBCHornerDer3d ( degree, 1, 4, 0, (double*)ctlpoints, t, \
    (double*)p, (double*)d1, (double*)d2, (double*)d3 )


boolean _mbs_FindBezPatchDiagFormd ( int degreeu, int degreev, int spdimen,
                                     CONST_ double *cpoints,
                                     int k, int l, double u, double v,
                                     double *dfcp, double *workspace );

boolean mbs_FindBezPatchDiagFormd ( int degreeu, int degreev, int spdimen,
                                    CONST_ double *cpoints,
                                    int k, int l, double u, double v,
                                    double *dfcp );

boolean _mbs_BCHornerDer3Pd ( int degreeu, int degreev, int spdimen,
                              CONST_ double *ctlpoints,
                              double u, double v,
                              double *p, double *pu, double *pv,
                              double *puu, double *puv, double *pvv,
                              double *puuu, double *puuv, double *puvv, double *pvvv,
                              double *workspace );

#define _mbs_BCHornerDer3P1d(degreeu,degreev,coeff,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv,workspace) \
  _mbs_BCHornerDer3Pd ( degreeu, degreev, 1, coeff, u, v, \
    p, pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv, workspace )
#define _mbs_BCHornerDer3P2d(degreeu,degreev,ctlpoints,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv,workspace) \
  _mbs_BCHornerDer3Pd ( degreeu, degreev, 2, (double*)ctlpoints, u, v, \
    (double*)p, (double*)pu, (double*)pv, (double*)puu, (double*)puv, \
    (double*)pvv, (double*)puuu, (double*)puuv, (double*)puvv, (double*)pvvv, workspace )
#define _mbs_BCHornerDer3P3d(degreeu,degreev,ctlpoints,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv,workspace) \
  _mbs_BCHornerDer3Pd ( degreeu, degreev, 3, (double*)ctlpoints, u, v, \
    (double*)p, (double*)pu, (double*)pv, (double*)puu, (double*)puv, \
    (double*)pvv, (double*)puuu, (double*)puuv, (double*)puvv, (double*)pvvv, workspace )
#define _mbs_BCHornerDer3P4d(degreeu,degreev,ctlpoints,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv,workspace) \
  _mbs_BCHornerDer3Pd ( degreeu, degreev, 4, (double*)ctlpoints, u, v, \
    (double*)p, (double*)pu, (double*)pv, (double*)puu, (double*)puv, \
    (double*)pvv, (double*)puuu, (double*)puuv, (double*)puvv, (double*)pvvv, workspace )

boolean mbs_BCHornerDer3Pd ( int degreeu, int degreev, int spdimen,
                             CONST_ double *ctlpoints,
                             double u, double v,
                             double *p, double *pu, double *pv,
                             double *puu, double *puv, double *pvv,
                             double *puuu, double *puuv, double *puvv, double *pvvv );

#define mbs_BCHornerDer3P1d(degreeu,degreev,coeff,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv) \
  mbs_BCHornerDer3Pd ( degreeu, degreev, 1, coeff, u, v, \
    p, pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv )
#define mbs_BCHornerDer3P2d(degreeu,degreev,ctlpoints,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv) \
  mbs_BCHornerDer3Pd ( degreeu, degreev, 2, (double*)ctlpoints, u, v, \
    (double*)p, (double*)pu, (double*)pv, (double*)puu, (double*)puv, \
    (double*)pvv, (double*)puuu, (double*)puuv, (double*)puvv, (double*)pvvv )
#define mbs_BCHornerDer3P3d(degreeu,degreev,ctlpoints,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv) \
  mbs_BCHornerDer3Pd ( degreeu, degreev, 3, (double*)ctlpoints, u, v, \
    (double*)p, (double*)pu, (double*)pv, (double*)puu, (double*)puv, \
    (double*)pvv, (double*)puuu, (double*)puuv, (double*)puvv, (double*)pvvv )
#define mbs_BCHornerDer3P4d(degreeu,degreev,ctlpoints,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv) \
  mbs_BCHornerDer3Pd ( degreeu, degreev, 4, (double*)ctlpoints, u, v, \
    (double*)p, (double*)pu, (double*)pv, (double*)puu, (double*)puv, \
    (double*)pvv, (double*)puuu, (double*)puuv, (double*)puvv, (double*)pvvv )


void mbs_deBoorBasisd ( int degree, int lastknot, const double *knots,
                        double t, int *fnz, int *nnz, double *bfv );


boolean mbs_multiBSCubicInterpd ( int lastinterpknot, double *interpknots,
                                  int ncurves, int spdimen, int xpitch,
                                  const double *x,
                                  int ypitch,
                                  char bcl, const double *ybcl,
                                  char bcr, const double *ybcr,
                                  int *lastbsknot, double *bsknots,
                                  int bspitch, double *ctlpoints );

boolean mbs_multiBSCubicClosedInterpd ( int lastinterpknot, double *interpknots,
                               int ncurves, int spdimen, int xpitch,
                               const double *x,
                               int *lastbsknot, double *bsknots,
                               int bspitch, double *ctlpoints );


boolean mbs_BCFrenetC2d ( int degree, const point2d *ctlpoints, double t,
                          point2d *cpoint, vector2d *fframe, double *curvature );

boolean mbs_BCFrenetC2Rd ( int degree, const point3d *ctlpoints, double t,
                           point2d *cpoint, vector2d *fframe, double *curvature );

boolean mbs_BCFrenetC3d ( int degree, const point3d *ctlpoints, double t,
                          point3d *cpoint, vector3d *fframe, double *curvatures );

boolean mbs_BCFrenetC3Rd ( int degree, const point4d *ctlpoints, double t,
                           point3d *cpoint, vector3d *fframe, double *curvatures );


boolean _mbs_FundFormsBP3d ( int degreeu, int degreev, const point3d *ctlpoints,
                             double u, double v,
                             double *firstform, double *secondform,
                             double *workspace );
boolean mbs_FundFormsBP3d ( int degreeu, int degreev, const point3d *ctlpoints,
                            double u, double v,
                            double *firstform, double *secondform );

boolean _mbs_GMCurvaturesBP3d ( int degreeu, int degreev, const point3d *ctlpoints,
                                double u, double v,
                                double *gaussian, double *mean, double *workspace );
boolean mbs_GMCurvaturesBP3d ( int degreeu, int degreev, const point3d *ctlpoints,
                               double u, double v,
                               double *gaussian, double *mean );

boolean _mbs_PrincipalDirectionsBP3d ( int degreeu, int degreev,
                                       const point3d *ctlpoints,
                                       double u, double v,
                                       double *k1, vector2d *v1,
                                       double *k2, vector2d *v2,
                                       double *workspace );
boolean mbs_PrincipalDirectionsBP3d ( int degreeu, int degreev,
                                      const point3d *ctlpoints,
                                      double u, double v,
                                      double *k1, vector2d *v1,
                                      double *k2, vector2d *v2 );

boolean _mbs_FundFormsBP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                              double u, double v,
                              double *firstform, double *secondform,
                              double *workspace );
boolean mbs_FundFormsBP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                             double u, double v,
                             double *firstform, double *secondform );

boolean _mbs_GMCurvaturesBP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                                 double u, double v,
                                 double *gaussian, double *mean, double *workspace );
boolean mbs_GMCurvaturesBP3Rd ( int degreeu, int degreev, const point4d *ctlpoints,
                                double u, double v,
                                double *gaussian, double *mean );

boolean _mbs_PrincipalDirectionsBP3Rd ( int degreeu, int degreev,
                                        const point4d *ctlpoints,
                                        double u, double v,
                                        double *k1, vector2d *v1,
                                        double *k2, vector2d *v2,
                                        double *workspace );
boolean mbs_PrincipalDirectionsBP3Rd ( int degreeu, int degreev,
                                       const point4d *ctlpoints,
                                       double u, double v,
                                       double *k1, vector2d *v1,
                                       double *k2, vector2d *v2 );


void mbs_multiBisectBezCurvesd ( int degree, int ncurves,
                                 int spdimen, int pitch,
                                 double *ctlp, double *ctlq );

#define mbs_BisectBC1d(degree,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degree, 1, 1, 0, ctlp, ctlq )
#define mbs_BisectBC2d(degree,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degree, 1, 2, 0, (double*)ctlp, (double*)ctlq )
#define mbs_BisectBC3d(degree,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degree, 1, 3, 0, (double*)ctlp, (double*)ctlq )
#define mbs_BisectBC4d(degree,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degree, 1, 4, 0, (double*)ctlp, (double*)ctlq )

#define mbs_BisectBP1ud(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degreeu, 1, (degreev+1), 0, ctlp, ctlq )
#define mbs_BisectBP1vd(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degreev, degreeu+1, 1, degreev+1, ctlp, ctlq )
#define mbs_BisectBP2ud(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degreeu, 1, 2*(degreev+1), 0, \
    (double*)ctlp, (double*)ctlq )
#define mbs_BisectBP2vd(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degreev, degreeu+1, 2, 2*(degreev+1), \
    (double*)ctlp, (double*)ctlq )
#define mbs_BisectBP3ud(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degreeu, 1, 3*(degreev+1), 0, \
    (double*)ctlp, (double*)ctlq )
#define mbs_BisectBP3vd(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degreev, degreeu+1, 3, 3*(degreev+1), \
    (double*)ctlp, (double*)ctlq )
#define mbs_BisectBP4ud(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degreeu, 1, 4*(degreev+1), 0, \
    (double*)ctlp, (double*)ctlq )
#define mbs_BisectBP4vd(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesd ( degreev, degreeu+1, 4, 4*(degreev+1), \
    (double*)ctlp, (double*)ctlq )


void mbs_multiDivideBezCurvesd ( int degree, int ncurves,
                                 int spdimen, int pitch,
                                 double t,
                                 double *ctlp, double *ctlq );

#define mbs_DivideBC1d(degree,t,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degree, 1, 1, 0, t, ctlp, ctlq )
#define mbs_DivideBC2d(degree,t,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degree, 1, 2, 0, t, (double*)ctlp, (double*)ctlq )
#define mbs_DivideBC3d(degree,t,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degree, 1, 3, 0, t, (double*)ctlp, (double*)ctlq )
#define mbs_DivideBC4d(degree,t,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degree, 1, 4, 0, t, (double*)ctlp, (double*)ctlq )

#define mbs_DivideBP1ud(degreeu,degreev,u,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degreeu, 1, (degreev+1), 0, u, ctlp, ctlq )
#define mbs_DivideBP1vd(degreeu,degreev,v,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degreev, degreeu+1, 1, degreev+1, v, ctlp, ctlq )
#define mbs_DivideBP2ud(degreeu,degreev,u,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degreeu, 1, 2*(degreev+1), 0, u, \
    (double*)ctlp, (double*)ctlq )
#define mbs_DivideBP2vd(degreeu,degreev,v,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degreev, degreeu+1, 2, 2*(degreev+1), v, \
    (double*)ctlp, (double*)ctlq )
#define mbs_DivideBP3ud(degreeu,degreev,u,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degreeu, 1, 3*(degreev+1), 0, u, \
    (double*)ctlp, (double*)ctlq )
#define mbs_DivideBP3vd(degreeu,degreev,v,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degreev, degreeu+1, 3, 3*(degreev+1), v, \
    (double*)ctlp, (double*)ctlq )
#define mbs_DivideBP4ud(degreeu,degreev,u,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degreeu, 1, 4*(degreev+1), 0, u, \
    (double*)ctlp, (double*)ctlq )
#define mbs_DivideBP4vd(degreeu,degreev,v,ctlp,ctlq) \
  mbs_multiDivideBezCurvesd ( degreev, degreeu+1, 4, 4*(degreev+1), v, \
    (double*)ctlp, (double*)ctlq )

#ifndef MULTIBS_H
void mbs_BezP3NormalDeg ( int degreeu, int degreev, int *ndegu, int *ndegv );
void mbs_BezP3RNormalDeg ( int degreeu, int degreev, int *ndegu, int *ndegv );
#endif

char mbs_BezP3Normald ( int degreeu, int degreev, const point3d *ctlpoints,
                         int *ndegu, int *ndegv, vector3d *ncp );

char mbs_BezP3RNormald ( int degreeu, int degreev, const point4d *ctlpoints,
                          int *ndegu, int *ndegv, vector3d *ncp );

boolean mbs_ApproxBSKnotsValidd ( int degree, int lastknot, const double *knots,
                                  int lastiknot, const double *iknots );

int mbs_ApproxBSBandmSized ( int degree, const double *knots,
                             int lastiknot, const double *iknots );

boolean mbs_ConstructApproxBSProfiled ( int degree, int lastknot,
                                        const double *knots,
                                        int lastiknot, const double *iknots,
                                        bandm_profile *prof );

boolean mbs_ConstructApproxBSMatrixd ( int degree, int lastknot,
                                       const double *knots,
                                       int lastiknot, const double *iknots,
                                       int *nrows, int *ncols,
                                       bandm_profile *prof,
                                       double *a );


boolean mbs_multiConstructApproxBSCd ( int degree, int lastknot,
                                       const double *knots,
                                       int lastpknot, const double *pknots,
                                       int ncurves, int spdimen,
                                       int ppitch, const double *ppoints,
                                       int bcpitch, double *ctlpoints );

#define mbs_ConstructApproxBSC1d(degree,lastknot,knots,lastpknot,pknots,\
    ppoints,ctlpoints) \
  mbs_multiConstructApproxBSCd (degree,lastknot,knots,lastpknot,pknots,1,1,0,\
    (double*)ppoints,0,(double*)ctlpoints)
#define mbs_ConstructApproxBSC2d(degree,lastknot,knots,lastpknot,pknots,\
    ppoints,ctlpoints) \
  mbs_multiConstructApproxBSCd (degree,lastknot,knots,lastpknot,pknots,1,2,0,\
    (double*)ppoints,0,(double*)ctlpoints)
#define mbs_ConstructApproxBSC3d(degree,lastknot,knots,lastpknot,pknots,\
    ppoints,ctlpoints) \
  mbs_multiConstructApproxBSCd (degree,lastknot,knots,lastpknot,pknots,1,3,0,\
    (double*)ppoints,0,(double*)ctlpoints)
#define mbs_ConstructApproxBSC4d(degree,lastknot,knots,lastpknot,pknots,\
    ppoints,ctlpoints) \
  mbs_multiConstructApproxBSCd (degree,lastknot,knots,lastpknot,pknots,1,4,0,\
    (double*)ppoints,0,(double*)ctlpoints)


boolean mbs_OsloKnotsCorrectd ( int lastuknot, const double *uknots,
                                int lastvknot, const double *vknots );

int mbs_BuildOsloMatrixProfiled ( int degree,
                                  int lastuknot, const double *uknots,
                                  int lastvknot, const double *vknots,
                                  bandm_profile *prof );
void mbs_BuildOsloMatrixd ( int degree, int lastuknot, const double *uknots,
                            const double *vknots,
                            const bandm_profile *prof, double *a );


boolean mbs_multiOsloInsertKnotsd ( int ncurves, int spdimen, int degree,
                                    int inlastknot, const double *inknots,
                                    int inpitch, double *inctlpoints,
                                    int outlastknot, const double *outknots,
                                    int outpitch, double *outctlpoints );

boolean mbs_multiOsloRemoveKnotsLSQd ( int ncurves, int spdimen, int degree,
                                       int inlastknot, const double *inknots,
                                       int inpitch, double *inctlpoints,
                                       int outlastknot, const double *outknots,
                                       int outpitch, double *outctlpoints );


boolean mbs_multiBSChangeLeftKnotsd ( int ncurves, int spdimen, int degree,
                                      double *knots, int pitch, double *ctlpoints,
                                      double *newknots );

boolean mbs_multiBSChangeRightKnotsd ( int ncurves, int spdimen, int degree, 
                                       int lastknot, double *knots,  
                                       int pitch, double *ctlpoints, 
                                       double *newknots );

#define mbs_BSChangeLeftKnotsC1d(degree,knots,coeff,newknots) \
  mbs_multiBSChangeLeftKnotsd(1,1,degree,knots,0,coeff,newknots)
#define mbs_BSChangeLeftKnotsC2d(degree,knots,ctlpoints,newknots) \
  mbs_multiBSChangeLeftKnotsd(1,2,degree,knots,0,(double*)ctlpoints,newknots)
#define mbs_BSChangeLeftKnotsC3d(degree,knots,ctlpoints,newknots) \
  mbs_multiBSChangeLeftKnotsd(1,3,degree,knots,0,(double*)ctlpoints,newknots)
#define mbs_BSChangeLeftKnotsC4d(degree,knots,ctlpoints,newknots) \
  mbs_multiBSChangeLeftKnotsd(1,4,degree,knots,0,(double*)ctlpoints,newknots)
#define mbs_BSChangeRightKnotsC1d(degree,lastknot,knots,coeff,newknots) \
  mbs_multiBSChangeRightKnotsd(1,1,degree,lastknot,knots,0,coeff,newknots)
#define mbs_BSChangeRightKnotsC2d(degree,lastknot,knots,ctlpoints,newknots) \
  mbs_multiBSChangeRightKnotsd(1,2,degree,lastknot,knots,0,(double*)ctlpoints,newknots)
#define mbs_BSChangeRightKnotsC3d(degree,lastknot,knots,ctlpoints,newknots) \
  mbs_multiBSChangeRightKnotsd(1,3,degree,lastknot,knots,0,(double*)ctlpoints,newknots)
#define mbs_BSChangeRightKnotsC4d(degree,lastknot,knots,ctlpoints,newknots) \
  mbs_multiBSChangeRightKnotsd(1,4,degree,lastknot,knots,0,(double*)ctlpoints,newknots)


typedef struct {
    int     ident;
    boolean closing;
    byte    cpdimen;
    short   degree;
    int     lastknot;
    double  *knots;
    double  *points;
  } mbs_polycurved;

typedef struct {
    double t;             /* this must be the first field */
    char  sign1, sign2;
  } mbs_signpoint1d;

int  mbs_TrimCVBoundSized ( int nelem, const mbs_polycurved *bound );

void *mbs_CompileTrimPatchBoundd ( int nelem, const mbs_polycurved *bound,
                                   void *buffer );

void mbs_FindBoundLineIntersectionsd ( const void *bound,
                                       const point2d *p0, double t0,
                                       const point2d* p1, double t1,
                                       mbs_signpoint1d *inters, int *ninters );

boolean mbs_DrawTrimBSPatchDomd ( int degu, int lastuknot, const double *uknots,
                                  int degv, int lastvknot, const double *vknots,
                                  int nelem, const mbs_polycurved *bound,
                                  int nu, double au, double bu,
                                  int nv, double av, double bv,
                                  int maxinters,
                                  void *usrptr,
                                  void (*NotifyLine)(void*,char,int,point2d*,point2d*),
                                  void (*DrawLine)(void*,point2d*,point2d*,int),
                                  void (*DrawCurve)(void*,int,int,const double*) );


boolean mbs_MonotonicPolylined ( int spdimen, int npoints, int pitch,
                                 const double *points, const double *v );

boolean mbs_MonotonicPolylineRd ( int spdimen, int npoints, int pitch,
                                  const double *points, const double *v );


boolean mbs_RasterizeBC2d ( int degree, const point2d *cpoints,
                            void (*output)(const xpoint *buf, int n),
                            boolean outlast );

boolean mbs_RasterizeBC2Rd ( int degree, const point3d *cpoints,
                             void (*output)(const xpoint *buf, int n),
                             boolean outlast );

boolean mbs_RasterizeBS2d ( int degree, int lastknot, const double *knots,
                            const point2d *cpoints,
                            void (*output)(const xpoint *buf, int n),
                            boolean outlast );

boolean mbs_RasterizeBS2Rd ( int degree, int lastknot, const double *knots,
                             const point3d *cpoints,
                             void (*output)(const xpoint *buf, int n),
                             boolean outlast );


boolean mbs_multiInterp2knHermiteBezd ( int ncurves, int spdimen, int degree,
                                        int nlbc, int lbcpitch, const double *lbc,
                                        int nrbc, int rbcpitch, const double *rbc,
                                        int pitch, double *ctlpoints );

boolean mbs_multiInterp2knHermiteBSd ( int ncurves, int spdimen, int degree,
                                       int lastknot, const double *knots,
                                       int nlbc, int lbcpitch, const double *lbc, 
                                       int nrbc, int rbcpitch, const double *rbc, 
                                       int pitch, double *ctlpoints );


void mbs_multiFindBezDerivatived ( int degree, int ncurves, int spdimen,
                                   int pitch, const double *ctlpoints,  
                                   int dpitch, double *dctlpoints );

#define mbs_FindBezDerivativeC1d(degree,coeff,dcoeff) \
  mbs_multiFindBezDerivatived ( degree, 1, 1, 0, coeff, 0, dcoeff )
#define mbs_FindBezDerivativeC2d(degree,ctlpoints,dctlpoints) \
  mbs_multiFindBezDerivatived ( degree, 1, 2, 0, (double*)ctlpoints, \
    0, (double*)dctlpoints )
#define mbs_FindBezDerivativeC3d(degree,ctlpoints,dctlpoints) \
  mbs_multiFindBezDerivatived ( degree, 1, 3, 0, (double*)ctlpoints, \
    0, (double*)dctlpoints )
#define mbs_FindBezDerivativeC4d(degree,ctlpoints,dctlpoints) \
  mbs_multiFindBezDerivatived ( degree, 1, 4, 0, (double*)ctlpoints, \
    0, (double*)dctlpoints )


void mbs_multiFindBSDerivatived ( int degree, int lastknot, const double *knots, 
                                  int ncurves, int spdimen,  
                                  int pitch, const double *ctlpoints,  
                                  int *lastdknot, double *dknots,  
                                  int dpitch, double *dctlpoints );

#define mbs_FindBSDerivativeC1d(degree,lastknot,knots,coeff, \
    lastdknot,dknots,dcoeff) \
  mbs_multiFindBSDerivatived ( degree, lastknot, knots, 1, 1, 0, coeff, \
    lastdknot, dknots, 0, dcoeff )
#define mbs_FindBSDerivativeC2d(degree,lastknot,knots,ctlpoints, \
    lastdknot,dknots,dctlpoints) \
  mbs_multiFindBSDerivatived ( degree, lastknot, knots, 1, 2, 0, \
    (double*)ctlpoints, lastdknot, dknots, 0, (double*)dctlpoints )
#define mbs_FindBSDerivativeC3d(degree,lastknot,knots,ctlpoints, \
    lastdknot,dknots,dctlpoints) \
  mbs_multiFindBSDerivatived ( degree, lastknot, knots, 1, 3, 0, \
    (double*)ctlpoints, lastdknot, dknots, 0, (double*)dctlpoints )
#define mbs_FindBSDerivativeC4d(degree,lastknot,knots,ctlpoints, \
    lastdknot,dknots,dctlpoints) \
  mbs_multiFindBSDerivatived ( degree, lastknot, knots, 1, 4, 0, \
    (double*)ctlpoints, lastdknot, dknots, 0, (double*)dctlpoints )


boolean mbs_FindBSCommonKnotSequenced ( int *degree, int *lastknot,
            double **knots, int nsequences, ... );
boolean mbs_multiAdjustBSCRepd ( int ncurves, int spdimen,
     int indegree, int inlastknot, const double *inknots,
     int inpitch, const double *inctlpoints,
     int outdegree, int outlastknot, CONST_ double *outknots,
     int outpitch, double *outctlpoints );

#define mbs_AdjustBSCRepC1d(indegree,inlastknot,inknots, \
    inctlpoints,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiAdjustBSCRepd (1,1,indegree,inlastknot,inknots,0, \
    inctlpoints,outdegree,outlastknot,outknots,0,outctlpoints)
#define mbs_AdjustBSCRepC2d(indegree,inlastknot,inknots, \
    inctlpoints,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiAdjustBSCRepd (1,2,indegree,inlastknot,inknots,0, \
    (double*)inctlpoints,outdegree,outlastknot,outknots,0,(double*)outctlpoints)
#define mbs_AdjustBSCRepC3d(indegree,inlastknot,inknots, \
    inctlpoints,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiAdjustBSCRepd (1,3,indegree,inlastknot,inknots,0, \
    (double*)inctlpoints,outdegree,outlastknot,outknots,0,(double*)outctlpoints)
#define mbs_AdjustBSCRepC4d(indegree,inlastknot,inknots, \
    inctlpoints,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiAdjustBSCRepd (1,4,indegree,inlastknot,inknots,0, \
    (double*)inctlpoints,outdegree,outlastknot,outknots,0,(double*)outctlpoints)


boolean mbs_multiAddBSCurvesd ( int ncurves, int spdimen,
                             int degree1, int lastknot1, CONST_ double *knots1,
                             int pitch1, CONST_ double *ctlpoints1,
                             int degree2, int lastknot2, CONST_ double *knots2,
                             int pitch2, CONST_ double *ctlpoints2,
                             int *sumdeg, int *sumlastknot, double *sumknots,
                             int sumpitch, double *sumctlpoints );

#define mbs_AddBSCurvesC1d(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiAddBSCurvesd (1,1,degree1,lastknot1,knots1,0,ctlpoints1, \
    degree2,lastknot2,knots2,0,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,sumctlpoints)
#define mbs_AddBSCurvesC2d(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiAddBSCurvesd (1,2,degree1,lastknot1,knots1,0,(double*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(double*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(double*)sumctlpoints)
#define mbs_AddBSCurvesC3d(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiAddBSCurvesd (1,3,degree1,lastknot1,knots1,0,(double*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(double*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(double*)sumctlpoints)
#define mbs_AddBSCurvesC4d(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiAddBSCurvesd (1,4,degree1,lastknot1,knots1,0,(double*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(double*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(double*)sumctlpoints)


boolean mbs_multiSubtractBSCurvesd ( int ncurves, int spdimen,
                             int degree1, int lastknot1, CONST_ double *knots1,
                             int pitch1, CONST_ double *ctlpoints1,
                             int degree2, int lastknot2, CONST_ double *knots2,
                             int pitch2, CONST_ double *ctlpoints2,
                             int *sumdeg, int *sumlastknot, double *sumknots,
                             int sumpitch, double *sumctlpoints );

#define mbs_SubtractBSCurvesC1d(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiSubtractBSCurvesd (1,1,degree1,lastknot1,knots1,0,ctlpoints1, \
    degree2,lastknot2,knots2,0,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,sumctlpoints)
#define mbs_SubtractBSCurvesC2d(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiSubtractBSCurvesd (1,2,degree1,lastknot1,knots1,0,(double*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(double*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(double*)sumctlpoints)
#define mbs_SubtractBSCurvesC3d(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiSubtractBSCurvesd (1,3,degree1,lastknot1,knots1,0,(double*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(double*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(double*)sumctlpoints)
#define mbs_SubtractBSCurvesC4d(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiSubtractBSCurvesd (1,4,degree1,lastknot1,knots1,0,(double*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(double*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(double*)sumctlpoints)


boolean mbs_FindPolynomialZerosd ( int degree, const double *coeff,
                                   int *nzeros, double *zeros, double eps );


boolean mbs_ClipBC2d ( int ncplanes, const vector3d *cplanes,
                       int degree, const point2d *cpoints,
                       boolean (*output) (int degree, const point2d *cpoints) );               
boolean mbs_ClipBC2Rd ( int ncplanes, const vector3d *cplanes,
                        int degree, const point3d *cpoints,   
                        boolean (*output) (int degree, const point3d *cpoints) );
boolean mbs_ClipBC3d ( int ncplanes, const vector4d *cplanes,
                       int degree, const point3d *cpoints,   
                       boolean (*output) (int degree, const point3d *cpoints) );
boolean mbs_ClipBC3Rd ( int ncplanes, const vector4d *cplanes,
                        int degree, const point4d *cpoints,
                        boolean (*output) (int degree, const point4d *cpoints) );

/* ///////////////////////////////////////////////////////////////////////// */
/* Bicubic polynomial Coons patches */
boolean _mbs_BezC1CoonsFindCornersd ( int spdimen,
                                     int degc00, const double *c00,
                                     int degc01, const double *c01,
                                     int degc10, const double *c10,
                                     int degc11, const double *c11,
                                     double *pcorners,
                                     double *workspace );
boolean mbs_BezC1CoonsFindCornersd ( int spdimen,
                                     int degc00, const double *c00,
                                     int degc01, const double *c01,
                                     int degc10, const double *c10,
                                     int degc11, const double *c11,
                                     double *pcorners );
boolean _mbs_BezC1CoonsToBezd ( int spdimen,
                                int degc00, const double *c00,
                                int degc01, const double *c01,
                                int degc10, const double *c10,
                                int degc11, const double *c11,
                                int degd00, const double *d00,
                                int degd01, const double *d01,
                                int degd10, const double *d10,
                                int degd11, const double *d11,
                                int *n, int *m, double *p,
                                double *workspace );
boolean mbs_BezC1CoonsToBezd ( int spdimen,
                               int degc00, const double *c00,
                               int degc01, const double *c01,
                               int degc10, const double *c10,
                               int degc11, const double *c11,
                               int degd00, const double *d00,
                               int degd01, const double *d01,
                               int degd10, const double *d10,
                               int degd11, const double *d11,
                               int *n, int *m, double *p );

boolean mbs_TabCubicHFuncDer2d ( double a, double b, int nkn, const double *kn,
                                 double *hfunc, double *dhfunc, double *ddhfunc );
boolean mbs_TabCubicHFuncDer3d ( double a, double b, int nkn, const double *kn, 
                                 double *hfunc, double *dhfunc, double *ddhfunc,
                                 double *dddhfunc );
boolean _mbs_TabBezCurveDer2d ( int spdimen, int degree, const double *cp,
                                int nkn, const double *kn,
                                int ppitch,
                                double *p, double *dp, double *ddp,
                                double *workspace );
boolean mbs_TabBezCurveDer2d ( int spdimen, int degree, const double *cp,
                               int nkn, const double *kn,
                               int ppitch,
                               double *p, double *dp, double *ddp );
void _mbs_TabBezC1Coonsd ( int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp,
                   double *workspace );
void _mbs_TabBezC1Coonsd ( int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp,
                   double *workspace );
boolean _mbs_TabBezC1CoonsDer2d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd10, const double *d10,
      int degd11, const double *d11,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *workspace );
boolean mbs_TabBezC1CoonsDer2d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd10, const double *d10,
      int degd11, const double *d11,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv );
boolean _mbs_TabBezC1CoonsDer3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd10, const double *d10,
      int degd11, const double *d11,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv,
      double *workspace );
boolean mbs_TabBezC1CoonsDer3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd10, const double *d10,
      int degd11, const double *d11,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv );

void _mbs_TabBezC1Coons0d ( int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp,
                   double *workspace );
boolean _mbs_TabBezC1Coons0Der2d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degd00, const double *d00,
      int degd01, const double *d01,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *workspace );
boolean mbs_TabBezC1Coons0Der2d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degd00, const double *d00,
      int degd01, const double *d01,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv );
boolean _mbs_TabBezC1Coons0Der3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degd00, const double *d00,
      int degd01, const double *d01,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv,
      double *workspace );
boolean mbs_TabBezC1Coons0Der3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degd00, const double *d00,
      int degd01, const double *d01,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv );

/* ///////////////////////////////////////////////////////////////////////// */
/* Biquintic polynomial Coons patches */
boolean _mbs_BezC2CoonsFindCornersd ( int spdimen,
                                      int degc00, const double *c00,
                                      int degc01, const double *c01,
                                      int degc02, const double *c02,
                                      int degc10, const double *c10,
                                      int degc11, const double *c11,
                                      int degc12, const double *c12,
                                      double *pcorners,
                                      double *workspace );
boolean mbs_BezC2CoonsFindCornersd ( int spdimen,
                                     int degc00, const double *c00,
                                     int degc01, const double *c01,
                                     int degc02, const double *c02,
                                     int degc10, const double *c10,
                                     int degc11, const double *c11,
                                     int degc12, const double *c12,
                                     double *pcorners );
boolean _mbs_BezC2CoonsToBezd ( int spdimen,
                                int degc00, const double *c00,
                                int degc01, const double *c01,
                                int degc02, const double *c02,
                                int degc10, const double *c10,
                                int degc11, const double *c11,
                                int degc12, const double *c12,
                                int degd00, const double *d00,
                                int degd01, const double *d01,
                                int degd02, const double *d02,
                                int degd10, const double *d10,
                                int degd11, const double *d11,
                                int degd12, const double *d12,
                                int *n, int *m, double *p,
                                double *workspace );
boolean mbs_BezC2CoonsToBezd ( int spdimen,
                               int degc00, const double *c00,
                               int degc01, const double *c01,
                               int degc02, const double *c02,
                               int degc10, const double *c10,
                               int degc11, const double *c11,
                               int degc12, const double *c12,
                               int degd00, const double *d00,
                               int degd01, const double *d01,
                               int degd02, const double *d02,
                               int degd10, const double *d10,
                               int degd11, const double *d11,
                               int degd12, const double *d12,
                               int *n, int *m, double *p );
boolean mbs_TabQuinticHFuncDer3d ( double a, double b, int nkn, const double *kn,
                                   double *hfunc, double *dhfunc,
                                   double *ddhfunc, double *dddhfunc );
boolean _mbs_TabBezCurveDer3d ( int spdimen, int degree, const double *cp,
                                int nkn, const double *kn,
                                int ppitch,
                                double *p, double *dp, double *ddp, double *dddp,
                                double *workspace );
boolean mbs_TabBezCurveDer3d ( int spdimen, int degree, const double *cp,
                               int nkn, const double *kn,
                               int ppitch,
                               double *p, double *dp, double *ddp, double *dddp );
void _mbs_TabBezC2Coonsd ( int spdimen, int nknu, int nknv,
                           const double *c, const double *d, const double *p,
                           const double *hu, const double *hv, double *pp,
                           double *workspace );

boolean _mbs_TabBezC2CoonsDer3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc02, const double *c02,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degc12, const double *c12,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd02, const double *d02,
      int degd10, const double *d10,
      int degd11, const double *d11,
      int degd12, const double *d12,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv,
      double *workspace );
boolean mbs_TabBezC2CoonsDer3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc02, const double *c02,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degc12, const double *c12,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd02, const double *d02,
      int degd10, const double *d10,
      int degd11, const double *d11,
      int degd12, const double *d12,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv );

void _mbs_TabBezC2Coons0d ( int spdimen, int nknu, int nknv,
                            const double *c, const double *d, const double *p,
                            const double *hu, const double *hv, double *pp,
                            double *workspace );
void _mbs_TabBezC2Coons0d ( int spdimen, int nknu, int nknv,
                            const double *c, const double *d, const double *p,
                            const double *hu, const double *hv, double *pp,
                            double *workspace );
boolean _mbs_TabBezC2Coons0Der3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc02, const double *c02,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd02, const double *d02,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv,
      double *workspace );
boolean mbs_TabBezC2Coons0Der3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc02, const double *c02,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd02, const double *d02,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv );

/* ///////////////////////////////////////////////////////////////////////// */
/* Bicubic B-spline Coons patches */          
boolean mbs_BSC1CoonsFindCornersd ( int spdimen,
            int degc00, int lastknotc00, const double *knotsc00, const double *c00,
            int degc01, int lastknotc01, const double *knotsc01, const double *c01,
            int degc10, int lastknotc10, const double *knotsc10, const double *c10,
            int degc11, int lastknotc11, const double *knotsc11, const double *c11,
            double *pcorners );
boolean mbs_BSC1CoonsToBSd ( int spdimen,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degc10, int lastknotc10, const double *knotsc10, const double *c10,
      int degc11, int lastknotc11, const double *knotsc11, const double *c11,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      int degd10, int lastknotd10, const double *knotsd10, const double *d10,
      int degd11, int lastknotd11, const double *knotsd11, const double *d11,
      int *degreeu, int *lastuknot, double *uknots,
      int *degreev, int *lastvknot, double *vknots, double *p );

void mbs_TabBSCurveDer2d ( int spdimen, int degree, int lastknot,
                           const double *knots, const double *cp,
                           int nkn, const double *kn, int ppitch,
                           double *p, double *dp, double *ddp );
boolean _mbs_TabBSC1Coonsd ( int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp );
boolean mbs_TabBSC1CoonsDer2d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degc10, int lastknotc10, const double *knotsc10, const double *c10,
      int degc11, int lastknotc11, const double *knotsc11, const double *c11,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      int degd10, int lastknotd10, const double *knotsd10, const double *d10,
      int degd11, int lastknotd11, const double *knotsd11, const double *d11,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv );
boolean mbs_TabBSC1CoonsDer3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degc10, int lastknotc10, const double *knotsc10, const double *c10,
      int degc11, int lastknotc11, const double *knotsc11, const double *c11,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      int degd10, int lastknotd10, const double *knotsd10, const double *d10,
      int degd11, int lastknotd11, const double *knotsd11, const double *d11,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,   
      double *puuu, double *puuv, double *puvv, double *pvvv );

boolean _mbs_TabBSC1Coons0d ( int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp );
boolean mbs_TabBSC1Coons0Der2d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv );
boolean mbs_TabBSC1Coons0Der3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,   
      double *puuu, double *puuv, double *puvv, double *pvvv );

/* ///////////////////////////////////////////////////////////////////////// */
/* Biquintic B-spline Coons patches */        
boolean mbs_BSC2CoonsFindCornersd ( int spdimen,
            int degc00, int lastknotc00, const double *knotsc00, const double *c00,
            int degc01, int lastknotc01, const double *knotsc01, const double *c01,
            int degc02, int lastknotc02, const double *knotsc02, const double *c02,
            int degc10, int lastknotc10, const double *knotsc10, const double *c10,
            int degc11, int lastknotc11, const double *knotsc11, const double *c11,
            int degc12, int lastknotc12, const double *knotsc12, const double *c12,
            double *pcorners );
boolean mbs_BSC2CoonsToBSd ( int spdimen,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degc02, int lastknotc02, const double *knotsc02, const double *c02,
      int degc10, int lastknotc10, const double *knotsc10, const double *c10,
      int degc11, int lastknotc11, const double *knotsc11, const double *c11,
      int degc12, int lastknotc12, const double *knotsc12, const double *c12,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      int degd02, int lastknotd02, const double *knotsd02, const double *d02,
      int degd10, int lastknotd10, const double *knotsd10, const double *d10,
      int degd11, int lastknotd11, const double *knotsd11, const double *d11,
      int degd12, int lastknotd12, const double *knotsd12, const double *d12,
      int *degreeu, int *lastuknot, double *uknots,
      int *degreev, int *lastvknot, double *vknots, double *p );
void mbs_TabBSCurveDer3d ( int spdimen, int degree, int lastknot,
                           const double *knots, const double *cp,
                           int nkn, const double *kn, int ppitch,
                           double *p, double *dp, double *ddp, double *dddp );
boolean _mbs_TabBSC2Coonsd ( int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp );
boolean mbs_TabBSC2CoonsDer3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degc02, int lastknotc02, const double *knotsc02, const double *c02,
      int degc10, int lastknotc10, const double *knotsc10, const double *c10,
      int degc11, int lastknotc11, const double *knotsc11, const double *c11,
      int degc12, int lastknotc12, const double *knotsc12, const double *c12,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      int degd02, int lastknotd02, const double *knotsd02, const double *d02,
      int degd10, int lastknotd10, const double *knotsd10, const double *d10,
      int degd11, int lastknotd11, const double *knotsd11, const double *d11,
      int degd12, int lastknotd12, const double *knotsd12, const double *d12,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv );
boolean _mbs_TabBSC2Coons0d ( int spdimen, int nknu, int nknv,
                   const double *c, const double *d, const double *p,
                   const double *hu, const double *hv, double *pp );
boolean mbs_TabBSC2Coons0Der3d ( int spdimen,
      int nknu, const double *knu, const double *hfuncu,
      const double *dhfuncu, const double *ddhfuncu, const double *dddhfuncu,
      int nknv, const double *knv, const double *hfuncv,
      const double *dhfuncv, const double *ddhfuncv, const double *dddhfuncv,
      int degc00, int lastknotc00, const double *knotsc00, const double *c00,
      int degc01, int lastknotc01, const double *knotsc01, const double *c01,
      int degc02, int lastknotc02, const double *knotsc02, const double *c02,
      int degd00, int lastknotd00, const double *knotsd00, const double *d00,
      int degd01, int lastknotd01, const double *knotsd01, const double *d01,
      int degd02, int lastknotd02, const double *knotsd02, const double *d02,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv,
      double *puuu, double *puuv, double *puvv, double *pvvv );


/* ///////////////////////////////////////////////////////////////////////// */
/* spherical product of two curves */
void mbs_SphericalProductd (
                int degree_eq, int lastknot_eq, const point2d *cpoints_eq,
                int degree_mer, int lastknot_mer, const point2d *cpoints_mer,
                int pitch, point3d *spr_cp );
void mbs_SphericalProductRd (
                int degree_eq, int lastknot_eq, const point3d *cpoints_eq,
                int degree_mer, int lastknot_mer, const point3d *cpoints_mer,
                int pitch, point4d *spr_cp );

/* ///////////////////////////////////////////////////////////////////////// */
/* Lane-Riesenfeld algorithm */
boolean mbs_multiLaneRiesenfeldd ( int spdimen, int ncurves, int degree,
                         int inlastknot, int inpitch, const double *incp,
                         int *outlastknot, int outpitch, double *outcp );

#define mbs_LaneRiesenfeldC1d(degree,inlastknot,incp,outlastknot,outcp) \
  mbs_multiLaneRiesenfeldd ( 1, 1, degree, inlastknot, 0, incp, \
    outlastknot, 0, outcp )
#define mbs_LaneRiesenfeldC2d(degree,inlastknot,incp,outlastknot,outcp) \
  mbs_multiLaneRiesenfeldd ( 2, 1, degree, inlastknot, 0, (double*)incp, \
    outlastknot, 0, (double*)outcp )
#define mbs_LaneRiesenfeldC3d(degree,inlastknot,incp,outlastknot,outcp) \
  mbs_multiLaneRiesenfeldd ( 3, 1, degree, inlastknot, 0, (double*)incp, \
    outlastknot, 0, (double*)outcp )
#define mbs_LaneRiesenfeldC4d(degree,inlastknot,incp,outlastknot,outcp) \
  mbs_multiLaneRiesenfeldd ( 4, 1, degree, inlastknot, 0, (double*)incp, \
    outlastknot, 0, (double*)outcp )


#ifdef __cplusplus
}
#endif

#endif /*MULTIBSD_H*/

