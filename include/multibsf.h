
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libmultibs library of C procedures -              */
/* processing B-spline and Bezier curves and surfaces                    */ 

#ifndef CONST_  /* a dirty trick to suppress many compiler warning messages */
#define CONST_ const
#endif

#ifndef MULTIBSF_H
#define MULTIBSF_H

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


int mbs_KnotMultiplicityf ( int lastknot, const float *knots, float t );

int mbs_FindKnotIntervalf ( int degree, int lastknot, const float *knots,
                            float t, int *mult );

float mbs_GrevilleAbscissaf ( int degree, float *knots, int i );

void mbs_TransformAffKnotsf ( int degree, int lastknot, const float *inknots,
                              float a, float b, float *outknots );

void mbs_multiReverseBSCurvef ( int degree, int lastknot, float *knots,
                                int ncurves, int spdimen, int pitch,
                                float *ctlpoints );
boolean mbs_ClosedKnotsCorrectf ( int degree, int lastknot, float *knots,
                                  float T, int K, float tol );

int mbs_SetKnotf ( int lastknot, float *knots,
                   int knotnum, int mult, float t );
int mbs_SetKnotClosedf ( int degree, int lastknot, float *knots, float T,
                         int knotnum, int mult, float t );

void _mbs_multideBoorKernelf ( int degree, const float *knots,
                               int ncurves, int spdimen,
                               int pitch, const float *ctlpoints,
                               float t, int k, int r, int lj,
                               int dpitch, float *d );

int mbs_multideBoorf ( int degree, int lastknot,
                       const float *knots,
                       int ncurves, int spdimen, int pitch,
                       const float *ctlpoints,
                       float t,
                       float *cpoints );

#define mbs_deBoorC1f(degree,lastknot,knots,coeff,t,value) \
  mbs_multideBoorf ( degree, lastknot, knots, 1, 1, 0, coeff, t, value )
#define mbs_deBoorC2f(degree,lastknot,knots,ctlpoints,t,cpoint) \
  mbs_multideBoorf ( degree, lastknot, knots, 1, 2, 0, (float*)ctlpoints, t, \
    (float*)cpoint )
#define mbs_deBoorC3f(degree,lastknot,knots,ctlpoints,t,cpoint) \
  mbs_multideBoorf ( degree, lastknot, knots, 1, 3, 0, (float*)ctlpoints, t, \
    (float*)cpoint )
#define mbs_deBoorC4f(degree,lastknot,knots,ctlpoints,t,cpoint) \
  mbs_multideBoorf ( degree, lastknot, knots, 1, 4, 0, (float*)ctlpoints, t, \
    (float*)cpoint )


void mbs_deBoorC2Rf ( int degree, int lastknot,
                      const float *knots, const point3f *ctlpoints,
                      float t, point2f *cpoint );

void mbs_deBoorC3Rf ( int degree, int lastknot,
                      const float *knots, const point4f *ctlpoints,
                      float t, point3f *cpoint );

void mbs_deBoorP3f ( int degreeu, int lastknotu, const float *knotsu,
                     int degreev, int lastknotv, const float *knotsv,
                     int pitch,
                     const point3f *ctlpoints,
                     float u, float v, point3f *ppoint );

void mbs_deBoorP3Rf ( int degreeu, int lastknotu, const float *knotsu,
                      int degreev, int lastknotv, const float *knotsv,
                      int pitch,
                      const point4f *ctlpoints,
                      float u, float v, point3f *ppoint );

void mbs_deBoorP4f ( int degreeu, int lastknotu, const float *knotsu,
                     int degreev, int lastknotv, const float *knotsv,
                     int pitch,
                     const point4f *ctlpoints,
                     float u, float v, point4f *ppoint );


int mbs_multideBoorDerf ( int degree, int lastknot,
                          const float *knots,
                          int ncurves, int spdimen, int pitch,
                          const float *ctlpoints,
                          float t,
                          float *cpoints, float *dervect );

#define mbs_deBoorDerC1f(degree,lastknot,knots,ctlpoints,t,cpoint,cder) \
  mbs_multideBoorDerf ( degree, lastknot, knots, 1, 1, 0, (float*)ctlpoints, \
    t, cpoint, cder )
#define mbs_deBoorDerC2f(degree,lastknot,knots,ctlpoints,t,cpoint,cder) \
  mbs_multideBoorDerf ( degree, lastknot, knots, 1, 2, 0, (float*)ctlpoints, \
    t, (float*)cpoint, (float*)cder )
#define mbs_deBoorDerC3f(degree,lastknot,knots,ctlpoints,t,cpoint,cder) \
  mbs_multideBoorDerf ( degree, lastknot, knots, 1, 3, 0, (float*)ctlpoints, \
    t, (float*)cpoint, (float*)cder )
#define mbs_deBoorDerC4f(degree,lastknot,knots,ctlpoints,t,cpoint,cder) \
  mbs_multideBoorDerf ( degree, lastknot, knots, 1, 4, 0, (float*)ctlpoints, \
    t, (float*)cpoint, (float*)cder )


int mbs_multideBoorDer2f ( int degree, int lastknot, const float *knots,
                           int ncurves, int spdimen,
                           int pitch, const float *ctlpoints,
                           float t, float *p, float *d1, float *d2 );

#define mbs_deBoorDer2C1f(degree,lastknot,knots,coeff,t,p,d1,d2) \
  mbs_multideBoorDer2f(degree,lastknot,knots,1,1,0,coeff,t,p,d1,d2)
#define mbs_deBoorDer2C2f(degree,lastknot,knots,ctlpoints,t,p,d1,d2) \
  mbs_multideBoorDer2f(degree,lastknot,knots,1,2,0,(float*)ctlpoints,t, \
    (float*)p,(float*)d1,(float*)d2)
#define mbs_deBoorDer2C3f(degree,lastknot,knots,ctlpoints,t,p,d1,d2) \
  mbs_multideBoorDer2f(degree,lastknot,knots,1,3,0,(float*)ctlpoints,t, \
    (float*)p,(float*)d1,(float*)d2)
#define mbs_deBoorDer2C4f(degree,lastknot,knots,ctlpoints,t,p,d1,d2) \
  mbs_multideBoorDer2f(degree,lastknot,knots,1,4,0,(float*)ctlpoints,t, \
    (float*)p,(float*)d1,(float*)d2)


int mbs_multideBoorDer3f ( int degree, int lastknot, const float *knots,
                           int ncurves, int spdimen,
                           int pitch, const float *ctlpoints, float t,
                           float *p, float *d1, float *d2, float *d3 );

#define mbs_deBoorDer3C1f(degree,lastknot,knots,coeff,t,p,d1,d2,d3) \
  mbs_multideBoorDer3f(degree,lastknot,knots,1,1,0,coeff,t,p,d1,d2,d3)
#define mbs_deBoorDer3C2f(degree,lastknot,knots,ctlpoints,t,p,d1,d2,d3) \
  mbs_multideBoorDer3f(degree,lastknot,knots,1,2,0,(float*)ctlpoints,t, \
    (float*)p,(float*)d1,(float*)d2,(float*)d3)
#define mbs_deBoorDer3C3f(degree,lastknot,knots,ctlpoints,t,p,d1,d2,d3) \
  mbs_multideBoorDer3f(degree,lastknot,knots,1,3,0,(float*)ctlpoints,t, \
    (float*)p,(float*)d1,(float*)d2,(float*)d3)
#define mbs_deBoorDer3C4f(degree,lastknot,knots,ctlpoints,t,p,d1,d2,d3) \
  mbs_multideBoorDer3f(degree,lastknot,knots,1,4,0,(float*)ctlpoints,t, \
    (float*)p,(float*)d1,(float*)d2,(float*)d3)


boolean mbs_deBoorDerPf ( int degreeu, int lastknotu, const float *knotsu,    
                          int degreev, int lastknotv, const float *knotsv,
                          int spdimen, int pitch, const float *ctlpoints,
                          float u, float v,
                          float *ppoint, float *uder, float *vder );

boolean mbs_deBoorDer2Pf ( int degreeu, int lastknotu, const float *knotsu,    
                           int degreev, int lastknotv, const float *knotsv,
                           int spdimen, int pitch, const float *ctlpoints,
                           float u, float v,
                           float *ppoint, float *uder, float *vder,
                           float *uuder, float *uvder, float *vvder );

boolean mbs_deBoorDer3Pf ( int degreeu, int lastknotu, const float *knotsu,    
                           int degreev, int lastknotv, const float *knotsv,
                           int spdimen, int pitch, const float *ctlpoints,
                           float u, float v,
                           float *ppoint, float *uder, float *vder,
                           float *uuder, float *uvder, float *vvder,
                           float *uuuder, float *uuvder, float *uvvder,
                           float *vvvder );


int mbs_multiKnotInsf ( int degree, int *lastknot,
                        float *knots,
                        int ncurves, int spdimen, int inpitch, int outpitch,
                        float *ctlpoints, float t );

#define mbs_KnotInsC1f(degree,lastknot,knots,coeff,t) \
  mbs_multiKnotInsf(degree,lastknot,knots,1,1,0,0,coeff, t)
#define mbs_KnotInsC2f(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsf(degree,lastknot,knots,1,2,0,0,(float*)ctlpoints, t)
#define mbs_KnotInsC3f(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsf(degree,lastknot,knots,1,3,0,0,(float*)ctlpoints, t)
#define mbs_KnotInsC4f(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsf(degree,lastknot,knots,1,4,0,0,(float*)ctlpoints, t)


int mbs_multiKnotInsClosedf ( int degree, int *lastknot, float *knots,
                              int ncurves, int spdimen, int inpitch, int outpitch,
                              float *ctlpoints, float t );

#define mbs_KnotInsClosedC1f(degree,lastknot,knots,coeff,t) \
  mbs_multiKnotInsClosedf(degree,lastknot,knots,1,1,0,0,coeff,t)
#define mbs_KnotInsClosedC2f(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsClosedf(degree,lastknot,knots,1,2,0,0,(float*)ctlpoints,t)
#define mbs_KnotInsClosedC3f(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsClosedf(degree,lastknot,knots,1,3,0,0,(float*)ctlpoints,t)
#define mbs_KnotInsClosedC4f(degree,lastknot,knots,ctlpoints,t) \
  mbs_multiKnotInsClosedf(degree,lastknot,knots,1,4,0,0,(float*)ctlpoints,t)


int mbs_multiKnotRemovef ( int degree, int *lastknot,
                           float *knots,
                           int ncurves, int spdimen, int inpitch, int outpitch,
                           float *ctlpoints, int knotnum );

#define mbs_KnotRemoveC1f(degree,lastknot,knots,coeff,knotnum) \
  mbs_multiKnotRemovef(degree,lastknot,knots,1,1,0,0,coeff,knotnum)
#define mbs_KnotRemoveC2f(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemovef(degree,lastknot,knots,1,2,0,0, \
    (float*)ctlpoints, knotnum)
#define mbs_KnotRemoveC3f(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemovef(degree,lastknot,knots,1,3,0,0, \
    (float*)ctlpoints, knotnum)
#define mbs_KnotRemoveC4f(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemovef(degree,lastknot,knots,1,4,0,0, \
    (float*)ctlpoints, knotnum)


int mbs_multiKnotRemoveClosedf ( int degree, int *lastknot,
                                 float *knots,
                                 int ncurves, int spdimen, int inpitch, int outpitch,
                                 float *ctlpoints, int knotnum );

#define mbs_KnotRemoveClosedC1f(degree,lastknot,knots,coeff,knotnum) \
  mbs_multiKnotRemoveClosedf (degree,lastknot,knots,1,1,0,0,coeff,knotnum)
#define mbs_KnotRemoveClosedC2f(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemoveClosedf (degree,lastknot,knots,1,2,0,0, \
    (float*)ctlpoints,knotnum)
#define mbs_KnotRemoveClosedC3f(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemoveClosedf (degree,lastknot,knots,1,3,0,0, \
    (float*)ctlpoints,knotnum)
#define mbs_KnotRemoveClosedC4f(degree,lastknot,knots,ctlpoints,knotnum) \
  mbs_multiKnotRemoveClosedf (degree,lastknot,knots,1,4,0,0, \
    (float*)ctlpoints,knotnum)


void mbs_multiRemoveSuperfluousKnotsf ( int ncurves, int spdimen, int degree,
                                        int *lastknot, float *knots,
                                        int inpitch, int outpitch,
                                        float *ctlpoints );


void mbs_multiMaxKnotInsf ( int ncurves, int spdimen, int degree,
                            int inlastknot, const float *inknots,
                            int inpitch, const float *inctlpoints,
                            int *outlastknot, float *outknots,
                            int outpitch, float *outctlpoints,
                            int *skipl, int *skipr );

#define mbs_MaxKnotInsC1f(degree,inlastknot,inknots,incoeff, \
                          outlastknot,outknots,outcoeff,skipl,skipr) \
  mbs_multiMaxKnotInsf(1,1,degree,inlastknot,inknots,0,incoeff, \
    outlastknot,outknots,0,outcoeff,skipl,skipr)
#define mbs_MaxKnotInsC2f(degree,inlastknot,inknots,inctlpoints, \
                          outlastknot,outknots,outctlpoints,skipl,skipr) \
  mbs_multiMaxKnotInsf(1,2,degree,inlastknot,inknots,0,(float*)inctlpoints, \
    outlastknot,outknots,0,(float*)outctlpoints,skipl,skipr)
#define mbs_MaxKnotInsC3f(degree,inlastknot,inknots,inctlpoints, \
                          outlastknot,outknots,outctlpoints,skipl,skipr) \
  mbs_multiMaxKnotInsf(1,3,degree,inlastknot,inknots,0,(float*)inctlpoints, \
    outlastknot,outknots,0,(float*)outctlpoints,skipl,skipr)
#define mbs_MaxKnotInsC4f(degree,inlastknot,inknots,inctlpoints, \
                          outlastknot,outknots,outctlpoints,skipl,skipr) \
  mbs_multiMaxKnotInsf(1,4,degree,inlastknot,inknots,0,(float*)inctlpoints, \
    outlastknot,outknots,0,(float*)outctlpoints,skipl,skipr)


void mbs_multiBSCurvesToBezf ( int spdimen, int ncurves,
                               int degree, int lastinknot, const float *inknots,
                               int inpitch, const float *inctlp,
                               int *kpcs, int *lastoutknot, float *outknots,
                               int outpitch, float *outctlp );

#define mbs_BSToBezC1f(degree,lastinknot,inknots,incoeff,kpcs, \
    lastoutknot,outknots,outcoeff) \
  mbs_multiBSCurvesToBezf(1,1,degree,lastinknot,inknots,0,incoeff, \
    kpcs,lastoutknot,outknots,0,outcoeff)
#define mbs_BSToBezC2f(degree,lastinknot,inknots,inctlp,kpcs, \
    lastoutknot,outknots,outctlp) \
  mbs_multiBSCurvesToBezf(2,1,degree,lastinknot,inknots,0,(float*)inctlp, \
    kpcs,lastoutknot,outknots,0,(float*)outctlp)
#define mbs_BSToBezC3f(degree,lastinknot,inknots,inctlp,kpcs, \
    lastoutknot,outknots,outctlp) \
  mbs_multiBSCurvesToBezf(3,1,degree,lastinknot,inknots,0,(float*)inctlp, \
    kpcs,lastoutknot,outknots,0,(float*)outctlp)
#define mbs_BSToBezC4f(degree,lastinknot,inknots,inctlp,kpcs, \
    lastoutknot,outknots,outctlp) \
  mbs_multiBSCurvesToBezf(4,1,degree,lastinknot,inknots,0,(float*)inctlp, \
    kpcs,lastoutknot,outknots,0,(float*)outctlp)


void mbs_BSPatchToBezf ( int spdimen,
                         int degreeu, int lastuknot, const float *uknots,
                         int degreev, int lastvknot, const float *vknots,
                         int inpitch, const float *inctlp,
                         int *kupcs, int *lastoutuknot, float *outuknots,
                         int *kvpcs, int *lastoutvknot, float *outvknots,
                         int outpitch, float *outctlp );


int mbs_NumKnotIntervalsf ( int degree, int lastknot, const float *knots );

int mbs_LastknotMaxInsf ( int degree, int lastknot, const float *knots, 
                          int *numknotintervals );

int mbs_BSProdRepSizef ( int degree1, int lastknot1, const float *knots1,
                         int degree2, int lastknot2, const float *knots2 );

int mbs_NumMaxKnotsf ( int degree, int lastknot, const float *knots );

void mbs_SetBSProdKnotsf ( int degree1, int lastknot1, const float *knots1,
                           int degree2, int lastknot2, const float *knots2,
                           int *degree, int *lastknot, float *knots );

void mbs_SetKnotPatternf ( int lastinknot, const float *inknots,
                           int multipl,
                           int *lastoutknot, float *outknots );

void mbs_multiBezScalef ( int degree, int narcs, int ncurves, int spdimen,  
                          int pitch, float *ctlpoints );

void mbs_multiBezUnscalef ( int degree, int narcs, int ncurves, int spdimen,  
                            int pitch, float *ctlpoints );

void mbs_multiMultBezCf ( int nscf, int degscf, int scfpitch,
                          const float *scfcoef,
                          int spdimen,
                          int nvecf, int degvecf, int vecfpitch,
                          const float *vecfcp,
                          int *degprod, int prodpitch, float *prodcp );

void mbs_multiMultBSCf ( int nscf, int degscf,
                         int scflastknot, const float *scfknots,
                         int scfpitch, const float *scfcoef,
                         int spdimen,
                         int nvecf, int degvecf,
                         int vecflastknot, const float *vecfknots,
                         int vecfpitch, const float *vecfcp,
                         int *degprod, int *prodlastknot, float *prodknots,
                         int prodpitch, float *prodcp );


void mbs_multiBCDegElevf ( int ncurves, int spdimen,
                           int inpitch, int indegree, const float *inctlpoints,
                           int deltadeg,
                           int outpitch, int *outdegree, float *outctlpoints );

#define mbs_BCDegElevC1f(indegree,incoeff,deltadeg,outdegree,outcoeff) \
  mbs_multiBCDegElevf ( 1, 1, 0, indegree, incoeff, deltadeg, \
    0, outdegree, outcoeff )
#define mbs_BCDegElevC2f(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegElevf ( 1, 2, 0, indegree, (float*)inctlpoints, deltadeg, \
    0, outdegree, (float*)outctlpoints )
#define mbs_BCDegElevC3f(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegElevf ( 1, 3, 0, indegree, (float*)inctlpoints, deltadeg, \
    0, outdegree, (float*)outctlpoints)
#define mbs_BCDegElevC4f(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegElevf ( 1, 4, 0, indegree, (float*)inctlpoints, deltadeg, \
    0, outdegree, (float*)outctlpoints )


void mbs_BCDegElevPf ( int spdimen,
                       int indegreeu, int indegreev, const float *inctlp,
                       int deltadegu, int deltadegv,
                       int *outdegreeu, int *outdegreev,
                       float *outctlp );

#define mbs_BCDegElevP1f(indegreeu,indegreev,incoeff,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outcoeff) \
  mbs_BCDegElevPf ( 1, indegreeu, indegreev, incoeff, deltadegu, deltadegv, \
    outdegreeu, outdegreev, outcoeff )
#define mbs_BCDegElevP2f(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegElevPf ( 2, indegreeu, indegreev, (float*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (float*)outctlp )
#define mbs_BCDegElevP3f(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegElevPf ( 3, indegreeu, indegreev, (float*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (float*)outctlp )
#define mbs_BCDegElevP4f(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegElevPf ( 4, indegreeu, indegreev, (float*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (float*)outctlp )


void mbs_multiBCDegRedf ( int ncurves, int spdimen,
                          int inpitch, int indegree, const float *inctlpoints,
                          int deltadeg,
                          int outpitch, int *outdegree, float *outctlpoints );

#define mbs_BCDegRedC1f(indegree,incoeff,deltadeg,outdegree,outcoeff) \
  mbs_multiBCDegRedf ( 1, 1, 0, indegree, incoeff, deltadeg, \
    0, outdegree, outcoeff )
#define mbs_BCDegRedC2f(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegRedf ( 1, 2, 0, indegree, (float*)inctlpoints, deltadeg, \
    0, outdegree, (float*)outctlpoints )
#define mbs_BCDegRedC3f(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegRedf ( 1, 3, 0, indegree, (float*)inctlpoints, deltadeg, \
    0, outdegree, (float*)outctlpoints)
#define mbs_BCDegRedC4f(indegree,inctlpoints,deltadeg,outdegree,outctlpoints) \
  mbs_multiBCDegRedf ( 1, 4, 0, indegree, (float*)inctlpoints, deltadeg, \
    0, outdegree, (float*)outctlpoints )


void mbs_BCDegRedPf ( int spdimen,
                      int indegreeu, int indegreev, const float *inctlp,
                      int deltadegu, int deltadegv,
                      int *outdegreeu, int *outdegreev,
                      float *outctlp );

#define mbs_BCDegRedP1f(indegreeu,indegreev,incoeff,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outcoeff) \
  mbs_BCDegRedPf ( 1, indegreeu, indegreev, incoeff, deltadegu, deltadegv, \
    outdegreeu, outdegreev, outcoeff )
#define mbs_BCDegRedP2f(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegRedPf ( 2, indegreeu, indegreev, (float*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (float*)outctlp )
#define mbs_BCDegRedP3f(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegRedPf ( 3, indegreeu, indegreev, (float*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (float*)outctlp )
#define mbs_BCDegRedP4f(indegreeu,indegreev,inctlp,deltadegu,deltadegv, \
    outdegreeu,outdegreev,outctlp) \
  mbs_BCDegRedPf ( 4, indegreeu, indegreev, (float*)inctlp, \
    deltadegu, deltadegv, outdegreeu, outdegreev, (float*)outctlp )


void mbs_multiBSDegElevf ( int ncurves, int spdimen,
                           int indegree, int inlastknot, const float *inknots,
                           int inpitch, const float *inctlpoints,
                           int deltadeg,
                           int *outdegree, int *outlastknot,
                           float *outknots, int outpitch, float *outctlpoints,
                           boolean freeend );

#define mbs_BSDegElevC1f(indegree,inlastknot,inknots,incoeff, \
    deltadeg,outdegree,outlastknot,outknots,outcoeff,freeend) \
  mbs_multiBSDegElevf(1,1,indegree,inlastknot,inknots,0,incoeff, \
    deltadeg,outdegree,outlastknot,outknots,0,outcoeff,freeend)
#define mbs_BSDegElevC2f(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints,freeend) \
  mbs_multiBSDegElevf(1,2,indegree,inlastknot,inknots,0,(float*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(float*)outctlpoints,freeend)
#define mbs_BSDegElevC3f(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints,freeend) \
  mbs_multiBSDegElevf(1,3,indegree,inlastknot,inknots,0,(float*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(float*)outctlpoints,freeend)
#define mbs_BSDegElevC4f(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints,freeend) \
  mbs_multiBSDegElevf(1,4,indegree,inlastknot,inknots,0,(float*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(float*)outctlpoints,freeend)

void mbs_multiBSDegElevClosedf ( int ncurves, int spdimen,
                         int indegree, int inlastknot, const float *inknots,
                         int inpitch, const float *inctlpoints,
                         int deltadeg,
                         int *outdegree, int *outlastknot,
                         float *outknots, int outpitch, float *outctlpoints );

#define mbs_BSDegElevClosedC1f(indegree,inlastknot,inknots,incoeff, \
    deltadeg,outdegree,outlastknot,outknots,outcoeff) \
  mbs_multiBSDegElevClosedf(1,1,indegree,inlastknot,inknots,0,incoeff, \
    deltadeg,outdegree,outlastknot,outknots,0,outcoeff)
#define mbs_BSDegElevClosedC2f(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiBSDegElevClosedf(1,2,indegree,inlastknot,inknots,0,(float*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(float*)outctlpoints)
#define mbs_BSDegElevClosedC3f(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiBSDegElevClosedf(1,3,indegree,inlastknot,inknots,0,(float*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(float*)outctlpoints)
#define mbs_BSDegElevClosedC4f(indegree,inlastknot,inknots,inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiBSDegElevClosedf(1,4,indegree,inlastknot,inknots,0,(float*)inctlpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(float*)outctlpoints)


boolean mbs_multiBSDegRedf ( int ncurves, int spdimen,
                             int indegree, int inlastknot, const float *inknots,
                             int inpitch, CONST_ float *inctlpoints,
                             int deltadeg,
                             int *outdegree, int *outlastknot, float *outknots,
                             int outpitch, float *outctlpoints );

#define mbs_BSDegRedC1f(indegree,inlastknot,inknots,incoeff,deltadeg, \
    outdegree,outlastknot,outknots,outcoeff) \
  mbs_multiBSDegRedf(1,1,indegree,inlastknot,inknots,0,incoeff,deltadeg, \
    outdegree,outlastknot,outknots,0,outcoeff)
#define mbs_BSDegRedC2f(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedf(1,2,indegree,inlastknot,inknots,0,(float*)incpoints,deltadeg, \
    outdegree,outlastknot,outknots,0,(float*)outcpoints)
#define mbs_BSDegRedC3f(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedf(1,3,indegree,inlastknot,inknots,0,(float*)incpoints,deltadeg, \
    outdegree,outlastknot,outknots,0,(float*)outcpoints)
#define mbs_BSDegRedC4f(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedf(1,4,indegree,inlastknot,inknots,0,(float*)incpoints,deltadeg, \
    outdegree,outlastknot,outknots,0,(float*)outcpoints)


boolean mbs_multiBSDegRedClosedf ( int ncurves, int spdimen,
                             int indegree, int inlastknot, const float *inknots,
                             int inpitch, CONST_ float *inctlpoints,
                             int deltadeg,
                             int *outdegree, int *outlastknot, float *outknots,
                             int outpitch, float *outctlpoints );

#define mbs_BSDegRedClosedC1f(indegree,inlastknot,inknots,incoeff,deltadeg, \
    outdegree,outlastknot,outknots,outcoeff) \
  mbs_multiBSDegRedClosedf(1,1,indegree,inlastknot,inknots,0,incoeff,deltadeg, \
    outdegree,outlastknot,outknots,0,outcoeff)
#define mbs_BSDegRedClosedC2f(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedClosedf(1,2,indegree,inlastknot,inknots,0,(float*)incpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(float*)outcpoints)
#define mbs_BSDegRedClosedC3f(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedClosedf(1,3,indegree,inlastknot,inknots,0,(float*)incpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(float*)outcpoints)
#define mbs_BSDegRedClosedC4f(indegree,inlastknot,inknots,incpoints,deltadeg, \
    outdegree,outlastknot,outknots,outcpoints) \
  mbs_multiBSDegRedClosedf(1,4,indegree,inlastknot,inknots,0,(float*)incpoints, \
    deltadeg,outdegree,outlastknot,outknots,0,(float*)outcpoints)


void mbs_multiBCHornerf ( int degree, int ncurves, int spdimen, int pitch,
                          const float *ctlpoints, float t, float *cpoints );

#define mbs_BCHornerC1f(degree,coeff,t,value) \
  mbs_multiBCHornerf ( degree, 1, 1, 0, coeff, t, value )
#define mbs_BCHornerC2f(degree,ctlpoints,t,cpoint) \
  mbs_multiBCHornerf ( degree, 1, 2, 0, (float*)ctlpoints, t, (float*)cpoint )
#define mbs_BCHornerC3f(degree,ctlpoints,t,cpoint) \
  mbs_multiBCHornerf ( degree, 1, 3, 0, (float*)ctlpoints, t, (float*)cpoint )
#define mbs_BCHornerC4f(degree,ctlpoints,t,cpoint) \
  mbs_multiBCHornerf ( degree, 1, 4, 0, (float*)ctlpoints, t, (float*)cpoint )

void mbs_BCHornerC2Rf ( int degree, const point3f *ctlpoints, float t,
                        point2f *cpoint );

void mbs_BCHornerC3Rf ( int degree, const point4f *ctlpoints, float t,
                        point3f *cpoint );


void mbs_BCHornerPf ( int degreeu, int degreev, int spdimen,
                      const float *ctlpoints,
                      float u, float v, float *ppoint );

#define mbs_BCHornerP1f(degreeu,degreev,coeff,u,v,ppoint) \
  mbs_BCHornerPf ( degreeu, degreev, 1, coeff, u, v, ppoint )
#define mbs_BCHornerP2f(degreeu,degreev,ctlpoints,u,v,ppoint ) \
  mbs_BCHornerPf ( degreeu, degreev, 2, (float*)ctlpoints, \
    u, v, (float*)ppoint )
#define mbs_BCHornerP3f(degreeu,degreev,ctlpoints,u,v,ppoint ) \
  mbs_BCHornerPf ( degreeu, degreev, 3, (float*)ctlpoints, \
    u, v, (float*)ppoint )
#define mbs_BCHornerP4f(degreeu,degreev,ctlpoints,u,v,ppoint ) \
  mbs_BCHornerPf ( degreeu, degreev, 4, (float*)ctlpoints, \
    u, v, (float*)ppoint )

void mbs_BCHornerP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                        float u, float v, point3f *p );


void mbs_multiBCHornerDerf ( int degree, int ncurves, int spdimen, int pitch,
                             const float *ctlpoints,
                             float t, float *p, float *d );

#define mbs_BCHornerDerC1f(degree,coeff,t,p,d) \
  mbs_multiBCHornerDerf ( degree, 1, 1, 0, coeff, t, p, d )
#define mbs_BCHornerDerC2f(degree,ctlpoints,t,p,d) \
  mbs_multiBCHornerDerf ( degree, 1, 2, 0, (float*)ctlpoints, t, \
    (float*)p, (float*)d )
#define mbs_BCHornerDerC3f(degree,ctlpoints,t,p,d) \
  mbs_multiBCHornerDerf ( degree, 1, 3, 0, (float*)ctlpoints, t, \
    (float*)p, (float*)d )
#define mbs_BCHornerDerC4f(degree,ctlpoints,t,p,d) \
  mbs_multiBCHornerDerf ( degree, 1, 4, 0, (float*)ctlpoints, t, \
    (float*)p, (float*)d )

void mbs_BCHornerDerC2Rf ( int degree, const point3f *ctlpoints, float t,
                           point2f *p, vector2f *d );
void mbs_BCHornerDerC3Rf ( int degree, const point4f *ctlpoints, float t,
                           point3f *p, vector3f *d );


void mbs_BCHornerDerPf ( int degreeu, int degreev, int spdimen,
                         const float *ctlpoints,
                         float u, float v,  
                         float *p, float *du, float *dv );

#define mbs_BCHornerDerP1f(degreeu,degreev,coeff,u,v,p,du,dv) \
  mbs_BCHornerDerPf ( degreeu, degreev, 1, coeff, u, v, p, du, dv )
#define mbs_BCHornerDerP2f(degreeu,degreev,ctlpoints,u,v,p,du,dv) \
  mbs_BCHornerDerPf ( degreeu, degreev, 2, (float*)ctlpoints, u, v, \
    (float*)p, (float*)du, (float*)dv )
#define mbs_BCHornerDerP3f(degreeu,degreev,ctlpoints,u,v,p,du,dv) \
  mbs_BCHornerDerPf ( degreeu, degreev, 3, (float*)ctlpoints, u, v, \
    (float*)p, (float*)du, (float*)dv )
#define mbs_BCHornerDerP4f(degreeu,degreev,ctlpoints,u,v,p,du,dv) \
  mbs_BCHornerDerPf ( degreeu, degreev, 4, (float*)ctlpoints, u, v, \
    (float*)p, (float*)du, (float*)dv )

void mbs_BCHornerDerP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                           float u, float v,
                           point3f *p, vector3f *du, vector3f *dv );


void mbs_BCHornerNvP3f ( int degreeu, int degreev, const point3f *ctlpoints,
                         float u, float v,
                         point3f *p, vector3f *nv );

void mbs_BCHornerNvP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                          float u, float v,
                          point3f *p, vector3f *nv );


void mbs_multiBCHornerDer2f ( int degree, int ncurves, int spdimen, int pitch,
                              const float *ctlpoints, 
                              float t, float *p, float *d1, float *d2 );

#define mbs_BCHornerDer2C1f(degree,coeff,t,p,d1,d2) \
  mbs_multiBCHornerDer2f ( degree, 1, 1, 0, coeff, t, p, d1, d2 )
#define mbs_BCHornerDer2C2f(degree,ctlpoints,t,p,d1,d2) \
  mbs_multiBCHornerDer2f ( degree, 1, 2, 0, (float*)ctlpoints, t, \
    (float*)p, (float*)d1, (float*)d2 )
#define mbs_BCHornerDer2C3f(degree,ctlpoints,t,p,d1,d2) \
  mbs_multiBCHornerDer2f ( degree, 1, 3, 0, (float*)ctlpoints, t, \
    (float*)p, (float*)d1, (float*)d2 )
#define mbs_BCHornerDer2C4f(degree,ctlpoints,t,p,d1,d2) \
  mbs_multiBCHornerDer2f ( degree, 1, 4, 0, (float*)ctlpoints, t, \
    (float*)p, (float*)d1, (float*)d2 )

void mbs_BCHornerDer2C2Rf ( int degree, const point3f *ctlpoints, float t,
                            point2f *p, vector2f *d1, vector2f *d2 );
void mbs_BCHornerDer2C3Rf ( int degree, const point4f *ctlpoints, float t,
                            point3f *p, vector3f *d1, vector3f *d2 );

void mbs_BCHornerDer2Pf ( int degreeu, int degreev, int spdimen,
                          const float *ctlpoints,
                          float u, float v,
                          float *p, float *du, float *dv,
                          float *duu, float *duv, float *dvv );

#define mbs_BCHornerDer2P1f(degreeu,degreev,coeff,u,v,p,du,dv,duu,duv,dvv) \
  mbs_BCHornerDer2Pf ( degreeu, degreev, 1, coeff, u, v, \
    p, du, dv, duu, duv, dvv )
#define mbs_BCHornerDer2P2f(degreeu,degreev,ctlpoints,u,v,p,du,dv,duu,duv,dvv) \
  mbs_BCHornerDer2Pf ( degreeu, degreev, 2, (float*)ctlpoints, u, v, \
    (float*)p, (float*)du, (float*)dv, (float*)duu, (float*)duv, (float*)dvv )
#define mbs_BCHornerDer2P3f(degreeu,degreev,ctlpoints,u,v,p,du,dv,duu,duv,dvv) \
  mbs_BCHornerDer2Pf ( degreeu, degreev, 3, (float*)ctlpoints, u, v, \
    (float*)p, (float*)du, (float*)dv, (float*)duu, (float*)duv, (float*)dvv )
#define mbs_BCHornerDer2P4f(degreeu,degreev,ctlpoints,u,v,p,du,dv,duu,duv,dvv) \
  mbs_BCHornerDer2Pf ( degreeu, degreev, 4, (float*)ctlpoints, u, v, \
    (float*)p, (float*)du, (float*)dv, (float*)duu, (float*)duv, (float*)dvv )

void mbs_BCHornerDer2P3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                            float u, float v,
                            point3f *p, vector3f *du, vector3f *dv,
                            vector3f *duu, vector3f *duv, vector3f *dvv );


void mbs_multiBCHornerDer3f ( int degree, int ncurves, int spdimen, int pitch,
                              const float *ctlpoints, float t,
                              float *p, float *d1, float *d2, float *d3 );

#define mbs_BCHornerDer3C1f(degree,coeff,t,p,d1,d2,d3) \
  mbs_multiBCHornerDer3f ( degree, 1, 1, 0, coeff, t, p, d1, d2, d3 )
#define mbs_BCHornerDer3C2f(degree,ctlpoints,t,p,d1,d2,d3) \
  mbs_multiBCHornerDer3f ( degree, 1, 2, 0, (float*)ctlpoints, t, \
    (float*)p, (float*)d1, (float*)d2, (float*)d3 )
#define mbs_BCHornerDer3C3f(degree,ctlpoints,t,p,d1,d2,d3) \
  mbs_multiBCHornerDer3f ( degree, 1, 3, 0, (float*)ctlpoints, t, \
    (float*)p, (float*)d1, (float*)d2, (float*)d3 )
#define mbs_BCHornerDer3C4f(degree,ctlpoints,t,p,d1,d2,d3) \
  mbs_multiBCHornerDer3f ( degree, 1, 4, 0, (float*)ctlpoints, t, \
    (float*)p, (float*)d1, (float*)d2, (float*)d3 )


void mbs_FindBezPatchDiagFormf ( int degreeu, int degreev, int spdimen,
                                 CONST_ float *cpoints,
                                 int k, int l, float u, float v,
                                 float *dfcp );

void mbs_BCHornerDer3Pf ( int degreeu, int degreev, int spdimen,
                          CONST_ float *ctlpoints,
                          float u, float v,
                          float *p, float *pu, float *pv,
                          float *puu, float *puv, float *pvv,
                          float *puuu, float *puuv, float *puvv, float *pvvv );

#define mbs_BCHornerDer3P1f(degreeu,degreev,coeff,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv) \
  mbs_BCHornerDer3Pf ( degreeu, degreev, 1, coeff, u, v, \
    p, pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv )
#define mbs_BCHornerDer3P2f(degreeu,degreev,ctlpoints,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv) \
  mbs_BCHornerDer3Pf ( degreeu, degreev, 2, (float*)ctlpoints, u, v, \
    (float*)p, (float*)pu, (float*)pv, (float*)puu, (float*)puv, \
    (float*)pvv, (float*)puuu, (float*)puuv, (float*)puvv, (float*)pvvv )
#define mbs_BCHornerDer3P3f(degreeu,degreev,ctlpoints,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv) \
  mbs_BCHornerDer3Pf ( degreeu, degreev, 3, (float*)ctlpoints, u, v, \
    (float*)p, (float*)pu, (float*)pv, (float*)puu, (float*)puv, \
    (float*)pvv, (float*)puuu, (float*)puuv, (float*)puvv, (float*)pvvv )
#define mbs_BCHornerDer3P4f(degreeu,degreev,ctlpoints,u,v, \
    p,pu,pv,puu,puv,pvv,puuu,puuv,puvv,pvvv) \
  mbs_BCHornerDer3Pf ( degreeu, degreev, 4, (float*)ctlpoints, u, v, \
    (float*)p, (float*)pu, (float*)pv, (float*)puu, (float*)puv, \
    (float*)pvv, (float*)puuu, (float*)puuv, (float*)puvv, (float*)pvvv )


void mbs_deBoorBasisf ( int degree, int lastknot, const float *knots,
                        float t, int *fnz, int *nnz, float *bfv );


void mbs_multiBSCubicInterpf ( int lastinterpknot, float *interpknots,
                               int ncurves, int spdimen, int xpitch,
                               const float *x,
                               int ypitch,
                               char bcl, const float *ybcl,
                               char bcr, const float *ybcr,
                               int *lastbsknot, float *bsknots,
                               int bspitch, float *ctlpoints );

void mbs_multiBSCubicClosedInterpf ( int lastinterpknot, float *interpknots,
                               int ncurves, int spdimen, int xpitch,
                               const float *x,
                               int *lastbsknot, float *bsknots,
                               int bspitch, float *ctlpoints );
                                                                                                                            

void mbs_BCFrenetC2f ( int degree, const point2f *ctlpoints, float t,
                       point2f *cpoint, vector2f *fframe, float *curvature );

void mbs_BCFrenetC2Rf ( int degree, const point3f *ctlpoints, float t,
                        point2f *cpoint, vector2f *fframe, float *curvature );

void mbs_BCFrenetC3f ( int degree, const point3f *ctlpoints, float t,
                       point3f *cpoint, vector3f *fframe, float *curvatures );

void mbs_BCFrenetC3Rf ( int degree, const point4f *ctlpoints, float t,
                        point3f *cpoint, vector3f *fframe, float *curvatures );


void mbs_FundFormsBP3f ( int degreeu, int degreev, const point3f *ctlpoints,
                         float u, float v,
                         float *firstform, float *secondform );

void mbs_GMCurvaturesBP3f ( int degreeu, int degreev, const point3f *ctlpoints,
                            float u, float v,
                            float *gaussian, float *mean );

void mbs_PrincipalDirectionsBP3f ( int degreeu, int degreev,
                                   const point3f *ctlpoints,
                                   float u, float v,
                                   float *k1, vector2f *v1,
                                   float *k2, vector2f *v2 );

void mbs_FundFormsBP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                          float u, float v,
                          float *firstform, float *secondform );

void mbs_GMCurvaturesBP3Rf ( int degreeu, int degreev, const point4f *ctlpoints,
                             float u, float v,
                             float *gaussian, float *mean );

void mbs_PrincipalDirectionsBP3Rf ( int degreeu, int degreev,
                                    const point4f *ctlpoints,
                                    float u, float v,
                                    float *k1, vector2f *v1,
                                    float *k2, vector2f *v2 );


void mbs_multiBisectBezCurvesf ( int degree, int ncurves,
                                 int spdimen, int pitch,
                                 float *ctlp, float *ctlq );

#define mbs_BisectBC1f(degree,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degree, 1, 1, 0, ctlp, ctlq )
#define mbs_BisectBC2f(degree,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degree, 1, 2, 0, (float*)ctlp, (float*)ctlq )
#define mbs_BisectBC3f(degree,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degree, 1, 3, 0, (float*)ctlp, (float*)ctlq )
#define mbs_BisectBC4f(degree,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degree, 1, 4, 0, (float*)ctlp, (float*)ctlq )

#define mbs_BisectBP1uf(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degreeu, 1, (degreev+1), 0, ctlp, ctlq )
#define mbs_BisectBP1vf(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degreev, degreeu+1, 1, degreev+1, ctlp, ctlq )
#define mbs_BisectBP2uf(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degreeu, 1, 2*(degreev+1), 0, \
    (float*)ctlp, (float*)ctlq )
#define mbs_BisectBP2vf(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degreev, degreeu+1, 2, 2*(degreev+1), \
    (float*)ctlp, (float*)ctlq )
#define mbs_BisectBP3uf(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degreeu, 1, 3*(degreev+1), 0, \
    (float*)ctlp, (float*)ctlq )
#define mbs_BisectBP3vf(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degreev, degreeu+1, 3, 3*(degreev+1), \
    (float*)ctlp, (float*)ctlq )
#define mbs_BisectBP4uf(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degreeu, 1, 4*(degreev+1), 0, \
    (float*)ctlp, (float*)ctlq )
#define mbs_BisectBP4vf(degreeu,degreev,ctlp,ctlq) \
  mbs_multiBisectBezCurvesf ( degreev, degreeu+1, 4, 4*(degreev+1), \
    (float*)ctlp, (float*)ctlq )


void mbs_multiDivideBezCurvesf ( int degree, int ncurves,
                                 int spdimen, int pitch,
                                 float t,
                                 float *ctlp, float *ctlq );

#define mbs_DivideBC1f(degree,t,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degree, 1, 1, 0, t, ctlp, ctlq )
#define mbs_DivideBC2f(degree,t,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degree, 1, 2, 0, t, (float*)ctlp, (float*)ctlq )
#define mbs_DivideBC3f(degree,t,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degree, 1, 3, 0, t, (float*)ctlp, (float*)ctlq )
#define mbs_DivideBC4f(degree,t,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degree, 1, 4, 0, t, (float*)ctlp, (float*)ctlq )

#define mbs_DivideBP1uf(degreeu,degreev,u,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degreeu, 1, (degreev+1), 0, u, ctlp, ctlq )
#define mbs_DivideBP1vf(degreeu,degreev,v,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degreev, degreeu+1, 1, degreev+1, v, ctlp, ctlq )
#define mbs_DivideBP2uf(degreeu,degreev,u,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degreeu, 1, 2*(degreev+1), 0, u, \
    (float*)ctlp, (float*)ctlq )
#define mbs_DivideBP2vf(degreeu,degreev,v,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degreev, degreeu+1, 2, 2*(degreev+1), v, \
    (float*)ctlp, (float*)ctlq )
#define mbs_DivideBP3uf(degreeu,degreev,u,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degreeu, 1, 3*(degreev+1), 0, u, \
    (float*)ctlp, (float*)ctlq )
#define mbs_DivideBP3vf(degreeu,degreev,v,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degreev, degreeu+1, 3, 3*(degreev+1), v, \
    (float*)ctlp, (float*)ctlq )
#define mbs_DivideBP4uf(degreeu,degreev,u,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degreeu, 1, 4*(degreev+1), 0, u, \
    (float*)ctlp, (float*)ctlq )
#define mbs_DivideBP4vf(degreeu,degreev,v,ctlp,ctlq) \
  mbs_multiDivideBezCurvesf ( degreev, degreeu+1, 4, 4*(degreev+1), v, \
    (float*)ctlp, (float*)ctlq )

#ifndef MULTIBS_H
void mbs_BezP3NormalDeg ( int degreeu, int degreev, int *ndegu, int *ndegv );
void mbs_BezP3RNormalDeg ( int degreeu, int degreev, int *ndegu, int *ndegv );
#endif

char mbs_BezP3Normalf ( int degreeu, int degreev, const point3f *ctlpoints,
                         int *ndegu, int *ndegv, vector3f *ncp );

char mbs_BezP3RNormalf ( int degreeu, int degreev, const point4f *ctlpoints,
                          int *ndegu, int *ndegv, vector3f *ncp );

boolean mbs_ApproxBSKnotsValidf ( int degree, int lastknot, const float *knots,
                                  int lastiknot, const float *iknots );

int mbs_ApproxBSBandmSizef ( int degree, const float *knots,
                             int lastiknot, const float *iknots );

boolean mbs_ConstructApproxBSProfilef ( int degree, int lastknot,
                                        const float *knots,
                                        int lastiknot, const float *iknots,
                                        bandm_profile *prof );

boolean mbs_ConstructApproxBSMatrixf ( int degree, int lastknot,
                                       const float *knots,
                                       int lastiknot, const float *iknots,
                                       int *nrows, int *ncols,
                                       bandm_profile *prof,
                                       float *a );


boolean mbs_multiConstructApproxBSCf ( int degree, int lastknot,
                                       const float *knots,
                                       int lastpknot, const float *pknots,
                                       int ncurves, int spdimen,
                                       int ppitch, const float *ppoints,
                                       int bcpitch, float *ctlpoints );

#define mbs_ConstructApproxBSC1f(degree,lastknot,knots,lastpknot,pknots,\
    ppoints,ctlpoints) \
  mbs_multiConstructApproxBSCf (degree,lastknot,knots,lastpknot,pknots,1,1,0,\
    (float*)ppoints,0,(float*)ctlpoints)
#define mbs_ConstructApproxBSC2f(degree,lastknot,knots,lastpknot,pknots,\
    ppoints,ctlpoints) \
  mbs_multiConstructApproxBSCf (degree,lastknot,knots,lastpknot,pknots,1,2,0,\
    (float*)ppoints,0,(float*)ctlpoints)
#define mbs_ConstructApproxBSC3f(degree,lastknot,knots,lastpknot,pknots,\
    ppoints,ctlpoints) \
  mbs_multiConstructApproxBSCf (degree,lastknot,knots,lastpknot,pknots,1,3,0,\
    (float*)ppoints,0,(float*)ctlpoints)
#define mbs_ConstructApproxBSC4f(degree,lastknot,knots,lastpknot,pknots,\
    ppoints,ctlpoints) \
  mbs_multiConstructApproxBSCf (degree,lastknot,knots,lastpknot,pknots,1,4,0,\
    (float*)ppoints,0,(float*)ctlpoints)


boolean mbs_OsloKnotsCorrectf ( int lastuknot, const float *uknots,
                                int lastvknot, const float *vknots );

int mbs_BuildOsloMatrixProfilef ( int degree,
                                  int lastuknot, const float *uknots,
                                  int lastvknot, const float *vknots,
                                  bandm_profile *prof );
void mbs_BuildOsloMatrixf ( int degree, int lastuknot, const float *uknots,
                            const float *vknots,
                            const bandm_profile *prof, float *a );


void mbs_multiOsloInsertKnotsf ( int ncurves, int spdimen, int degree,
                                 int inlastknot, const float *inknots,
                                 int inpitch, float *inctlpoints,
                                 int outlastknot, const float *outknots,
                                 int outpitch, float *outctlpoints );

void mbs_multiOsloRemoveKnotsLSQf ( int ncurves, int spdimen, int degree,
                                    int inlastknot, const float *inknots,
                                    int inpitch, float *inctlpoints,
                                    int outlastknot, const float *outknots,
                                    int outpitch, float *outctlpoints );


void mbs_multiBSChangeLeftKnotsf ( int ncurves, int spdimen, int degree,
                                   float *knots, int pitch, float *ctlpoints,
                                   float *newknots );

void mbs_multiBSChangeRightKnotsf ( int ncurves, int spdimen, int degree, 
                                    int lastknot, float *knots,  
                                    int pitch, float *ctlpoints, 
                                    float *newknots );

#define mbs_BSChangeLeftKnotsC1f(degree,knots,coeff,newknots) \
  mbs_multiBSChangeLeftKnotsf(1,1,degree,knots,0,coeff,newknots)
#define mbs_BSChangeLeftKnotsC2f(degree,knots,ctlpoints,newknots) \
  mbs_multiBSChangeLeftKnotsf(1,2,degree,knots,0,(float*)ctlpoints,newknots)
#define mbs_BSChangeLeftKnotsC3f(degree,knots,ctlpoints,newknots) \
  mbs_multiBSChangeLeftKnotsf(1,3,degree,knots,0,(float*)ctlpoints,newknots)
#define mbs_BSChangeLeftKnotsC4f(degree,knots,ctlpoints,newknots) \
  mbs_multiBSChangeLeftKnotsf(1,4,degree,knots,0,(float*)ctlpoints,newknots)
#define mbs_BSChangeRightKnotsC1f(degree,lastknot,knots,coeff,newknots) \
  mbs_multiBSChangeRightKnotsf(1,1,degree,lastknot,knots,0,coeff,newknots)
#define mbs_BSChangeRightKnotsC2f(degree,lastknot,knots,ctlpoints,newknots) \
  mbs_multiBSChangeRightKnotsf(1,2,degree,lastknot,knots,0,(float*)ctlpoints,newknots)
#define mbs_BSChangeRightKnotsC3f(degree,lastknot,knots,ctlpoints,newknots) \
  mbs_multiBSChangeRightKnotsf(1,3,degree,lastknot,knots,0,(float*)ctlpoints,newknots)
#define mbs_BSChangeRightKnotsC4f(degree,lastknot,knots,ctlpoints,newknots) \
  mbs_multiBSChangeRightKnotsf(1,4,degree,lastknot,knots,0,(float*)ctlpoints,newknots)


typedef struct {
    boolean closing;
    byte    spdimen;
    short   degree;
    short   lastknot;
    float   *knots;
    float   *points;
  } polycurvef;

typedef struct {
    float t;             /* this must be the first field */
    char  sign1, sign2;
  } signpoint1f;

int  mbs_TrimCVBoundSizef ( int nelem, const polycurvef *bound );

void *mbs_CompileTrimPatchBoundf ( int nelem, const polycurvef *bound,
                                   void *buffer );

void mbs_FindBoundLineIntersectionsf ( const void *bound,
                                       const point2f *p0, float t0,
                                       const point2f* p1, float t1,
                                       signpoint1f *inters, int *ninters );

void mbs_DrawTrimBSPatchDomf ( int degu, int lastuknot, const float *uknots,
                               int degv, int lastvknot, const float *vknots,
                               int nelem, const polycurvef *bound,
                               int nu, float au, float bu,
                               int nv, float av, float bv,
                               int maxinters,
                               void (*NotifyLine)(char,int,point2f*,point2f*),
                               void (*DrawLine)(point2f*,point2f*,int),
                               void (*DrawCurve)(int,int,const float*) );


boolean mbs_MonotonicPolylinef ( int spdimen, int npoints, int pitch,
                                 const float *points, const float *v );

boolean mbs_MonotonicPolylineRf ( int spdimen, int npoints, int pitch,
                                  const float *points, const float *v );


void mbs_RasterizeBC2f ( int degree, const point2f *cpoints,
                         void (*output)(const xpoint *buf, int n),
                         boolean outlast );

void mbs_RasterizeBC2Rf ( int degree, const point3f *cpoints,
                          void (*output)(const xpoint *buf, int n),
                          boolean outlast );

void mbs_RasterizeBS2f ( int degree, int lastknot, const float *knots,
                         const point2f *cpoints,
                         void (*output)(const xpoint *buf, int n),
                         boolean outlast );

void mbs_RasterizeBS2Rf ( int degree, int lastknot, const float *knots,
                          const point3f *cpoints,
                          void (*output)(const xpoint *buf, int n),
                          boolean outlast );


void mbs_multiInterp2knHermiteBezf ( int ncurves, int spdimen, int degree,
                                     int nlbc, int lbcpitch, const float *lbc, 
                                     int nrbc, int rbcpitch, const float *rbc, 
                                     int pitch, float *ctlpoints );

void mbs_multiInterp2knHermiteBSf ( int ncurves, int spdimen, int degree,
                                    int lastknot, const float *knots,
                                    int nlbc, int lbcpitch, const float *lbc, 
                                    int nrbc, int rbcpitch, const float *rbc, 
                                    int pitch, float *ctlpoints );


void mbs_multiFindBezDerivativef ( int degree, int ncurves, int spdimen,  
                                   int pitch, const float *ctlpoints,  
                                   int dpitch, float *dctlpoints );

#define mbs_FindBezDerivativeC1f(degree,coeff,dcoeff) \
  mbs_multiFindBezDerivativef ( degree, 1, 1, 0, coeff, 0, dcoeff )
#define mbs_FindBezDerivativeC2f(degree,ctlpoints,dctlpoints) \
  mbs_multiFindBezDerivativef ( degree, 1, 2, 0, (float*)ctlpoints, \
    0, (float*)dctlpoints )
#define mbs_FindBezDerivativeC3f(degree,ctlpoints,dctlpoints) \
  mbs_multiFindBezDerivativef ( degree, 1, 3, 0, (float*)ctlpoints, \
    0, (float*)dctlpoints )
#define mbs_FindBezDerivativeC4f(degree,ctlpoints,dctlpoints) \
  mbs_multiFindBezDerivativef ( degree, 1, 4, 0, (float*)ctlpoints, \
    0, (float*)dctlpoints )


void mbs_multiFindBSDerivativef ( int degree, int lastknot, const float *knots, 
                                  int ncurves, int spdimen,  
                                  int pitch, const float *ctlpoints,  
                                  int *lastdknot, float *dknots,  
                                  int dpitch, float *dctlpoints );

#define mbs_FindBSDerivativeC1f(degree,lastknot,knots,coeff, \
    lastdknot,dknots,dcoeff) \
  mbs_multiFindBSDerivativef ( degree, lastknot, knots, 1, 1, 0, coeff, \
    lastdknot, dknots, 0, dcoeff )
#define mbs_FindBSDerivativeC2f(degree,lastknot,knots,ctlpoints, \
    lastdknot,dknots,dctlpoints) \
  mbs_multiFindBSDerivativef ( degree, lastknot, knots, 1, 2, 0, \
    (float*)ctlpoints, lastdknot, dknots, 0, (float*)dctlpoints )
#define mbs_FindBSDerivativeC3f(degree,lastknot,knots,ctlpoints, \
    lastdknot,dknots,dctlpoints) \
  mbs_multiFindBSDerivativef ( degree, lastknot, knots, 1, 3, 0, \
    (float*)ctlpoints, lastdknot, dknots, 0, (float*)dctlpoints )
#define mbs_FindBSDerivativeC4f(degree,lastknot,knots,ctlpoints, \
    lastdknot,dknots,dctlpoints) \
  mbs_multiFindBSDerivativef ( degree, lastknot, knots, 1, 4, 0, \
    (float*)ctlpoints, lastdknot, dknots, 0, (float*)dctlpoints )


boolean mbs_FindBSCommonKnotSequencef ( int *degree, int *lastknot,
            float **knots, int nsequences, ... );
boolean mbs_multiAdjustBSCRepf ( int ncurves, int spdimen,
     int indegree, int inlastknot, const float *inknots,
     int inpitch, const float *inctlpoints,
     int outdegree, int outlastknot, CONST_ float *outknots,
     int outpitch, float *outctlpoints );

#define mbs_AdjustBSCRepC1f(indegree,inlastknot,inknots, \
    inctlpoints,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiAdjustBSCRepf (1,1,indegree,inlastknot,inknots,0, \
    inctlpoints,outdegree,outlastknot,outknots,0,outctlpoints)
#define mbs_AdjustBSCRepC2f(indegree,inlastknot,inknots, \
    inctlpoints,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiAdjustBSCRepf (1,2,indegree,inlastknot,inknots,0, \
    (float*)inctlpoints,outdegree,outlastknot,outknots,0,(float*)outctlpoints)
#define mbs_AdjustBSCRepC3f(indegree,inlastknot,inknots, \
    inctlpoints,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiAdjustBSCRepf (1,3,indegree,inlastknot,inknots,0, \
    (float*)inctlpoints,outdegree,outlastknot,outknots,0,(float*)outctlpoints)
#define mbs_AdjustBSCRepC4f(indegree,inlastknot,inknots, \
    inctlpoints,outdegree,outlastknot,outknots,outctlpoints) \
  mbs_multiAdjustBSCRepf (1,4,indegree,inlastknot,inknots,0, \
    (float*)inctlpoints,outdegree,outlastknot,outknots,0,(float*)outctlpoints)


void mbs_multiAddBSCurvesf ( int ncurves, int spdimen,
                             int degree1, int lastknot1, CONST_ float *knots1,
                             int pitch1, CONST_ float *ctlpoints1,
                             int degree2, int lastknot2, CONST_ float *knots2,
                             int pitch2, CONST_ float *ctlpoints2,
                             int *sumdeg, int *sumlastknot, float *sumknots,
                             int sumpitch, float *sumctlpoints );

#define mbs_AddBSCurvesC1f(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiAddBSCurvesf (1,1,degree1,lastknot1,knots1,0,ctlpoints1, \
    degree2,lastknot2,knots2,0,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,sumctlpoints)
#define mbs_AddBSCurvesC2f(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiAddBSCurvesf (1,2,degree1,lastknot1,knots1,0,(float*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(float*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(float*)sumctlpoints)
#define mbs_AddBSCurvesC3f(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiAddBSCurvesf (1,3,degree1,lastknot1,knots1,0,(float*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(float*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(float*)sumctlpoints)
#define mbs_AddBSCurvesC4f(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiAddBSCurvesf (1,4,degree1,lastknot1,knots1,0,(float*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(float*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(float*)sumctlpoints)


void mbs_multiSubtractBSCurvesf ( int ncurves, int spdimen,
                             int degree1, int lastknot1, CONST_ float *knots1,
                             int pitch1, CONST_ float *ctlpoints1,
                             int degree2, int lastknot2, CONST_ float *knots2,
                             int pitch2, CONST_ float *ctlpoints2,
                             int *sumdeg, int *sumlastknot, float *sumknots,  
                             int sumpitch, float *sumctlpoints );

#define mbs_SubtractBSCurvesC1f(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiSubtractBSCurvesf (1,1,degree1,lastknot1,knots1,0,ctlpoints1, \
    degree2,lastknot2,knots2,0,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,sumctlpoints)
#define mbs_SubtractBSCurvesC2f(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiSubtractBSCurvesf (1,2,degree1,lastknot1,knots1,0,(float*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(float*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(float*)sumctlpoints)
#define mbs_SubtractBSCurvesC3f(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiSubtractBSCurvesf (1,3,degree1,lastknot1,knots1,0,(float*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(float*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(float*)sumctlpoints)
#define mbs_SubtractBSCurvesC4f(degree1,lastknot1,knots1,ctlpoints1, \
    degree2,lastknot2,knots2,ctlpoints2, \
    sumdeg,sumlastknot,sumknots,sumctlpoints) \
  mbs_multiSubtractBSCurvesf (1,4,degree1,lastknot1,knots1,0,(float*)ctlpoints1, \
    degree2,lastknot2,knots2,0,(float*)ctlpoints2, \
    sumdeg,sumlastknot,sumknots,0,(float*)sumctlpoints)


boolean mbs_FindPolynomialZerosf ( int degree, const float *coeff,
                                   int *nzeros, float *zeros, float eps );


void mbs_ClipBC2f ( int ncplanes, const vector3f *cplanes,
                    int degree, const point2f *cpoints,
                    void (*output) (int degree, const point2f *cpoints) );
void mbs_ClipBC2Rf ( int ncplanes, const vector3f *cplanes,
                     int degree, const point3f *cpoints,   
                     void (*output) (int degree, const point3f *cpoints) );
void mbs_ClipBC3f ( int ncplanes, const vector4f *cplanes,
                    int degree, const point3f *cpoints,
                    void (*output) (int degree, const point3f *cpoints) );
void mbs_ClipBC3Rf ( int ncplanes, const vector4f *cplanes,
                     int degree, const point4f *cpoints,
                     void (*output) (int degree, const point4f *cpoints) );

/* ///////////////////////////////////////////////////////////////////////// */
/* Bicubic polynomial Coons patches */
void mbs_BezC1CoonsFindCornersf ( int spdimen,
                                  int degc00, const float *c00,
                                  int degc01, const float *c01,
                                  int degc10, const float *c10,
                                  int degc11, const float *c11,
                                  float *pcorners );
boolean mbs_BezC1CoonsToBezf ( int spdimen,
                               int degc00, const float *c00,
                               int degc01, const float *c01,
                               int degc10, const float *c10,
                               int degc11, const float *c11,
                               int degd00, const float *d00,
                               int degd01, const float *d01,
                               int degd10, const float *d10,
                               int degd11, const float *d11,
                               int *n, int *m, float *p );

void mbs_TabCubicHFuncDer2f ( float a, float b, int nkn, const float *kn,
                              float *hfunc, float *dhfunc, float *ddhfunc );
void mbs_TabCubicHFuncDer3f ( float a, float b, int nkn, const float *kn, 
                              float *hfunc, float *dhfunc, float *ddhfunc,
                              float *dddhfunc );
void mbs_TabBezCurveDer2f ( int spdimen, int degree, const float *cp,
                            int nkn, const float *kn,
                            int ppitch,
                            float *p, float *dp, float *ddp );
boolean _mbs_TabBezC1Coonsf ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp );
boolean mbs_TabBezC1CoonsDer2f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degc10, const float *c10,
      int degc11, const float *c11,
      int degd00, const float *d00,
      int degd01, const float *d01,
      int degd10, const float *d10,
      int degd11, const float *d11,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv );
boolean mbs_TabBezC1CoonsDer3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degc10, const float *c10,
      int degc11, const float *c11,
      int degd00, const float *d00,
      int degd01, const float *d01,
      int degd10, const float *d10,
      int degd11, const float *d11,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv );

boolean _mbs_TabBezC1Coons0f (
                   int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp );
boolean mbs_TabBezC1Coons0Der2f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degd00, const float *d00,
      int degd01, const float *d01,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv );
boolean mbs_TabBezC1Coons0Der3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degd00, const float *d00,
      int degd01, const float *d01,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv );

/* ///////////////////////////////////////////////////////////////////////// */
/* Biquintic polynomial Coons patches */
void mbs_BezC2CoonsFindCornersf ( int spdimen,
                                  int degc00, const float *c00,
                                  int degc01, const float *c01,
                                  int degc02, const float *c02,
                                  int degc10, const float *c10,
                                  int degc11, const float *c11,
                                  int degc12, const float *c12,
                                  float *pcorners );
boolean mbs_BezC2CoonsToBezf ( int spdimen,
                               int degc00, const float *c00,
                               int degc01, const float *c01,
                               int degc02, const float *c02,
                               int degc10, const float *c10,
                               int degc11, const float *c11,
                               int degc12, const float *c12,
                               int degd00, const float *d00,
                               int degd01, const float *d01,
                               int degd02, const float *d02,
                               int degd10, const float *d10,
                               int degd11, const float *d11,
                               int degd12, const float *d12,
                               int *n, int *m, float *p );
void mbs_TabQuinticHFuncDer3f ( float a, float b, int nkn, const float *kn,
                                float *hfunc, float *dhfunc,
                                float *ddhfunc, float *dddhfunc );
void mbs_TabBezCurveDer3f ( int spdimen, int degree, const float *cp,
                            int nkn, const float *kn,
                            int ppitch,
                            float *p, float *dp, float *ddp, float *dddp );
boolean _mbs_TabBezC2Coonsf (
                   int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp );

boolean mbs_TabBezC2CoonsDer3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degc02, const float *c02,
      int degc10, const float *c10,
      int degc11, const float *c11,
      int degc12, const float *c12,
      int degd00, const float *d00,
      int degd01, const float *d01,
      int degd02, const float *d02,
      int degd10, const float *d10,
      int degd11, const float *d11,
      int degd12, const float *d12,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv );

boolean _mbs_TabBezC2Coons0f (
                   int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp );
boolean mbs_TabBezC2Coons0Der3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degc02, const float *c02,
      int degd00, const float *d00,
      int degd01, const float *d01,
      int degd02, const float *d02,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv );


/* ///////////////////////////////////////////////////////////////////////// */
/* Bicubic B-spline Coons patches */
void mbs_BSC1CoonsFindCornersf ( int spdimen,
          int degc00, int lastknotc00, const float *knotsc00, const float *c00,
          int degc01, int lastknotc01, const float *knotsc01, const float *c01,
          int degc10, int lastknotc10, const float *knotsc10, const float *c10,
          int degc11, int lastknotc11, const float *knotsc11, const float *c11,
          float *pcorners );
boolean mbs_BSC1CoonsToBSf ( int spdimen,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degc10, int lastknotc10, const float *knotsc10, const float *c10,
      int degc11, int lastknotc11, const float *knotsc11, const float *c11,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      int degd10, int lastknotd10, const float *knotsd10, const float *d10,
      int degd11, int lastknotd11, const float *knotsd11, const float *d11,
      int *degreeu, int *lastuknot, float *uknots,
      int *degreev, int *lastvknot, float *vknots, float *p );

void mbs_TabBSCurveDer2f ( int spdimen, int degree, int lastknot,
                           const float *knots, const float *cp,
                           int nkn, const float *kn, int ppitch,
                           float *p, float *dp, float *ddp );
boolean _mbs_TabBSC1Coonsf ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp );
boolean mbs_TabBSC1CoonsDer2f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degc10, int lastknotc10, const float *knotsc10, const float *c10,
      int degc11, int lastknotc11, const float *knotsc11, const float *c11,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      int degd10, int lastknotd10, const float *knotsd10, const float *d10,
      int degd11, int lastknotd11, const float *knotsd11, const float *d11,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv );
boolean mbs_TabBSC1CoonsDer3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degc10, int lastknotc10, const float *knotsc10, const float *c10,
      int degc11, int lastknotc11, const float *knotsc11, const float *c11,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      int degd10, int lastknotd10, const float *knotsd10, const float *d10,
      int degd11, int lastknotd11, const float *knotsd11, const float *d11,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv );

boolean _mbs_TabBSC1Coons0f ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp );
boolean mbs_TabBSC1Coons0Der2f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv );
boolean mbs_TabBSC1Coons0Der3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv );

/* ///////////////////////////////////////////////////////////////////////// */
/* Biquintic B-spline Coons patches */          
void mbs_BSC2CoonsFindCornersf ( int spdimen,
          int degc00, int lastknotc00, const float *knotsc00, const float *c00,
          int degc01, int lastknotc01, const float *knotsc01, const float *c01,
          int degc02, int lastknotc02, const float *knotsc02, const float *c02,
          int degc10, int lastknotc10, const float *knotsc10, const float *c10,
          int degc11, int lastknotc11, const float *knotsc11, const float *c11,
          int degc12, int lastknotc12, const float *knotsc12, const float *c12,
          float *pcorners );
boolean mbs_BSC2CoonsToBSf ( int spdimen,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degc02, int lastknotc02, const float *knotsc02, const float *c02,
      int degc10, int lastknotc10, const float *knotsc10, const float *c10,
      int degc11, int lastknotc11, const float *knotsc11, const float *c11,
      int degc12, int lastknotc12, const float *knotsc12, const float *c12,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      int degd02, int lastknotd02, const float *knotsd02, const float *d02,
      int degd10, int lastknotd10, const float *knotsd10, const float *d10,
      int degd11, int lastknotd11, const float *knotsd11, const float *d11,
      int degd12, int lastknotd12, const float *knotsd12, const float *d12,
      int *degreeu, int *lastuknot, float *uknots,
      int *degreev, int *lastvknot, float *vknots, float *p );
void mbs_TabBSCurveDer3f ( int spdimen, int degree, int lastknot,
                           const float *knots, const float *cp,
                           int nkn, const float *kn, int ppitch,
                           float *p, float *dp, float *ddp, float *dddp );
boolean _mbs_TabBSC2Coonsf ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp );
boolean mbs_TabBSC2CoonsDer3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degc02, int lastknotc02, const float *knotsc02, const float *c02,
      int degc10, int lastknotc10, const float *knotsc10, const float *c10,
      int degc11, int lastknotc11, const float *knotsc11, const float *c11,
      int degc12, int lastknotc12, const float *knotsc12, const float *c12,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      int degd02, int lastknotd02, const float *knotsd02, const float *d02,
      int degd10, int lastknotd10, const float *knotsd10, const float *d10,
      int degd11, int lastknotd11, const float *knotsd11, const float *d11,
      int degd12, int lastknotd12, const float *knotsd12, const float *d12,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv );
boolean _mbs_TabBSC2Coons0f ( int spdimen, int nknu, int nknv,
                   const float *c, const float *d, const float *p,
                   const float *hu, const float *hv, float *pp );
boolean mbs_TabBSC2Coons0Der3f ( int spdimen,
      int nknu, const float *knu, const float *hfuncu,
      const float *dhfuncu, const float *ddhfuncu, const float *dddhfuncu,
      int nknv, const float *knv, const float *hfuncv,
      const float *dhfuncv, const float *ddhfuncv, const float *dddhfuncv,
      int degc00, int lastknotc00, const float *knotsc00, const float *c00,
      int degc01, int lastknotc01, const float *knotsc01, const float *c01,
      int degc02, int lastknotc02, const float *knotsc02, const float *c02,
      int degd00, int lastknotd00, const float *knotsd00, const float *d00,
      int degd01, int lastknotd01, const float *knotsd01, const float *d01,
      int degd02, int lastknotd02, const float *knotsd02, const float *d02,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv,
      float *puuu, float *puuv, float *puvv, float *pvvv );


/* ///////////////////////////////////////////////////////////////////////// */
/* spherical product of two curves */
void mbs_SphericalProductf (
                int degree_eq, int lastknot_eq, const point2f *cpoints_eq,
                int degree_mer, int lastknot_mer, const point2f *cpoints_mer,
                int pitch, point3f *spr_cp );
void mbs_SphericalProductRf (
                int degree_eq, int lastknot_eq, const point3f *cpoints_eq,
                int degree_mer, int lastknot_mer, const point3f *cpoints_mer,
                int pitch, point4f *spr_cp );

/* ///////////////////////////////////////////////////////////////////////// */
/* Lane-Riesenfeld algorithm */
boolean mbs_multiLaneRiesenfeldf ( int spdimen, int ncurves, int degree,
                         int inlastknot, int inpitch, const float *incp,
                         int *outlastknot, int outpitch, float *outcp );

#define mbs_LaneRiesenfeldC1f(degree,inlastknot,incp,outlastknot,outcp) \
  mbs_multiLaneRiesenfeldf ( 1, 1, degree, inlastknot, 0, incp, \
    outlastknot, 0, outcp )
#define mbs_LaneRiesenfeldC2f(degree,inlastknot,incp,outlastknot,outcp) \
  mbs_multiLaneRiesenfeldf ( 2, 1, degree, inlastknot, 0, (float*)incp, \
    outlastknot, 0, (float*)outcp )
#define mbs_LaneRiesenfeldC3f(degree,inlastknot,incp,outlastknot,outcp) \
  mbs_multiLaneRiesenfeldf ( 3, 1, degree, inlastknot, 0, (float*)incp, \
    outlastknot, 0, (float*)outcp )
#define mbs_LaneRiesenfeldC4f(degree,inlastknot,incp,outlastknot,outcp) \
  mbs_multiLaneRiesenfeldf ( 4, 1, degree, inlastknot, 0, (float*)incp, \
    outlastknot, 0, (float*)outcp )


#ifdef __cplusplus
}
#endif

#endif /*MULTIBSF_H*/

