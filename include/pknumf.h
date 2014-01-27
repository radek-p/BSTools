
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* Header file for the libpknum library of C procedures - numerical methods */

#ifndef CONST_  /* a dirty trick to suppress many compiler warning messages */
#define CONST_ const
#endif

#ifndef PKNUMF_H
#define PKNUMF_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
typedef struct bandm_profile {
  int firstnz;  /* index of the first row with a nonzero  */
                /* coefficient in the column */
  int ind;      /* index of this coefficient in the array */
} bandm_profile;

typedef struct { float x, y; } complexf;
#endif

/* ///////////////////////////////////////////////////////////////////////// */
double pkn_ScalarProductf ( int spdimen, const float *a, const float *b );
double pkn_SecondNormf ( int spdimen, const float *b );

double pkn_detf ( int n, float *a );

void pkn_AddMatrixf ( int nrows, int rowlen,
                      int inpitch1, const float *indata1,
                      int inpitch2, const float *indata2,
                      int outpitch, float *outdata );
void pkn_AddMatrixMf ( int nrows, int rowlen,
                       int inpitch1, const float *indata1,
                       int inpitch2, const float *indata2,
                       double a,
                       int outpitch, float *outdata );
void pkn_SubtractMatrixf ( int nrows, int rowlen,
                           int inpitch1, const float *indata1,
                           int inpitch2, const float *indata2,
                           int outpitch, float *outdata );
void pkn_MatrixMDifferencef ( int nrows, int rowlen,
                              int inpitch1, const float *indata1,
                              int inpitch2, const float *indata2,
                              double a,
                              int outpitch, float *outdata );
void pkn_MatrixLinCombf ( int nrows, int rowlen,
                          int inpitch1, const float *indata1,
                          double a,
                          int inpitch2, const float *indata2,
                          double b,
                          int outpitch, float *outdata );

void pkn_MultMatrixNumf ( int nrows, int rowlen,
                          int inpitch, const float *indata, 
                          double a,
                          int outpitch, float *outdata );
void pkn_MultArrayf ( int nrows, int rowlen, int pitch_a, CONST_ float *a,
                      int pitch_b, CONST_ float *b,
                      int pitch_c, float *c );

void pkn_MultMatrixf ( int nrows_a, int rowlen_a, int pitch_a, CONST_ float *a,
                       int rowlen_b, int pitch_b, CONST_ float *b,
                       int pitch_c, float *c );
void pkn_MultMatrixAddf ( int nrows_a, int rowlen_a, int pitch_a, CONST_ float *a,
                          int rowlen_b, int pitch_b, CONST_ float *b,
                          int pitch_c, float *c );
void pkn_MultMatrixSubf ( int nrows_a, int rowlen_a, int pitch_a, CONST_ float *a,
                          int rowlen_b, int pitch_b, CONST_ float *b,
                          int pitch_c, float *c );

void pkn_MultTMatrixf ( int nrows_a, int rowlen_a, int pitch_a, const float *a,
                        int rowlen_b, int pitch_b, const float *b,
                        int pitch_c, float *c );
void pkn_MultTMatrixAddf ( int nrows_a, int rowlen_a, int pitch_a,
                           const float *a,
                           int rowlen_b, int pitch_b, const float *b,
                           int pitch_c, float *c );
void pkn_MultTMatrixSubf ( int nrows_a, int rowlen_a, int pitch_a,
                           const float *a,
                           int rowlen_b, int pitch_b, const float *b,
                           int pitch_c, float *c );
void pkn_MultMatrixTf ( int nrows_a, int rowlen_a, int pitch_a, const float *a,
                        int nrows_b, int pitch_b, const float *b,
                        int pitch_c, float *c );
void pkn_MultMatrixTAddf ( int nrows_a, int rowlen_a, int pitch_a, const float *a,
                           int nrows_b, int pitch_b, const float *b,
                           int pitch_c, float *c );
void pkn_MultMatrixTSubf ( int nrows_a, int rowlen_a, int pitch_a, const float *a,
                           int nrows_b, int pitch_b, const float *b,
                           int pitch_c, float *c );

/* ///////////////////////////////////////////////////////////////////////// */
void pkn_FindGivensRotationf ( float a, float b, float *c, float *s );
void pkn_FindGivensRotXif ( float a, float b, float *c, float *s, float *xi );
void pkn_FindXiGivensRotf ( float xi, float *c, float *s );
void pkn_ApplyGivensRotationf ( float c, float s, float *a, float *b );
void pkn_ApplySymGivensRotationf ( float c, float s,
                                   float *d, float *e, float *f );

/* ///////////////////////////////////////////////////////////////////////// */
void pkn_MVectorSumf ( int m, int n, float *sum, ... );
void pkn_MVectorLinCombf ( int m, int n, float *sum, ... );

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_GaussDecomposePLUQf ( int n, float *a, int *P, int *Q );
void pkn_multiSolvePLUQf ( int n, const float *lu, const int *P, const int *Q,
                           int spdimen, int pitch, float *b );
boolean pkn_multiGaussSolveLinEqf ( int n, const float *a,
                                    int spdimen, int pitch, float *b );
boolean pkn_GaussInvertMatrixf ( int n, float *a );

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_QRDecomposeMatrixf ( int nrows, int ncols, float *a, float *aa );
void pkn_multiReflectVectorf ( int nrows, int ncols,
                               const float *a, const float *aa,
                               int spdimen, int pitch, float *b );
void pkn_multiInvReflectVectorf ( int nrows, int ncols,
                                  const float *a, const float *aa,
                                  int spdimen, int pitch, float *b );
void pkn_multiMultUTVectorf ( int nrows, const float *a,
                              int spdimen, int bpitch, float *b,
                              int xpitch, float *x );
void pkn_multiMultInvUTVectorf ( int nrows, const float *a,
                                 int spdimen, int bpitch, float *b,
                                 int xpitch, float *x );
void pkn_multiMultTrUTVectorf ( int nrows, const float *a,
                                int spdimen, int bpitch, float *b,   
                                int xpitch, float *x );
void pkn_multiMultInvTrUTVectorf ( int nrows, const float *a,
                                   int spdimen, int bpitch, float *b,   
                                   int xpitch, float *x );
boolean pkn_multiSolveRLSQf ( int nrows, int ncols, float *a,
                              int spdimen, int bpitch, float *b,
                              int xpitch, float *x );
void pkn_QRGetReflectionf ( int nrows, int ncols,
                            const float *a, const float *aa,
                            int nrefl, float *w, float *gamma );


/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
void pkn_BandmFindQRMSizes ( int ncols, const bandm_profile *aprof,
                             int *qsize, int *rsize );
#endif

void pkn_BandmQRDecomposeMatrixf ( int nrows, int ncols,
                                  const bandm_profile *aprof, const float *a,
                                  bandm_profile *qprof, float *q,
                                  bandm_profile *rprof, float *r );
void pkn_multiBandmReflectVectorf ( int ncols,
                                    const bandm_profile *qprof, const float *q,
                                    int spdimen, float *b );
void pkn_multiBandmInvReflectVectorf ( int ncols,
                                       const bandm_profile *qprof,
                                       const float *q,
                                       int spdimen, float *b );
void pkn_multiBandmMultVectorf ( int nrows, int ncols,
                                 const bandm_profile *aprof, const float *a,
                                 int spdimen, const float *x, float *y );
void pkn_multiBandmMultInvUTMVectorf ( int nrows,
                                       const bandm_profile *rprof,
                                       const float *r,
                                       int spdimen, const float *x, float *y );
void pkn_multiBandmMultTrVectorf ( int ncols,
                                   const bandm_profile *aprof, const float *a,
                                   int spdimen, const float *x, float *y );
void pkn_multiBandmMultInvTrUTMVectorf ( int nrows,
                                         const bandm_profile *rprof,
                                         const float *r,
                                         int spdimen, const float *x, float *y );

boolean pkn_multiBandmSolveRLSQf ( int nrows, int ncols,
                                   const bandm_profile *aprof, const float *a,
                                   int nrsides, int spdimen,
                                   int bpitch, const float *b,
                                   int xpitch, float *x );
boolean pkn_multiBandmSolveDLSQf ( int nrows, int ncols,
                                   const bandm_profile *atprof, const float *at,
                                   int nrsides, int spdimen,
                                   int bpitch, const float *b,
                                   int x0pitch, const float *x0,
                                   int xpitch, float *x );
boolean pkn_multiBandmSolveCRLSQf ( int nrows, int ncols,
                                    const bandm_profile *aprof, const float *a,
                                    int nconstr, int cpitch, const float *c,
                                    int nrsides, int spdimen,
                                    int bpitch, const float *b,
                                    int dpitch, const float *d,
                                    int xpitch, float *x );

void pkn_PrintBandmf ( int ncols, const bandm_profile *aprof, const float *a );
void pkn_PrintBandmRowSumf ( int ncols, const bandm_profile *aprof, const float *a );
void pkn_PrintMatf ( int nrows, int ncols, const float *a );

#ifndef PKNUM_H
void pkn_PrintProfile ( int ncols, const bandm_profile *prof );
#endif

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
#define pkn_LowerTrMatIndex(i,j) \
  ( (i)*((i)+1)/2+(j) )
#define pkn_SymMatIndex(i,j) \
  ( (i) >= (j) ? (i)*((i)+1)/2+(j) : (j)*((j)+1)/2+(i) )
#endif

boolean pkn_CholeskyDecompf ( int n, float *a );
void pkn_SymMatrixMultf ( int n, CONST_ float *a, int spdimen,
                          int bpitch, CONST_ float *b,
                          int xpitch, float *x );
void pkn_LowerTrMatrixMultf ( int n, CONST_ float *l, int spdimen,
                              int bpitch, CONST_ float *b,
                              int xpitch, float *x );
void pkn_UpperTrMatrixMultf ( int n, CONST_ float *l, int spdimen,
                              int bpitch, CONST_ float *b,
                              int xpitch, float *x );
void pkn_LowerTrMatrixSolvef ( int n, const float *l, int spdimen,
                               int bpitch, const float *b,
                               int xpitch, float *x );
void pkn_UpperTrMatrixSolvef ( int n, const float *l, int spdimen,
                               int bpitch, const float *b,
                               int xpitch, float *x );

void pkn_MatrixLowerTrMultf ( int m, int n, int bpitch, CONST_ float *b,
                              CONST_ float *l, int xpitch, float *x );
void pkn_MatrixUpperTrMultf ( int m, int n, int bpitch, CONST_ float *b,
                              CONST_ float *l, int xpitch, float *x );
void pkn_MatrixLowerTrSolvef ( int m, int n, int bpitch, CONST_ float *b,
                               CONST_ float *l, int xpitch, float *x );
void pkn_MatrixUpperTrSolvef ( int m, int n, int bpitch, CONST_ float *b,
                               CONST_ float *l, int xpitch, float *x );
void pkn_MatrixLowerTrMultAddf ( int m, int n, int bpitch, CONST_ float *b,
                                 CONST_ float *l, int xpitch, float *x );
void pkn_MatrixUpperTrMultAddf ( int m, int n, int bpitch, CONST_ float *b,
                                 CONST_ float *l, int xpitch, float *x );
boolean pkn_MatrixLowerTrSolveAddf ( int m, int n, int bpitch, CONST_ float *b,
                                     CONST_ float *l, int xpitch, float *x );
boolean pkn_MatrixUpperTrSolveAddf ( int m, int n, int bpitch, CONST_ float *b,
                                     CONST_ float *l, int xpitch, float *x );
void pkn_MatrixLowerTrMultSubf ( int m, int n, int bpitch, CONST_ float *b,
                                 CONST_ float *l, int xpitch, float *x );
void pkn_MatrixUpperTrMultSubf ( int m, int n, int bpitch, CONST_ float *b,
                                 CONST_ float *l, int xpitch, float *x );
boolean pkn_MatrixLowerTrSolveSubf ( int m, int n, int bpitch, CONST_ float *b,
                                     CONST_ float *l, int xpitch, float *x );
boolean pkn_MatrixUpperTrSolveSubf ( int m, int n, int bpitch, CONST_ float *b,
                                     CONST_ float *l, int xpitch, float *x );

void pkn_SymToFullMatrixf ( int n, const float *syma,
                            int pitch, float *fulla );
void pkn_FullToSymMatrixf ( int n, int pitch, const float *fulla,
                            float *syma );
#define pkn_FullToLTrMatrixf(n,pitch,fulla,ltra) \
  pkn_FullToSymMatrixf(n,pitch,fulla,ltra)
void pkn_LTrToFullMatrixf ( int n, const float *ltra,
                            int pitch, float *fulla );
void pkn_UTrToFullMatrixf ( int n, const float *utra,
                            int pitch, float *fulla );
void pkn_FullToUTrMatrixf ( int n, int pitch, const float *fulla,
                            float *utra );

boolean pkn_ComputeQSQTf ( int m, const float *s,   
                           int n, const float *a, const float *aa,
                           float *b );
boolean pkn_ComputeQTSQf ( int m, const float *s,
                           int n, const float *a, const float *aa,
                           float *b );

void pkn_SymMatSubAATf ( int n, float *b, int m, int pitch_a, CONST_ float *a );

void pkn_SymMatFindEigenvalueIntervalf ( int n, float *a,
                                         float *lmin, float *lmax );
boolean pkn_SymMatFindEigenvaluesf ( int n, float *a, float *eigenval );

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
#define pkn_LHessenbergMatIndex(i,j) \
  (pkn_LowerTrMatIndex(((i)+1),(j))-1)
#define pkn_UHessenbergMatIndex(i,j) \
  (pkn_LowerTrMatIndex(((j)+1),(i))-1)
#endif

boolean pkn_QRDecompUHessenbergf ( int n, float *ah );
void pkn_multiSolveQRUHessenbergf ( int n, float *qrh,
                                    int spdimen, int pitch, float *b );
boolean pkn_multiSolveUHessenbergLinEqf ( int n, float *ah,
                                          int spdimen, int pitch, float *b );

/* ///////////////////////////////////////////////////////////////////////// */
float pkn_Illinoisf ( float (*f) (float), float a, float b, float eps,
                      boolean *error );

boolean pkn_SolveSqEqf ( float p, float q, float *x1, float *x2 );

float pkn_GoldenRatf ( float (*f) (float), float a, float b, float eps,     
                       boolean *error );


/* ///////////////////////////////////////////////////////////////////////// */
void pkn_Setup2DerA11Matrixf ( float xu, float yu, float xv, float yv,
                               float *A11 );

void pkn_Setup2DerA21Matrixf ( float xuu, float yuu, float xuv,
                               float yuv, float xvv, float yvv,
                               float *A21 );
void pkn_Setup2DerA22Matrixf ( float xu, float yu, float xv, float yv,
                               float *A22 );

void pkn_Setup2DerA31Matrixf ( float xuuu, float yuuu, float xuuv, float yuuv,
                               float xuvv, float yuvv, float xvvv, float yvvv,
                               float *A31 );
void pkn_Setup2DerA32Matrixf ( float xu, float yu, float xv, float yv,
                               float xuu, float yuu, float xuv,
                               float yuv, float xvv, float yvv,
                               float *A32 );
void pkn_Setup2DerA33Matrixf ( float xu, float yu, float xv, float yv,
                               float *A33 );

void pkn_Setup2DerA41Matrixf ( float xuuuu, float yuuuu, float xuuuv, float yuuuv,
                               float xuuvv, float yuuvv, float xuvvv, float yuvvv,
                               float xvvvv, float yvvvv,
                               float *A41 );
void pkn_Setup2DerA42Matrixf ( float xu, float yu, float xv, float yv,
                               float xuu, float yuu, float xuv,
                               float yuv, float xvv, float yvv,
                               float xuuu, float yuuu, float xuuv, float yuuv,
                               float xuvv, float yuvv, float xvvv, float yvvv,
                               float *A42 );
void pkn_Setup2DerA43Matrixf ( float xu, float yu, float xv, float yv,
                               float xuu, float yuu, float xuv,
                               float yuv, float xvv, float yvv,
                               float *A43 );
void pkn_Setup2DerA44Matrixf ( float xu, float yu, float xv, float yv,
                               float *A44 );


boolean pkn_Comp2Derivatives1f ( float xu, float yu, float xv, float yv,
                                 int spdimen, const float *gx, const float *gy,
                                 float *hu, float *hv );
boolean pkn_Comp2Derivatives2f ( float xu, float yu, float xv, float yv,
                                 float xuu, float yuu, float xuv,
                                 float yuv, float xvv, float yvv,
                                 int spdimen,
                                 const float *gx, const float *gy,
                                 const float *gxx, const float *gxy,
                                 const float *gyy,
                                 float *huu, float *huv, float *hvv );
boolean pkn_Comp2Derivatives3f ( float xu, float yu, float xv, float yv,
                                 float xuu, float yuu, float xuv,
                                 float yuv, float xvv, float yvv,
                                 float xuuu, float yuuu, float xuuv, float yuuv,
                                 float xuvv, float yuvv, float xvvv, float yvvv,
                                 int spdimen,
                                 const float *gx, const float *gy,
                                 const float *gxx, const float *gxy,
                                 const float *gyy,
                                 const float *gxxx, const float *gxxy,
                                 const float *gxyy, const float *gyyy,
                                 float *huuu, float *huuv,
                                 float *huvv, float *hvvv );
boolean pkn_Comp2Derivatives4f ( float xu, float yu, float xv, float yv,
                                 float xuu, float yuu, float xuv,
                                 float yuv, float xvv, float yvv,
                                 float xuuu, float yuuu, float xuuv, float yuuv,
                                 float xuvv, float yuvv, float xvvv, float yvvv,
                                 float xuuuu, float yuuuu, float xuuuv,
                                 float yuuuv, float xuuvv, float yuuvv,
                                 float xuvvv, float yuvvv,
                                 float xvvvv, float yvvvv,
                                 int spdimen,
                                 const float *gx, const float *gy,
                                 const float *gxx, const float *gxy,
                                 const float *gyy,
                                 const float *gxxx, const float *gxxy,
                                 const float *gxyy, const float *gyyy,
                                 const float *gxxxx, const float *gxxxy,
                                 const float *gxxyy, const float *gxyyy,
                                 const float *gyyyy,
                                 float *huuuu, float *huuuv, float *huuvv,
                                 float *huvvv, float *hvvvv );

boolean pkn_Comp2iDerivatives1f ( float xu, float yu, float xv, float yv,
                                  int spdimen, const float *hu, const float *hv,
                                  float *gx, float *gy );
boolean pkn_Comp2iDerivatives2f ( float xu, float yu, float xv, float yv,
                                  float xuu, float yuu, float xuv,
                                  float yuv, float xvv, float yvv,
                                  int spdimen,
                                  const float *hu, const float *hv,
                                  const float *huu, const float *huv,
                                  const float *hvv,
                                  float *gx, float *gy,
                                  float *gxx, float *gxy, float *gyy );
boolean pkn_Comp2iDerivatives3f ( float xu, float yu, float xv, float yv,
                                  float xuu, float yuu, float xuv,
                                  float yuv, float xvv, float yvv,
                                  float xuuu, float yuuu, float xuuv, float yuuv,
                                  float xuvv, float yuvv, float xvvv, float yvvv,
                                  int spdimen,
                                  const float *hu, const float *hv,
                                  const float *huu, const float *huv,
                                  const float *hvv,
                                  const float *huuu, const float *huuv,
                                  const float *huvv, const float *hvvv,
                                  float *gx, float *gy,
                                  float *gxx, float *gxy, float *gyy,
                                  float *gxxx, float *gxxy,
                                  float *gxyy, float *gyyy );
boolean pkn_Comp2iDerivatives4f ( float xu, float yu, float xv, float yv,
                                  float xuu, float yuu, float xuv,
                                  float yuv, float xvv, float yvv,
                                  float xuuu, float yuuu, float xuuv, float yuuv,
                                  float xuvv, float yuvv, float xvvv, float yvvv,
                                  float xuuuu, float yuuuu, float xuuuv,
                                  float yuuuv, float xuuvv, float yuuvv,
                                  float xuvvv, float yuvvv,
                                  float xvvvv, float yvvvv,
                                  int spdimen,
                                  const float *hu, const float *hv,
                                  const float *huu, const float *huv,
                                  const float *hvv,
                                  const float *huuu, const float *huuv,
                                  const float *huvv, const float *hvvv,
                                  const float *huuuu, const float *huuuv,
                                  const float *huuvv, const float *huvvv,
                                  const float *hvvvv,
                                  float *gx, float *gy,
                                  float *gxx, float *gxy, float *gyy,
                                  float *gxxx, float *gxxy,
                                  float *gxyy, float *gyyy,
                                  float *gxxxx, float *gxxxy, float *gxxyy,
                                  float *gxyyy, float *gyyyy );


boolean pkn_f2iDerivatives1f ( float xu, float yu, float xv, float yv,
                              float *gx, float *gy );
boolean pkn_f2iDerivatives2f ( float xu, float yu, float xv, float yv,
                               float xuu, float yuu, float xuv,
                               float yuv, float xvv, float yvv,
                               float *gx, float *gy,
                               float *gxx, float *gxy, float *gyy );
boolean pkn_f2iDerivatives3f ( float xu, float yu, float xv, float yv,
                               float xuu, float yuu, float xuv,
                               float yuv, float xvv, float yvv,
                               float xuuu, float yuuu, float xuuv, float yuuv,
                               float xuvv, float yuvv, float xvvv, float yvvv,
                               float *gx, float *gy,
                               float *gxx, float *gxy, float *gyy,
                               float *gxxx, float *gxxy,
                               float *gxyy, float *gyyy );
boolean pkn_f2iDerivatives4f ( float xu, float yu, float xv, float yv,
                               float xuu, float yuu, float xuv,
                               float yuv, float xvv, float yvv,
                               float xuuu, float yuuu, float xuuv, float yuuv,
                               float xuvv, float yuvv, float xvvv, float yvvv,
                               float xuuuu, float yuuuu, float xuuuv,
                               float yuuuv, float xuuvv, float yuuvv,
                               float xuvvv, float yuvvv,
                               float xvvvv, float yvvvv,
                               float *gx, float *gy,
                               float *gxx, float *gxy, float *gyy,
                               float *gxxx, float *gxxy,
                               float *gxyy, float *gyyy,
                               float *gxxxx, float *gxxxy, float *gxxyy,
                               float *gxyyy, float *gyyyy );

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
int pkn_Block1ArraySize ( int k, int r, int s );
int pkn_Block1FindBlockPos ( int k, int r, int s, int i, int j );
int pkn_Block1FindElemPos ( int k, int r, int s, int i, int j );
#endif

boolean pkn_Block1CholeskyDecompMf ( int k, int r, int s, float *A );
void pkn_Block1LowerTrMSolvef ( int k, int r, int s, CONST_ float *A,
                                int spdimen, int xpitch, float *x );
void pkn_Block1UpperTrMSolvef ( int k, int r, int s, CONST_ float *A,
                                int spdimen, int xpitch, float *x );
boolean pkn_Block1SymMatrixMultf ( int k, int r, int s, CONST_ float *A,
                                   int spdimen, int xpitch, float *x,
                                   int ypitch, float *y );


/* ///////////////////////////////////////////////////////////////////////// */
int pkn_FFTPermutef ( int n, int rowlen, int pitch, complexf *a );
boolean pkn_FFTf ( int n, int rowlen, int pitch, complexf *a );
boolean pkn_InvFFTf ( int n, int rowlen, int pitch, complexf *a );


/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
int pkn_Block2ArraySize ( int k, int r, int s, int t );
int pkn_Block2FindBlockPos ( int k, int r, int s, int t, int i, int j );
int pkn_Block2FindElemPos ( int k, int r, int s, int t, int i, int j );
#endif

boolean pkn_Block2CholeskyDecompMf ( int k, int r, int s, int t, float *A );
void pkn_Block2LowerTrMSolvef ( int k, int r, int s, int t, CONST_ float *L,
                                int spdimen, int xpitch, float *x );
void pkn_Block2UpperTrMSolvef ( int k, int r, int s, int t, CONST_ float *L,
                                int spdimen, int xpitch, float *x );
void pkn_Block2SymMatrixMultf ( int k, int r, int s, int t, CONST_ float *A,
                                int spdimen, int xpitch, CONST_ float *x,
                                int ypitch, float *y );

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
int pkn_Block3ArraySize ( int k, int r, int s );
int pkn_Block3FindBlockPos ( int k, int r, int s, int i, int j );
int pkn_Block3FindElemPos ( int k, int r, int s, int i, int j );
#endif

boolean pkn_Block3CholeskyDecompMf ( int k, int r, int s, float *A );
void pkn_Block3LowerTrMSolvef ( int k, int r, int s, CONST_ float *L,
                                int spdimen, int xpitch, float *x );
void pkn_Block3UpperTrMSolvef ( int k, int r, int s, CONST_ float *L,
                                int spdimen, int xpitch, float *x );

void pkn_Block3SymMatrixMultf ( int k, int r, int s, CONST_ float *A,
                                int spdimen, int xpitch, CONST_ float *x,
                                int ypitch, float *y );

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKNUM_H
int pkn_NRBArraySize ( int n, const int *prof );
#endif

boolean pkn_NRBFindRowsf ( int n, const int *prof, CONST_ float *a,
                           float **row );
boolean pkn_NRBSymCholeskyDecompf ( int n, const int *prof, float *a,
                                    float **row, boolean *abort );
boolean pkn_NRBSymMultf ( int n, const int *prof, CONST_ float *a,
                          float **row,
                          int spdimen, int xpitch, CONST_ float *x,
                          int ypitch, float *y );
boolean pkn_NRBLowerTrMultf ( int n, const int *prof, CONST_ float *a,
                              float **row,
                              int spdimen, int xpitch, CONST_ float *x,
                              int ypitch, float *y );
boolean pkn_NRBUpperTrMultf ( int n, const int *prof, CONST_ float *a,
                              float **row,
                              int spdimen, int xpitch, CONST_ float *x,
                              int ypitch, float *y );
boolean pkn_NRBLowerTrSolvef ( int n, const int *prof, CONST_ float *l,
                               float **row,
                               int spdimen, int bpitch, CONST_ float *b,
                               int xpitch, float *x );
boolean pkn_NRBUpperTrSolvef ( int n, const int *prof, CONST_ float *l,
                               float **row,
                               int spdimen, int bpitch, CONST_ float *b,
                               int xpitch, float *x );

boolean pkn_NRBSymFindEigenvalueIntervalf ( int n, const int *prof,
                                            float *a, float **row,
                                            float *amin, float *amax );

boolean pkn_NRBComputeQTSQf ( int n, int *prof, float *Amat, float **Arows,
                              int w, float *Bmat, float *bb,
                              int *qaprof, float **QArows );
boolean pkn_NRBComputeQSQTf ( int n, int *prof, float *Amat, float **Arows,
                              int w, float *Bmat, float *bb,
                              int *qaprof, float **QArows );

boolean pkn_NRBComputeQTSQblf ( int n, int *prof, float *Amat, float **Arows,
                                int w, float *Bmat, float *bb,
                                int *qa11prof, float **QA11rows,
                                int *qa22prof, float **QA22rows,
                                float **QA21 );
boolean pkn_NRBComputeQSQTblf ( int n, int *prof, float *Amat, float **Arows,
                                int w, float *Bmat, float *bb,
                                int *qa11prof, float **QA11rows,
                                int *qa22prof, float **QA22rows,
                                float **QA21 );

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_PCGf ( int n, void *usrdata, float *b, float *x,
            boolean (*multAx)( int n, void *usrdata, const float *x, float *Ax ),
            boolean (*multQIx)( int n, void *usrdata, const float *x, float *Qix ),
            int maxit, float eps, float delta, int *itm );

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_MultSPMVectorf ( int nrows, int ncols, int nnz,
                             const index2 *ai, const float *ac,
                             int spdimen, const float *x,
                             float *y );
boolean pkn_MultSPMTVectorf ( int nrows, int ncols, int nnz,
                              const index2 *ai, const float *ac,
                              int spdimen, const float *x,
                              float *y );

boolean pkn_MultSPsubMVectorf ( int nrows, int ncols, int nnz,
                                const index3 *ai, const float *ac,
                                int spdimen, const float *x,
                                float *y );
boolean pkn_MultSPsubMTVectorf ( int nrows, int ncols, int nnz,
                                 const index3 *ai, const float *ac,
                                 int spdimen, const float *x,
                                 float *y );

/* fast sparse matrix multiplication - uses data produced by */
/* one of the pkn_SPMFind*nnz* or pkn_SPsubMFind*nnz* procedures */
void pkn_SPMFastMultMMf ( float *ac, float *bc,
                          int nnzab, int *abpos, index2 *aikbkj,
                          float *abc );

/* slower (but still decent) matrix multiplication procedures */
boolean pkn_SPMmultMMCf ( int nra, int nca, int ncb,
                          unsigned int nnza, index2 *ai, float *ac,
                          unsigned int *apermut, int *acols, boolean ca,
                          unsigned int nnzb, index2 *bi, float *bc,
                          unsigned int *bpermut, int *bcols, boolean cb,
                          index2 *abi, float *abc );
boolean pkn_SPMmultMMTCf ( int nra, int nca, int nrb,
                           unsigned int nnza, index2 *ai, float *ac,
                           unsigned int *apermut, int *acols, boolean ca,
                           unsigned int nnzb, index2 *bi, float *bc,
                           unsigned int *bpermut, int *brows, boolean rb,
                           index2 *abi, float *abc );
boolean pkn_SPMmultMTMCf ( int nra, int nca, int ncb,
                           unsigned int nnza, index2 *ai, float *ac,
                           unsigned int *apermut, int *arows, boolean ra,
                           unsigned int nnzb, index2 *bi, float *bc,
                           unsigned int *bpermut, int *bcols, boolean cb,
                           index2 *abi, float *abc );

boolean pkn_SPsubMmultMMCf ( int nra, int nca, int ncb,
                             unsigned int nnza, index3 *ai, float *ac,
                             unsigned int *apermut, int *acols, boolean ca,
                             unsigned int nnzb, index3 *bi, float *bc,
                             unsigned int *bpermut, int *bcols, boolean cb,
                             index2 *abi, float *abc );
boolean pkn_SPsubMmultMMTCf ( int nra, int nca, int nrb,
                              unsigned int nnza, index3 *ai, float *ac,
                              unsigned int *apermut, int *acols, boolean ca,
                              unsigned int nnzb, index3 *bi, float *bc,
                              unsigned int *bpermut, int *brows, boolean rb,
                              index2 *abi, float *abc );
boolean pkn_SPsubMmultMTMCf ( int nra, int nca, int ncb,
                              unsigned int nnza, index3 *ai, float *ac,
                              unsigned int *apermut, int *arows, boolean ra,
                              unsigned int nnzb, index3 *bi, float *bc,
                              unsigned int *bpermut, int *bcols, boolean cb,
                              index2 *abi, float *abc );

/* ///////////////////////////////////////////////////////////////////////// */
boolean pkn_QuadRectanglesf ( float a, float b, int n,
                              float *qknots, float *qcoeff );
boolean pkn_QuadSimpsonf ( float a, float b, int n,
                           float *qknots, float *qcoeff );

boolean pkn_QuadGaussLegendre4f ( float a, float b, int n,
                                  float *qknots, float *qcoeff );
boolean pkn_QuadGaussLegendre6f ( float a, float b, int n,
                                  float *qknots, float *qcoeff );
boolean pkn_QuadGaussLegendre8f ( float a, float b, int n,
                                  float *qknots, float *qcoeff );
boolean pkn_QuadGaussLegendre10f ( float a, float b, int n,
                                   float *qknots, float *qcoeff );
boolean pkn_QuadGaussLegendre12f ( float a, float b, int n,
                                   float *qknots, float *qcoeff );
boolean pkn_QuadGaussLegendre14f ( float a, float b, int n,
                                   float *qknots, float *qcoeff );
boolean pkn_QuadGaussLegendre16f ( float a, float b, int n,
                                   float *qknots, float *qcoeff );
boolean pkn_QuadGaussLegendre18f ( float a, float b, int n,
                                   float *qknots, float *qcoeff );
boolean pkn_QuadGaussLegendre20f ( float a, float b, int n,
                                   float *qknots, float *qcoeff );

/* ///////////////////////////////////////////////////////////////////////// */
#ifndef PKN_LMT_CONTINUE
/* return values of pkn_NLMIterf */
#define PKN_LMT_ERROR           -1
#define PKN_LMT_CONTINUE_LM      0
#define PKN_LMT_CONTINUE_LM_P    1
#define PKN_LMT_CONTINUE_N       2
#define PKN_LMT_FOUND_MINIMUM    3
#define PKN_LMT_FOUND_ZEROGRAD   4
#define PKN_LMT_FOUND_ZEROGRAD_P 5 
#define PKN_LMT_FOUND_BARRIER    6
#define PKN_LMT_CROSSED_LIMIT    7
#define PKN_LMT_NO_PROGRESS      8
#endif

#ifndef PKN_SD_CONTINUE
/* return values of pkn_SDIterd */
#define PKN_SD_ERROR            -1
#define PKN_SD_CONTINUE          0
#define PKN_SD_FOUND_ZEROGRAD    1
#define PKN_SD_FOUND_BARRIER     2
#define PKN_SD_CROSSED_LIMIT     3
#define PKN_SD_NO_PROGRESS       4
#endif

typedef boolean (*pkn_NLMTevalfuncf)( int n, void *usrdata,
                                      float *x, float *f );
typedef boolean (*pkn_NLMTevalfuncgf)( int n, void *usrdata,
                                       float *x, float *f, float *g );
typedef boolean (*pkn_NLMTevalfuncghf)( int n, void *usrdata,
                                        float *x, float *f, float *g, float *h );
typedef boolean (*pkn_NLMTtransfuncf)( int n, void *usrdata, float *x );
typedef boolean (*pkn_NLMTtunnelfuncf)( int n, void *usrdata,
                                        float *x0, float *x1, boolean *went_out );

int pkn_NLMIterf ( int n, void *usrdata, float *x,
                   pkn_NLMTevalfuncf funcf, pkn_NLMTevalfuncgf funcfg,
                   pkn_NLMTevalfuncghf funcfgh,
                   pkn_NLMTtransfuncf trans, pkn_NLMTtunnelfuncf tunnel,
                   float lowerbound, float eps, float delta,
                   float *nu );

int pkn_SDIterf ( int n, void *usrdata, float *x,
                  pkn_NLMTevalfuncf funcf, pkn_NLMTevalfuncgf funcfg,
                  pkn_NLMTtransfuncf trans, pkn_NLMTtunnelfuncf tunnel,
                  float lowerbound, float eps, float delta,
                  float *nu );

boolean _pkn_DivideIntervalf ( float *ga, float *gc, float *gd, float *gb,
                               float *fa, float *fc, float *fd, float *fb );

#ifdef __cplusplus
}
#endif

#endif /* PKNUMF_H*/

