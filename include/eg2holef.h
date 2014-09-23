
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* this header file is for application use */

#ifndef EG2HOLEF_H
#define EG2HOLEF_H

#ifndef EGHOLEF_H
#include "egholef.h"
#endif

#ifdef __cplusplus   
extern "C" {
#endif

#ifndef CONST_  /* a dirty trick to suppress many compiler warning messages */
#define CONST_ const
#endif

#define G2H_FINALDEG  9

#define G2H_OMCDEG 7  /* degree of common boundary curves */

/* the following constants are related with the spline basis */
#define G2H_S_MAX_NK  4
#define G2H_S_MAX_M1  3
#define G2H_S_MAX_M2  7

/* ///////////////////////////////////////////////////////////////////////// */
/* Below are the possible queries and answers for the option procedure. */
/* It is always safe to return the "default" answer with no data. */

#define G2H_DEFAULT                    0 /* default answer for all queries */

#define G2HQUERY_CENTRAL_POINT         1
#define G2H_CENTRAL_POINT_ALT          1
#define G2H_CENTRAL_POINT_GIVEN        2

#define G2HQUERY_CENTRAL_DERIVATIVES1  2
#define G2H_CENTRAL_DERIVATIVES1_ALT   1
#define G2H_CENTRAL_DERIVATIVES1_GIVEN 2

#define G2HQUERY_DOMAIN_CURVES         3
#define G2H_DOMAIN_CURVES_DEG4         1

#define G2HQUERY_BASIS                 4
#define G2H_USE_RESTRICTED_BASIS       1

#define G2HQUERY_QUADRATURE            5
#define G2H_QUADRATURE_GAUSS_LEGENDRE  1

/* ///////////////////////////////////////////////////////////////////////// */
/* core procedures */

void g2h_SetOptionProcf ( GHoleDomainf *domain,
    int (*OptionProc)( GHoleDomainf *domain, int query, int qn,
                       int *ndata, int **idata, float **fdata ) );

boolean g2h_ComputeBasisf ( GHoleDomainf *domain );

boolean g2h_ComputeFormMatrixf ( GHoleDomainf *domain );
boolean g2h_DecomposeMatrixf ( GHoleDomainf *domain );
boolean g2h_FillHolef ( GHoleDomainf *domain,
                        int spdimen, CONST_ float *hole_cp,
                        float *acoeff, void *usrptr,
                        void (*outpatch) ( int n, int m, const float *cp,
                                           void *usrptr ) );

boolean g2h_ComputeExtFormMatrixf ( GHoleDomainf *domain );
boolean g2h_DecomposeExtMatrixf ( GHoleDomainf *domain );
boolean g2h_ExtFillHolef ( GHoleDomainf *domain,
                           int spdimen, CONST_ float *hole_cp,
                           float *acoeff, void *usrptr,
                           void (*outpatch) ( int n, int m, const float *cp,
                                              void *usrptr ) );

int g2h_V0SpaceDimf ( GHoleDomainf *domain );
int g2h_ExtV0SpaceDimf ( GHoleDomainf *domain );
boolean g2h_GetBPDerivativesf ( GHoleDomainf *domain, int cno, float *val );
boolean g2h_GetBFuncPatchf ( GHoleDomainf *domain, int fn, int pn, float *bp );

boolean g2h_SetConstraintMatrixf ( GHoleDomainf *domain,
                                   int nconstr, const float *cmat );
boolean g2h_FillHoleConstrf ( GHoleDomainf *domain,
                              int spdimen, CONST_ float *hole_cp,
                              int nconstr, CONST_ float *constr,
                              float *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const float *cp,
                                                 void *usrptr ) );

boolean g2h_SetAltConstraintMatrixf ( GHoleDomainf *domain, int spdimen,
                                      int nconstr, const float *cmat );
boolean g2h_FillHoleAltConstrf ( GHoleDomainf *domain,
                              int spdimen, CONST_ float *hole_cp,
                              int naconstr, CONST_ float *constr,
                              float *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const float *cp,
                                                 void *usrptr ) );

boolean g2h_SetExtConstraintMatrixf ( GHoleDomainf *domain,
                                      int nconstr, const float *cmat );
boolean g2h_ExtFillHoleConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int nconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) );

boolean g2h_SetExtAltConstraintMatrixf ( GHoleDomainf *domain, int spdimen,
                                      int naconstr, const float *acmat );
boolean g2h_ExtFillHoleAltConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int naconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) );

float g2h_FunctionalValuef ( GHoleDomainf *domain, int spdimen,
                             CONST_ float *hole_cp, CONST_ float *acoeff );
float g2h_ExtFunctionalValuef ( GHoleDomainf *domain, int spdimen,
                                CONST_ float *hole_cp, CONST_ float *acoeff );
boolean g2h_NLFunctionalValuef ( GHoleDomainf *domain,   
                                 const point3f *hole_cp, const vector3f *acoeff,
                                 float *funcval );
boolean g2h_NLExtFunctionalValuef ( GHoleDomainf *domain,
                                    const point3f *hole_cp, const vector3f *acoeff,
                                    float *funcval );

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_ComputeNLNormalf ( GHoleDomainf *domain, 
                               const point3f *hole_cp,
                               vector3f *anv );

boolean g2h_NLFillHolef ( GHoleDomainf *domain, const point3f *hole_cp,
                          float *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const point3f *cp,
                                             void *usrptr ) );
boolean g2h_NLFillHoleConstrf ( GHoleDomainf *domain, const point3f *hole_cp,   
                    int nconstr, const vector3f *constr,
                    float *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3f *cp,
                                       void *usrptr ) );
boolean g2h_NLFillHoleAltConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                    int nconstr, const float *constr,
                    float *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3f *cp,
                                       void *usrptr ) );

boolean g2h_NLExtFillHolef ( GHoleDomainf *domain, const point3f *hole_cp,
                             float *acoeff, void *usrptr,
                             void (*outpatch) ( int n, int m, const point3f *cp,
                                                void *usrptr ) );
boolean g2h_NLExtFillHoleConstrf ( GHoleDomainf *domain,
                     const point3f *hole_cp,
                     int nconstr, const vector3f *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int m, const point3f *cp,
                                        void *usrptr ) );
boolean g2h_NLExtFillHoleAltConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                         int naconstr, const float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const point3f *cp,
                                            void *usrptr ) ); 

boolean g2h_NLSplFillHolef ( GHoleDomainf *domain, const point3f *hole_cp,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) );
boolean g2h_NLSplFillHoleConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                     int nconstr, const vector3f *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) );
boolean g2h_NLSplFillHoleAltConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                     int nconstr, const float *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_GetFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                   CONST_ float *hole_cp, float *acoeff,
                                   void (*outcurve) ( int n, const float *cp ) );
boolean g2h_GetExtFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                      CONST_ float *hole_cp, float *acoeff,
                                      void (*outcurve) ( int n, const float *cp ) );
boolean g2h_GetSplFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                      CONST_ float *hole_cp, float *acoeff,
                                      void (*outcurve) ( int n, int lkn,
                                               const float *kn, const float *cp ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* spline basis procedures */
boolean g2h_ComputeSplBasisf ( GHoleDomainf *domain, int nk, int m1, int m2 );
boolean g2h_ComputeSplFormMatrixf ( GHoleDomainf *domain );
boolean g2h_DecomposeSplMatrixf ( GHoleDomainf *domain );
boolean g2h_SplFillHolef ( GHoleDomainf *domain,
               int spdimen, CONST_ float *hole_cp,
               float *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const float *knu,
                                  int m, int lknv, const float *knv,            
                                  const float *cp, void *usrptr ) );

int g2h_SplV0SpaceDimf ( GHoleDomainf *domain );
boolean g2h_SetSplConstraintMatrixf ( GHoleDomainf *domain,
                                      int nconstr, const float *cmat );
boolean g2h_SplFillHoleConstrf ( GHoleDomainf *domain,
               int spdimen, CONST_ float *hole_cp,
               int nconstr, CONST_ float *constr,
               float *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const float *knu,
                                  int m, int lknv, const float *knv,
                                  const float *cp, void *usrptr ) );

boolean g2h_SetSplAltConstraintMatrixf ( GHoleDomainf *domain, int spdimen,
                                         int naconstr, const float *acmat );
boolean g2h_SplFillHoleAltConstrf ( GHoleDomainf *domain,
               int spdimen, CONST_ float *hole_cp,
               int naconstr, CONST_ float *constr,
               float *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const float *knu,
                                  int m, int lknv, const float *knv,
                                  const float *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* drawing procedures */
void g2h_DrawDomAuxPatchesf ( GHoleDomainf *domain,
               void (*drawpatch) ( int n, int m, const point2f *cp ) );
void g2h_DrawBasAuxPatchesf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int m, const float *cp ) );
boolean g2h_DrawJFunctionf ( GHoleDomainf *domain, int k, int l,
                             void (*drawpoly) ( int deg, const float *f ) );
void g2h_DrawDiPatchesf ( GHoleDomainf *domain,
                      void (*drawpatch) ( int n, int m, const point2f *cp ) );
void g2h_ExtractPartitionf ( GHoleDomainf *domain,
                             int *hole_k, int *hole_m,
                             float *partition,
                             float *part_delta,
                             float *spart_alpha,
                             float *spart_malpha,
                             float *spart_salpha,
                             float *spart_knot,
                             float *alpha0,
                             boolean *spart_sgn,
                             boolean *spart_both );
void g2h_ExtractCentralPointf ( GHoleDomainf *domain,
                                point2f *centp, vector2f *centder );
void g2h_DrawBasAFunctionf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3f *cp ) );
void g2h_DrawBasBFunctionf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3f *cp ) );
boolean g2h_DrawBasCNetf ( GHoleDomainf *domain, int fn,
               void (*drawnet) ( int n, int m, const point3f *cp ) );
void g2h_DrawBFAomcf ( GHoleDomainf *domain, int fn, 
                       void (*drawpoly)(int degree, const float *coeff) );
void g2h_DrawBFBomcf ( GHoleDomainf *domain, int fn,
                       void (*drawpoly)(int degree, const float *coeff) );
void g2h_DrawFinalSurfBCf ( GHoleDomainf *domain,   
                            int spdimen, const float *hole_cp,
                            const float *acoeff, 
                            void (*drawcurve)(int degree, int spdimen,
                                              const float *cp) );
void g2h_ExtDrawFinalSurfBCf ( GHoleDomainf *domain,
                               int spdimen, const float *hole_cp,
                               const float *acoeff, 
                               void (*drawcurve)(int degree, int spdimen,
                                                 const float *cp) );
void g2h_DrawMatricesf ( GHoleDomainf *domain,
                         void (*drawmatrix)(int nfa, int nfb,
                                            float *amat, float *bmat) );
void g2h_DrawExtMatricesf ( GHoleDomainf *domain,
                            void (*drawmatrix)(int k, int r, int s,
                                               float *Aii, float *Bi) );
int g2h_DrawBFcpnf ( int hole_k, unsigned char *bfcpn );

boolean g2h_GetABasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                    float *bfpc );
boolean g2h_GetBBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                    float *bfpc );

/* ///////////////////////////////////////////////////////////////////////// */
/* spline basis drawing procedures */
void g2h_DrawSplBasFuncNumf ( GHoleDomainf *domain,
                        int *nfunc_a, int *nfunc_b, int *nfunc_c, int *nfunc_d );

void g2h_DrawSplBasAuxPatchesf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int lknu, const float *knu,
                                   int m, int lknv, const float *knv,
                                   const point3f *cp ) );

void g2h_DrawSplBasFunctionf ( GHoleDomainf *domain, int fn,
             void (*drawpatch) ( int n, int lknu, const float *knu,
                                 int m, int lknv, const float *knv,
                                 const point3f *cp ) );

void g2h_DrawSplBFAomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) );
void g2h_DrawSplBFBomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) );
void g2h_DrawSplBFDomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) );
void g2h_DrawSplFinalSurfBCf ( GHoleDomainf *domain,
                               int spdimen, const float *hole_cp,
                               const float *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                             int lastknot, const float *knots,
                                             const float *cp) );

void g2h_DrawSplMatricesf ( GHoleDomainf *domain,
                            void (*drawmatrix)( int k, int r, int s, int t,
                                                float *A, float *B ) );

boolean g2h_GetSplABasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc );
boolean g2h_GetSplBBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc );
boolean g2h_GetSplDBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc );

/* ///////////////////////////////////////////////////////////////////////// */
int g2h_SymPatchMatrixSize ( int hole_k );
boolean g2h_GetSymPatchMatrixf ( GHoleDomainf *domain, float *patchmatrix );
boolean g2h_GetExtSymPatchMatrixf ( GHoleDomainf *domain, float *patchmatrix );
boolean g2h_MatrixFillSymHolef ( int hole_k, const float *patchmatrix,
                                 int spdimen, const float *hole_cp, void *usrptr,
                                 void (*outpatch) ( int n, int m, const float *cp,
                                                    void *usrptr ) );
/* ///////////////////////////////////////////////////////////////////////// */
int g2h_GetErrorCodef ( GHoleDomainf *domain, char **ErrorString );

#ifdef __cplusplus
}
#endif

#endif

