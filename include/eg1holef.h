
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* this header file is for application use */

#ifndef EG1HOLEF_H
#define EG1HOLEF_H

#ifndef EGHOLEF_H
#include "egholef.h"
#endif

#ifdef __cplusplus   
extern "C" {
#endif

#define G1H_FINALDEG  5
#define G1H_OMCDEG 4  /*degree of common boundary curves */

/* the following constants are related with the spline basis */
#define G1H_S_MAX_NK  4
#define G1H_S_MAX_M1  2
#define G1H_S_MAX_M2  4

/* ///////////////////////////////////////////////////////////////////////// */
/* Below are the possible queries and answers for the option procedure. */
/* It is always safe to return the "default" answer with no data. */

#define G1H_DEFAULT                    0 /* default answer for all queries */

#define G1HQUERY_CENTRAL_POINT         1
#define G1H_CENTRAL_POINT_ALT          1
#define G1H_CENTRAL_POINT_GIVEN        2

#define G1HQUERY_CENTRAL_DERIVATIVES1  2
#define G1H_CENRTAL_DERIVATIVES1_ALT   1
#define G1H_CENTRAL_DERIVATIVES1_GIVEN 2

#define G1HQUERY_DOMAIN_CURVES         3

#define G1HQUERY_BASIS                 4
#define G1H_USE_RESTRICTED_BASIS       1

#define G1HQUERY_QUADRATURE            5
#define G1H_QUADRATURE_GAUSS_LEGENDRE  1

#define G1HQUERY_Q2_FORM_CONSTANT      6
#define G1H_Q2_USE_SUPPLIED_CONSTANT   1

/* ///////////////////////////////////////////////////////////////////////// */
/* core procedures */

void g1h_SetOptionProcf ( GHoleDomainf *domain,
    int (*OptionProc)( GHoleDomainf *domain, int query, int qn,
                       int *ndata, int **idata, float **fdata ) );

boolean g1h_ComputeBasisf ( GHoleDomainf *domain );

boolean g1h_ComputeFormMatrixf ( GHoleDomainf *domain );
boolean g1h_DecomposeMatrixf ( GHoleDomainf *domain );
boolean g1h_FillHolef ( GHoleDomainf *domain,
                        int spdimen, CONST_ float *hole_cp,
                        float *acoeff, void *usrptr,
                        void (*outpatch) ( int n, int m, const float *cp,
                                           void *usrptr ) );

boolean g1h_ComputeExtFormMatrixf ( GHoleDomainf *domain );
boolean g1h_DecomposeExtMatrixf ( GHoleDomainf *domain );
boolean g1h_ExtFillHolef ( GHoleDomainf *domain,
                           int spdimen, CONST_ float *hole_cp,
                           float *acoeff, void *usrptr,
                           void (*outpatch) ( int n, int m, const float *cp,
                                              void *usrptr ) );

int g1h_V0SpaceDimf ( GHoleDomainf *domain );
int g1h_ExtV0SpaceDimf ( GHoleDomainf *domain );
boolean g1h_GetBPDerivativesf ( GHoleDomainf *domain, int cno, float *val );
boolean g1h_GetBFuncPatchf ( GHoleDomainf *domain, int fn, int pn, float *bp );

boolean g1h_SetConstraintMatrixf ( GHoleDomainf *domain,
                                   int nconstr, const float *cmat );
boolean g1h_FillHoleConstrf ( GHoleDomainf *domain,
                              int spdimen, CONST_ float *hole_cp,
                              int nconstr, CONST_ float *constr,
                              float *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const float *cp,
                                                 void *usrptr ) );

boolean g1h_SetAltConstraintMatrixf ( GHoleDomainf *domain, int spdimen,
                                      int nconstr, const float *cmat );
boolean g1h_FillHoleAltConstrf ( GHoleDomainf *domain,
                              int spdimen, CONST_ float *hole_cp,
                              int naconstr, CONST_ float *constr,
                              float *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const float *cp,
                                                 void *usrptr ) );

boolean g1h_SetExtConstraintMatrixf ( GHoleDomainf *domain,
                                      int nconstr, const float *cmat );
boolean g1h_ExtFillHoleConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int nconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) );

boolean g1h_SetExtAltConstraintMatrixf ( GHoleDomainf *domain, int spdimen,
                                      int naconstr, const float *acmat );
boolean g1h_ExtFillHoleAltConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int naconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) );

float g1h_FunctionalValuef ( GHoleDomainf *domain, int spdimen,
                             const float *hole_cp, const float *acoeff );
float g1h_ExtFunctionalValuef ( GHoleDomainf *domain, int spdimen,
                                CONST_ float *hole_cp, CONST_ float *acoeff );
boolean g1h_NLFunctionalValuef ( GHoleDomainf *domain,
                                 const point3f *hole_cp, const vector3f *acoeff,
                                 float *funcval );
boolean g1h_NLExtFunctionalValuef ( GHoleDomainf *domain,
                                    const point3f *hole_cp, const vector3f *acoeff,
                                    float *funcval );

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_GetFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                   CONST_ float *hole_cp, float *acoeff,
                                   void (*outcurve) ( int n, const float *cp ) );
boolean g1h_GetExtFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                      CONST_ float *hole_cp, float *acoeff,
                                      void (*outcurve) ( int n, const float *cp ) );
boolean g1h_GetSplFinalPatchCurvesf ( GHoleDomainf *domain, int spdimen,
                                      CONST_ float *hole_cp, float *acoeff,
                                      void (*outcurve) ( int n, int lkn,
                                               const float *kn, const float *cp ) );

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_ComputeNLNormalf ( GHoleDomainf *domain, 
                               const point3f *hole_cp,
                               vector3f *anv );

boolean g1h_NLFillHolef ( GHoleDomainf *domain, const point3f *hole_cp,
                          float *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const point3f *cp,
                                             void *usrptr ) );
boolean g1h_NLFillHoleConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                    int nconstr, const vector3f *constr,
                    float *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3f *cp,
                                       void *usrptr ) );
boolean g1h_NLFillHoleAltConstrf ( GHoleDomainf *domain, const point3f *hole_cp,   
                    int nconstr, const float *constr,
                    float *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3f *cp,
                                       void *usrptr ) );

boolean g1h_NLExtFillHolef ( GHoleDomainf *domain, const point3f *hole_cp,
                             float *acoeff, void *usrptr,
                             void (*outpatch) ( int n, int m, const point3f *cp,
                                                void *usrptr ) );
boolean g1h_NLExtFillHoleConstrf ( GHoleDomainf *domain,
                     const point3f *hole_cp,
                     int nconstr, const vector3f *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int m, const point3f *cp,
                                        void *usrptr ) );
boolean g1h_NLExtFillHoleAltConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                         int naconstr, const float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const point3f *cp,
                                            void *usrptr ) );

boolean g1h_NLSplFillHolef ( GHoleDomainf *domain, const point3f *hole_cp,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) );
boolean g1h_NLSplFillHoleConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                     int nconstr, const vector3f *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) );
boolean g1h_NLSplFillHoleAltConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                     int nconstr, const float *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* G1, quasi G2 hole filling procedures */
void g1h_DestroyQ2PrivateDataf ( GHoleDomainf *domain );

boolean g1h_Q2ComputeFormMatrixf ( GHoleDomainf *domain );
boolean g1h_Q2DecomposeMatrixf ( GHoleDomainf *domain );
boolean g1h_Q2FillHolef ( GHoleDomainf *domain,
                          int spdimen, CONST_ float *hole_cp,
                          float *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const float *cp,
                                             void *usrptr ) );
boolean g1h_Q2FillHoleConstrf ( GHoleDomainf *domain,
                                int spdimen, CONST_ float *hole_cp,
                                int nconstr, CONST_ float *constr,
                                float *acoeff, void *usrptr,
                                void (*outpatch) ( int n, int m, const float *cp,
                                                   void *usrptr ) );
boolean g1h_Q2FillHoleAltConstrf ( GHoleDomainf *domain,
                              int spdimen, CONST_ float *hole_cp,
                              int naconstr, CONST_ float *constr,
                              float *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const float *cp,
                                                 void *usrptr ) );

boolean g1h_Q2NLFillHoleConstrf ( GHoleDomainf *domain, CONST_ point3f *hole_cp,
                                  int nconstr, CONST_ vector3f *constr,
                                  float *acoeff, void *usrptr,
                                  void (*outpatch) ( int n, int m,
                                                     const point3f *cp,
                                                     void *usrptr ) );

boolean g1h_Q2NLFillHoleAltConstrf ( GHoleDomainf *domain, CONST_ point3f *hole_cp,
                                     int naconstr, CONST_ float *constr,
                                     float *acoeff, void *usrptr,
                                     void (*outpatch) ( int n, int m,
                                                        const point3f *cp,
                                                        void *usrptr ) );

boolean g1h_Q2ExtComputeFormMatrixf ( GHoleDomainf *domain );
boolean g1h_Q2ExtDecomposeMatrixf ( GHoleDomainf *domain );
boolean g1h_Q2ExtFillHolef ( GHoleDomainf *domain,
                             int spdimen, CONST_ float *hole_cp,
                             float *acoeff, void *usrptr,
                             void (*outpatch) ( int n, int m, const float *cp,
                                                void *usrptr ) );

boolean g1h_Q2ExtFillHoleConstrf ( GHoleDomainf *domain,
                         int spdimen, CONST_ float *hole_cp,
                         int nconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) );
boolean g1h_Q2ExtFillHoleAltConstrf ( GHoleDomainf *domain, 
                         int spdimen, CONST_ float *hole_cp,
                         int naconstr, CONST_ float *constr,
                         float *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const float *cp,
                                            void *usrptr ) );

boolean g1h_Q2NLExtFillHoleConstrf ( GHoleDomainf *domain, CONST_ point3f *hole_cp,
                                  int nconstr, CONST_ vector3f *constr,
                                  float *acoeff, void *usrptr,
                                  void (*outpatch) ( int n, int m,
                                                     const point3f *cp,
                                                     void *usrptr ) );

boolean g1h_Q2NLExtFillHoleAltConstrf ( GHoleDomainf *domain, CONST_ point3f *hole_cp,
                                     int naconstr, CONST_ float *constr,
                                     float *acoeff, void *usrptr,
                                     void (*outpatch) ( int n, int m,
                                                        const point3f *cp,
                                                        void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* spline basis procedures */
boolean g1h_ComputeSplBasisf ( GHoleDomainf *domain, int nk, int m1, int m2 );
boolean g1h_ComputeSplFormMatrixf ( GHoleDomainf *domain );
boolean g1h_DecomposeSplMatrixf ( GHoleDomainf *domain );
boolean g1h_SplFillHolef ( GHoleDomainf *domain,
               int spdimen, CONST_ float *hole_cp,
               float *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const float *knu,
                                  int m, int lknv, const float *knv,
                                  const float *cp, void *usrptr ) );

int g1h_SplV0SpaceDimf ( GHoleDomainf *domain );

boolean g1h_Q2SplComputeFormMatrixf ( GHoleDomainf *domain );
boolean g1h_Q2SplDecomposeMatrixf ( GHoleDomainf *domain );
boolean g1h_Q2SplFillHolef ( GHoleDomainf *domain,
               int spdimen, CONST_ float *hole_cp,
               float *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const float *knu,
                                  int m, int lknv, const float *knv,
                                  const float *cp, void *usrptr ) );

boolean g1h_SetSplConstraintMatrixf ( GHoleDomainf *domain,
                                      int nconstr, const float *cmat );
boolean g1h_SplFillHoleConstrf ( GHoleDomainf *domain,
               int spdimen, CONST_ float *hole_cp,
               int nconstr, CONST_ float *constr,
               float *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const float *knu,
                                  int m, int lknv, const float *knv,
                                  const float *cp, void *usrptr ) );
boolean g1h_Q2SplFillHoleConstrf ( GHoleDomainf *domain,
               int spdimen, CONST_ float *hole_cp,
               int nconstr, CONST_ float *constr,
               float *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const float *knu,
                                  int m, int lknv, const float *knv,
                                  const float *cp, void *usrptr ) );

boolean g1h_SetSplAltConstraintMatrixf ( GHoleDomainf *domain, int spdimen,
                                         int naconstr, const float *acmat );
boolean g1h_SplFillHoleAltConstrf ( GHoleDomainf *domain,
               int spdimen, CONST_ float *hole_cp,
               int naconstr, CONST_ float *constr,
               float *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const float *knu,
                                  int m, int lknv, const float *knv,
                                  const float *cp, void *usrptr ) );
boolean g1h_Q2SplFillHoleAltConstrf ( GHoleDomainf *domain,
               int spdimen, CONST_ float *hole_cp,
               int naconstr, CONST_ float *constr,
               float *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const float *knu,
                                  int m, int lknv, const float *knv,
                                  const float *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* G1 quasi G2 nonlinear construction procedures */
boolean g1h_Q2NLFillHolef ( GHoleDomainf *domain,
                    const point3f *hole_cp,
                    float *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3f *cp,
                                       void *usrptr ) );

boolean g1h_Q2NLExtFillHolef ( GHoleDomainf *domain,
                    const point3f *hole_cp,
                    float *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3f *cp,
                                       void *usrptr ) );

boolean g1h_Q2NLSplFillHolef ( GHoleDomainf *domain,
                    const point3f *hole_cp,
                    float *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int lknu, const float *knu,
                                       int m, int lknv, const float *knv,
                                       const point3f *cp, void *usrptr ) );
boolean g1h_Q2NLSplFillHoleConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                     int nconstr, const vector3f *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) );
boolean g1h_Q2NLSplFillHoleAltConstrf ( GHoleDomainf *domain, const point3f *hole_cp,
                     int nconstr, const float *constr,
                     float *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const float *knu,
                                        int m, int lknv, const float *knv,
                                        const point3f *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* drawing procedures */
void g1h_DrawDomAuxPatchesf ( GHoleDomainf *domain,
               void (*drawpatch) ( int n, int m, const point2f *cp ) );
void g1h_DrawBasAuxPatchesf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int m, const float *cp ) );
boolean g1h_DrawJFunctionf ( GHoleDomainf *domain, int k, int l,
                             void (*drawpoly) ( int deg, const float *f ) );
void g1h_DrawDiPatchesf ( GHoleDomainf *domain,
                      void (*drawpatch) ( int n, int m, const point2f *cp ) );
void g1h_ExtractPartitionf ( GHoleDomainf *domain,
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
void g1h_ExtractCentralPointf ( GHoleDomainf *domain,
                                point2f *centp, vector2f *centder );
void g1h_DrawBasAFunctionf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3f *cp ) );
void g1h_DrawBasBFunctionf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3f *cp ) );
boolean g1h_DrawBasCNetf ( GHoleDomainf *domain, int fn,
               void (*drawnet) ( int n, int m, const point3f *cp ) );
void g1h_DrawBFAomcf ( GHoleDomainf *domain, int fn,
                       void (*drawpoly)(int degree, const float *coeff) );
void g1h_DrawBFBomcf ( GHoleDomainf *domain, int fn,
                       void (*drawpoly)(int degree, const float *coeff) );
void g1h_DrawFinalSurfBCf ( GHoleDomainf *domain,
                            int spdimen, const float *hole_cp,
                            const float *acoeff,
                            void (*drawcurve)(int degree, int spdimen,
                                              const float *cp) );
void g1h_ExtDrawFinalSurfBCf ( GHoleDomainf *domain,
                               int spdimen, const float *hole_cp,
                               const float *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                                 const float *cp) );
void g1h_DrawMatricesf ( GHoleDomainf *domain,
                         void (*drawmatrix)(int nfa, int nfb,
                                            float *amat, float *bmat) );
void g1h_DrawExtMatricesf ( GHoleDomainf *domain,
                            void (*drawmatrix)(int k, int r, int s, float *Aii,
                                               float *Bi) );
int g1h_DrawBFcpnf ( int hole_k, unsigned char *bfcpn );

void g1h_Q2DrawMatricesf ( GHoleDomainf *domain,
                           void (*drawmatrix)(int nfa, int nfb,
                                              float *amat, float *bmat) );
void g1h_Q2DrawExtMatricesf ( GHoleDomainf *domain,
                              void (*drawmatrix)(int k, int r, int s,
                                                 float *Aii, float *Bi) );
boolean g1h_GetABasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                    float *bfpc );
boolean g1h_GetBBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                    float *bfpc );

/* ///////////////////////////////////////////////////////////////////////// */
/* spline basis drawing procedures */
void g1h_DrawSplBasFuncNumf ( GHoleDomainf *domain,
                        int *nfunc_a, int *nfunc_b, int *nfunc_c, int *nfunc_d );

void g1h_DrawSplBasAuxPatchesf ( GHoleDomainf *domain, int fn,
               void (*drawpatch) ( int n, int lknu, const float *knu,
                                   int m, int lknv, const float *knv,
                                   const point3f *cp ) );

void g1h_DrawSplBasFunctionf ( GHoleDomainf *domain, int fn,
             void (*drawpatch) ( int n, int lknu, const float *knu,
                                 int m, int lknv, const float *knv,
                                 const point3f *cp ) );

void g1h_DrawSplBFAomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) );
void g1h_DrawSplBFBomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) );
void g1h_DrawSplBFDomcf ( GHoleDomainf *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const float *knots,
                                           const float *coeff) );
void g1h_DrawSplFinalSurfBCf ( GHoleDomainf *domain,
                               int spdimen, const float *hole_cp,
                               const float *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                             int lastknot, const float *knots,
                                             const float *cp) );

void g1h_DrawSplMatricesf ( GHoleDomainf *domain,
                            void (*drawmatrix)( int k, int r, int s, int t,
                                                float *A, float *B ) );

void g1h_Q2DrawSplMatricesf ( GHoleDomainf *domain,
                              void (*drawmatrix)(int k, int r, int s,
                                                 float *Aii, float *Bi) );

boolean g1h_GetSplABasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc );
boolean g1h_GetSplBBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc );
boolean g1h_GetSplDBasisFPatchCurvef ( GHoleDomainf *domain, int fn, int i,
                                       int *lkn, float *kn, float *bfpc );

/* ///////////////////////////////////////////////////////////////////////// */
int g1h_SymPatchMatrixSize ( int hole_k );
boolean g1h_GetSymPatchMatrixf ( GHoleDomainf *domain, float *patchmatrix );
boolean g1h_GetExtSymPatchMatrixf ( GHoleDomainf *domain, float *patchmatrix );
boolean g1h_Q2GetSymPatchMatrixf ( GHoleDomainf *domain, float *patchmatrix );
boolean g1h_Q2GetExtSymPatchMatrixf ( GHoleDomainf *domain, float *patchmatrix );
boolean g1h_MatrixFillSymHolef ( int hole_k, const float *patchmatrix,
                                 int spdimen, const float *hole_cp, void *usrptr,
                                 void (*outpatch) ( int n, int m, const float *cp,
                                                    void *usrptr ) );
/* ///////////////////////////////////////////////////////////////////////// */
int g1h_GetErrorCodef ( GHoleDomainf *domain, char **ErrorString );

#ifdef __cplusplus
}
#endif

#endif

