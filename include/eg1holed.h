
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* this header file is for application use */

#ifndef EG1HOLED_H
#define EG1HOLED_H

#ifndef EGHOLED_H
#include "egholed.h"
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

void g1h_SetOptionProcd ( GHoleDomaind *domain,
    int (*OptionProc)( GHoleDomaind *domain, int query, int qn,
                       int *ndata, int **idata, double **fdata ) );

boolean g1h_ComputeBasisd ( GHoleDomaind *domain );

boolean g1h_ComputeFormMatrixd ( GHoleDomaind *domain );
boolean g1h_DecomposeMatrixd ( GHoleDomaind *domain );
boolean g1h_FillHoled ( GHoleDomaind *domain,
                        int spdimen, CONST_ double *hole_cp,
                        double *acoeff, void *usrptr,
                        void (*outpatch) ( int n, int m, const double *cp,
                                           void *usrptr ) );

boolean g1h_ComputeExtFormMatrixd ( GHoleDomaind *domain );
boolean g1h_DecomposeExtMatrixd ( GHoleDomaind *domain );
boolean g1h_ExtFillHoled ( GHoleDomaind *domain,
                           int spdimen, CONST_ double *hole_cp,
                           double *acoeff, void *usrptr,
                           void (*outpatch) ( int n, int m, const double *cp,
                                              void *usrptr ) );

int g1h_V0SpaceDimd ( GHoleDomaind *domain );
int g1h_ExtV0SpaceDimd ( GHoleDomaind *domain );
boolean g1h_GetBPDerivativesd ( GHoleDomaind *domain, int cno, double *val );
boolean g1h_GetBFuncPatchd ( GHoleDomaind *domain, int fn, int pn, double *bp );

boolean g1h_SetConstraintMatrixd ( GHoleDomaind *domain,
                                   int nconstr, const double *cmat );
boolean g1h_FillHoleConstrd ( GHoleDomaind *domain,
                              int spdimen, CONST_ double *hole_cp,
                              int nconstr, CONST_ double *constr,
                              double *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const double *cp,
                                                 void *usrptr ) );

boolean g1h_SetAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                      int nconstr, const double *cmat );
boolean g1h_FillHoleAltConstrd ( GHoleDomaind *domain,
                              int spdimen, CONST_ double *hole_cp,
                              int naconstr, CONST_ double *constr,
                              double *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const double *cp,
                                                 void *usrptr ) );

boolean g1h_SetExtConstraintMatrixd ( GHoleDomaind *domain,
                                      int nconstr, const double *cmat );
boolean g1h_ExtFillHoleConstrd ( GHoleDomaind *domain,
                         int spdimen, CONST_ double *hole_cp,
                         int nconstr, CONST_ double *constr,
                         double *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const double *cp,
                                            void *usrptr ) );

boolean g1h_SetExtAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                      int naconstr, const double *acmat );
boolean g1h_ExtFillHoleAltConstrd ( GHoleDomaind *domain,
                         int spdimen, CONST_ double *hole_cp,
                         int naconstr, CONST_ double *constr,
                         double *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const double *cp,
                                            void *usrptr ) );

double g1h_FunctionalValued ( GHoleDomaind *domain, int spdimen,
                             const double *hole_cp, const double *acoeff );
double g1h_ExtFunctionalValued ( GHoleDomaind *domain, int spdimen,
                                CONST_ double *hole_cp, CONST_ double *acoeff );
boolean g1h_NLFunctionalValued ( GHoleDomaind *domain,
                                 const point3d *hole_cp, const vector3d *acoeff,
                                 double *funcval );
boolean g1h_NLExtFunctionalValued ( GHoleDomaind *domain,
                                    const point3d *hole_cp, const vector3d *acoeff,
                                    double *funcval );

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_GetFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                   CONST_ double *hole_cp, double *acoeff,
                                   void (*outcurve) ( int n, const double *cp ) );
boolean g1h_GetExtFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                      CONST_ double *hole_cp, double *acoeff,
                                      void (*outcurve) ( int n, const double *cp ) );
boolean g1h_GetSplFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                      CONST_ double *hole_cp, double *acoeff,
                                      void (*outcurve) ( int n, int lkn,
                                               const double *kn, const double *cp ) );

/* ///////////////////////////////////////////////////////////////////////// */
boolean g1h_ComputeNLNormald ( GHoleDomaind *domain, 
                               const point3d *hole_cp,
                               vector3d *anv );

boolean g1h_NLFillHoled ( GHoleDomaind *domain, const point3d *hole_cp,
                          double *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const point3d *cp,
                                             void *usrptr ) );
boolean g1h_NLFillHoleConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                    int nconstr, const vector3d *constr,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) );
boolean g1h_NLFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,   
                    int nconstr, const double *constr,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) );

boolean g1h_NLExtFillHoled ( GHoleDomaind *domain, const point3d *hole_cp,
                             double *acoeff, void *usrptr,
                             void (*outpatch) ( int n, int m, const point3d *cp,
                                                void *usrptr ) );
boolean g1h_NLExtFillHoleConstrd ( GHoleDomaind *domain,
                     const point3d *hole_cp,
                     int nconstr, const vector3d *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int m, const point3d *cp,
                                        void *usrptr ) );
boolean g1h_NLExtFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                         int naconstr, const double *constr,
                         double *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const point3d *cp,
                                            void *usrptr ) );

boolean g1h_NLSplFillHoled ( GHoleDomaind *domain, const point3d *hole_cp,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const double *knu,
                                        int m, int lknv, const double *knv,
                                        const point3d *cp, void *usrptr ) );
boolean g1h_NLSplFillHoleConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                     int nconstr, const vector3d *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const double *knu,
                                        int m, int lknv, const double *knv,
                                        const point3d *cp, void *usrptr ) );
boolean g1h_NLSplFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                     int nconstr, const double *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const double *knu,
                                        int m, int lknv, const double *knv,
                                        const point3d *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* G1, quasi G2 hole filling procedures */
void g1h_DestroyQ2PrivateDatad ( GHoleDomaind *domain );

boolean g1h_Q2ComputeFormMatrixd ( GHoleDomaind *domain );
boolean g1h_Q2DecomposeMatrixd ( GHoleDomaind *domain );
boolean g1h_Q2FillHoled ( GHoleDomaind *domain,
                          int spdimen, CONST_ double *hole_cp,
                          double *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const double *cp,
                                             void *usrptr ) );
boolean g1h_Q2FillHoleConstrd ( GHoleDomaind *domain,
                                int spdimen, CONST_ double *hole_cp,
                                int nconstr, CONST_ double *constr, 
                                double *acoeff, void *usrptr,
                                void (*outpatch) ( int n, int m, const double *cp,
                                                   void *usrptr ) );
boolean g1h_Q2FillHoleAltConstrd ( GHoleDomaind *domain,
                              int spdimen, CONST_ double *hole_cp,
                              int naconstr, CONST_ double *constr,
                              double *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const double *cp,
                                                 void *usrptr ) );

boolean g1h_Q2NLFillHoleConstrd ( GHoleDomaind *domain, CONST_ point3d *hole_cp,
                                  int nconstr, CONST_ vector3d *constr,
                                  double *acoeff, void *usrptr,
                                  void (*outpatch) ( int n, int m,
                                                     const point3d *cp,
                                                     void *usrptr ) );

boolean g1h_Q2NLFillHoleAltConstrd ( GHoleDomaind *domain, CONST_ point3d *hole_cp,
                                     int naconstr, CONST_ double *constr,
                                     double *acoeff, void *usrptr,
                                     void (*outpatch) ( int n, int m,
                                                        const point3d *cp,
                                                        void *usrptr ) );

boolean g1h_Q2ExtComputeFormMatrixd ( GHoleDomaind *domain );
boolean g1h_Q2ExtDecomposeMatrixd ( GHoleDomaind *domain );
boolean g1h_Q2ExtFillHoled ( GHoleDomaind *domain,
                             int spdimen, CONST_ double *hole_cp,
                             double *acoeff, void *usrptr,
                             void (*outpatch) ( int n, int m, const double *cp,
                                                void *usrptr ) );

boolean g1h_Q2ExtFillHoleConstrd ( GHoleDomaind *domain,
                         int spdimen, CONST_ double *hole_cp,
                         int nconstr, CONST_ double *constr,
                         double *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const double *cp,
                                            void *usrptr ) );
boolean g1h_Q2ExtFillHoleAltConstrd ( GHoleDomaind *domain, 
                         int spdimen, CONST_ double *hole_cp,
                         int naconstr, CONST_ double *constr,
                         double *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const double *cp,
                                            void *usrptr ) );

boolean g1h_Q2NLExtFillHoleConstrd ( GHoleDomaind *domain, CONST_ point3d *hole_cp,
                                  int nconstr, CONST_ vector3d *constr,
                                  double *acoeff, void *usrptr,
                                  void (*outpatch) ( int n, int m,
                                                     const point3d *cp,
                                                     void *usrptr ) );

boolean g1h_Q2NLExtFillHoleAltConstrd ( GHoleDomaind *domain, CONST_ point3d *hole_cp,
                                     int naconstr, CONST_ double *constr,
                                     double *acoeff, void *usrptr,
                                     void (*outpatch) ( int n, int m,
                                                        const point3d *cp,
                                                        void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* spline basis procedures */
boolean g1h_ComputeSplBasisd ( GHoleDomaind *domain, int nk, int m1, int m2 );
boolean g1h_ComputeSplFormMatrixd ( GHoleDomaind *domain );
boolean g1h_DecomposeSplMatrixd ( GHoleDomaind *domain );
boolean g1h_SplFillHoled ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) );

int g1h_SplV0SpaceDimd ( GHoleDomaind *domain );

boolean g1h_Q2SplComputeFormMatrixd ( GHoleDomaind *domain );
boolean g1h_Q2SplDecomposeMatrixd ( GHoleDomaind *domain );
boolean g1h_Q2SplFillHoled ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) );

boolean g1h_SetSplConstraintMatrixd ( GHoleDomaind *domain,
                                      int nconstr, const double *cmat );
boolean g1h_SplFillHoleConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int nconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) );
boolean g1h_Q2SplFillHoleConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int nconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) );

boolean g1h_SetSplAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                         int naconstr, const double *acmat );
boolean g1h_SplFillHoleAltConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int naconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) );
boolean g1h_Q2SplFillHoleAltConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int naconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* G1 quasi G2 nonlinear construction procedures */
boolean g1h_Q2NLFillHoled ( GHoleDomaind *domain,
                    const point3d *hole_cp,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) );

boolean g1h_Q2NLExtFillHoled ( GHoleDomaind *domain,
                    const point3d *hole_cp,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) );

boolean g1h_Q2NLSplFillHoled ( GHoleDomaind *domain,
                    const point3d *hole_cp,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int lknu, const double *knu,
                                       int m, int lknv, const double *knv,
                                       const point3d *cp, void *usrptr ) );
boolean g1h_Q2NLSplFillHoleConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                     int nconstr, const vector3d *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const double *knu,
                                        int m, int lknv, const double *knv,
                                        const point3d *cp, void *usrptr ) );
boolean g1h_Q2NLSplFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                     int nconstr, const double *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const double *knu,
                                        int m, int lknv, const double *knv,
                                        const point3d *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* drawing procedures */
void g1h_DrawDomAuxPatchesd ( GHoleDomaind *domain,
               void (*drawpatch) ( int n, int m, const point2d *cp ) );
void g1h_DrawBasAuxPatchesd ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int m, const double *cp ) );
boolean g1h_DrawJFunctiond ( GHoleDomaind *domain, int k, int l,
                             void (*drawpoly) ( int deg, const double *f ) );
void g1h_DrawDiPatchesd ( GHoleDomaind *domain,
                      void (*drawpatch) ( int n, int m, const point2d *cp ) );
void g1h_ExtractPartitiond ( GHoleDomaind *domain,
                             int *hole_k, int *hole_m,
                             double *partition,
                             double *part_delta,
                             double *spart_alpha,
                             double *spart_malpha,
                             double *spart_salpha,
                             double *spart_knot,
                             double *alpha0,
                             boolean *spart_sgn,
                             boolean *spart_both );
void g1h_ExtractCentralPointd ( GHoleDomaind *domain,
                                point2d *centp, vector2d *centder );
void g1h_DrawBasAFunctiond ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3d *cp ) );
void g1h_DrawBasBFunctiond ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3d *cp ) );
boolean g1h_DrawBasCNetd ( GHoleDomaind *domain, int fn,
               void (*drawnet) ( int n, int m, const point3d *cp ) );
void g1h_DrawBFAomcd ( GHoleDomaind *domain, int fn, 
                       void (*drawpoly)(int degree, const double *coeff) );
void g1h_DrawBFBomcd ( GHoleDomaind *domain, int fn,
                       void (*drawpoly)(int degree, const double *coeff) );
void g1h_DrawFinalSurfBCd ( GHoleDomaind *domain,   
                            int spdimen, const double *hole_cp,
                            const double *acoeff, 
                            void (*drawcurve)(int degree, int spdimen,
                                              const double *cp) );
void g1h_ExtDrawFinalSurfBCd ( GHoleDomaind *domain,
                               int spdimen, const double *hole_cp,
                               const double *acoeff, 
                               void (*drawcurve)(int degree, int spdimen,
                                                 const double *cp) );
void g1h_DrawMatricesd ( GHoleDomaind *domain,
                         void (*drawmatrix)(int nfa, int nfb,
                                            double *amat, double *bmat) );
void g1h_DrawExtMatricesd ( GHoleDomaind *domain,
                            void (*drawmatrix)(int k, int r, int s, double *Aii,
                                               double *Bi) );
int g1h_DrawBFcpnd ( int hole_k, unsigned char *bfcpn );

void g1h_Q2DrawMatricesd ( GHoleDomaind *domain,
                           void (*drawmatrix)(int nfa, int nfb,
                                              double *amat, double *bmat) );
void g1h_Q2DrawExtMatricesd ( GHoleDomaind *domain,
                              void (*drawmatrix)(int k, int r, int s,
                                                 double *Aii, double *Bi) );
boolean g1h_GetABasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                    double *bfpc );
boolean g1h_GetBBasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                    double *bfpc );

/* ///////////////////////////////////////////////////////////////////////// */
/* spline basis drawing procedures */
void g1h_DrawSplBasFuncNumd ( GHoleDomaind *domain,
                        int *nfunc_a, int *nfunc_b, int *nfunc_c, int *nfunc_d );

void g1h_DrawSplBasAuxPatchesd ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int lknu, const double *knu,
                                   int m, int lknv, const double *knv,
                                   const point3d *cp ) );

void g1h_DrawSplBasFunctiond ( GHoleDomaind *domain, int fn,
             void (*drawpatch) ( int n, int lknu, const double *knu,
                                 int m, int lknv, const double *knv,
                                 const point3d *cp ) );

void g1h_DrawSplBFAomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) );
void g1h_DrawSplBFBomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) );
void g1h_DrawSplBFDomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) );
void g1h_DrawSplFinalSurfBCd ( GHoleDomaind *domain,
                               int spdimen, const double *hole_cp,
                               const double *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                             int lastknot, const double *knots,
                                             const double *cp) );

void g1h_DrawSplMatricesd ( GHoleDomaind *domain,
                            void (*drawmatrix)( int k, int r, int s, int t,
                                                double *A, double *B ) );

void g1h_Q2DrawSplMatricesd ( GHoleDomaind *domain,
                              void (*drawmatrix)(int k, int r, int s,
                                                 double *Aii, double *Bi) );

boolean g1h_GetSplABasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                       int *lkn, double *kn, double *bfpc );
boolean g1h_GetSplBBasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                       int *lkn, double *kn, double *bfpc );
boolean g1h_GetSplDBasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                       int *lkn, double *kn, double *bfpc );

/* ///////////////////////////////////////////////////////////////////////// */
int g1h_SymPatchMatrixSize ( int hole_k );
boolean g1h_GetSymPatchMatrixd ( GHoleDomaind *domain, double *patchmatrix );
boolean g1h_GetExtSymPatchMatrixd ( GHoleDomaind *domain, double *patchmatrix );
boolean g1h_Q2GetSymPatchMatrixd ( GHoleDomaind *domain, double *patchmatrix );
boolean g1h_Q2GetExtSymPatchMatrixd ( GHoleDomaind *domain, double *patchmatrix );
boolean g1h_MatrixFillSymHoled ( int hole_k, const double *patchmatrix,
                                 int spdimen, const double *hole_cp, void *usrptr,
                                 void (*outpatch) ( int n, int m, const double *cp,
                                                    void *usrptr ) );
/* ///////////////////////////////////////////////////////////////////////// */
int g1h_GetErrorCoded ( GHoleDomaind *domain, char **ErrorString );

#ifdef __cplusplus
}
#endif

#endif

