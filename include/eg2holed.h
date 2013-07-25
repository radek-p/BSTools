
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2010                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* this header file is for application use */

#ifndef EG2HOLED_H
#define EG2HOLED_H

#ifndef EGHOLED_H
#include "egholed.h"
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
#define G2H_CENRTAL_DERIVATIVES1_ALT   1
#define G2H_CENTRAL_DERIVATIVES1_GIVEN 2

#define G2HQUERY_DOMAIN_CURVES         3
#define G2H_DOMAIN_CURVES_DEG4         1

#define G2HQUERY_BASIS                 4
#define G2H_USE_RESTRICTED_BASIS       1

#define G2HQUERY_QUADRATURE            5
#define G2H_QUADRATURE_GAUSS_LEGENDRE  1

/* ///////////////////////////////////////////////////////////////////////// */
/* core procedures */

void g2h_SetOptionProcd ( GHoleDomaind *domain,
    int (*OptionProc)( GHoleDomaind *domain, int query, int qn,
                       int *ndata, int **idata, double **fdata ) );

boolean g2h_ComputeBasisd ( GHoleDomaind *domain );

boolean g2h_ComputeFormMatrixd ( GHoleDomaind *domain );
boolean g2h_DecomposeMatrixd ( GHoleDomaind *domain );
boolean g2h_FillHoled ( GHoleDomaind *domain,
                        int spdimen, CONST_ double *hole_cp,
                        double *acoeff, void *usrptr,
                        void (*outpatch) ( int n, int m, const double *cp,
                                           void *usrptr ) );

boolean g2h_ComputeExtFormMatrixd ( GHoleDomaind *domain );
boolean g2h_DecomposeExtMatrixd ( GHoleDomaind *domain );
boolean g2h_ExtFillHoled ( GHoleDomaind *domain,
                           int spdimen, CONST_ double *hole_cp,
                           double *acoeff, void *usrptr,
                           void (*outpatch) ( int n, int m, const double *cp,
                                              void *usrptr ) );

int g2h_V0SpaceDimd ( GHoleDomaind *domain );
int g2h_ExtV0SpaceDimd ( GHoleDomaind *domain );
boolean g2h_GetBPDerivativesd ( GHoleDomaind *domain, int cno, double *val );
boolean g2h_GetBFuncPatchd ( GHoleDomaind *domain, int fn, int pn, double *bp );

boolean g2h_SetConstraintMatrixd ( GHoleDomaind *domain,
                                   int nconstr, const double *cmat );
boolean g2h_FillHoleConstrd ( GHoleDomaind *domain,
                              int spdimen, CONST_ double *hole_cp,
                              int nconstr, CONST_ double *constr,
                              double *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const double *cp,
                                                 void *usrptr ) );

boolean g2h_SetAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                      int nconstr, const double *cmat );
boolean g2h_FillHoleAltConstrd ( GHoleDomaind *domain,
                              int spdimen, CONST_ double *hole_cp,
                              int naconstr, CONST_ double *constr,
                              double *acoeff, void *usrptr,
                              void (*outpatch) ( int n, int m, const double *cp,
                                                 void *usrptr ) );

boolean g2h_SetExtConstraintMatrixd ( GHoleDomaind *domain,
                                      int nconstr, const double *cmat );
boolean g2h_ExtFillHoleConstrd ( GHoleDomaind *domain,
                         int spdimen, CONST_ double *hole_cp,
                         int nconstr, CONST_ double *constr,
                         double *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const double *cp,
                                            void *usrptr ) );

boolean g2h_SetExtAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                      int naconstr, const double *acmat );
boolean g2h_ExtFillHoleAltConstrd ( GHoleDomaind *domain,
                         int spdimen, CONST_ double *hole_cp,
                         int naconstr, CONST_ double *constr,
                         double *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const double *cp,
                                            void *usrptr ) );

double g2h_FunctionalValued ( GHoleDomaind *domain, int spdimen,
                             CONST_ double *hole_cp, CONST_ double *acoeff );
double g2h_ExtFunctionalValued ( GHoleDomaind *domain, int spdimen,
                                CONST_ double *hole_cp, CONST_ double *acoeff );
boolean g2h_NLFunctionalValued ( GHoleDomaind *domain,   
                                 const point3d *hole_cp, const vector3d *acoeff,
                                 double *funcval );
boolean g2h_NLExtFunctionalValued ( GHoleDomaind *domain,
                                    const point3d *hole_cp, const vector3d *acoeff,
                                    double *funcval );

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_ComputeNLNormald ( GHoleDomaind *domain, 
                               const point3d *hole_cp,
                               vector3d *anv );

boolean g2h_NLFillHoled ( GHoleDomaind *domain, const point3d *hole_cp,
                          double *acoeff, void *usrptr,
                          void (*outpatch) ( int n, int m, const point3d *cp,
                                             void *usrptr ) );
boolean g2h_NLFillHoleConstrd ( GHoleDomaind *domain, const point3d *hole_cp,   
                    int nconstr, const vector3d *constr,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) );
boolean g2h_NLFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                    int nconstr, const double *constr,
                    double *acoeff, void *usrptr,
                    void (*outpatch) ( int n, int m, const point3d *cp,
                                       void *usrptr ) );

boolean g2h_NLExtFillHoled ( GHoleDomaind *domain, const point3d *hole_cp,
                             double *acoeff, void *usrptr,
                             void (*outpatch) ( int n, int m, const point3d *cp,
                                                void *usrptr ) );
boolean g2h_NLExtFillHoleConstrd ( GHoleDomaind *domain,
                     const point3d *hole_cp,
                     int nconstr, const vector3d *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int m, const point3d *cp,
                                        void *usrptr ) );
boolean g2h_NLExtFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                         int naconstr, const double *constr,
                         double *acoeff, void *usrptr,
                         void (*outpatch) ( int n, int m, const point3d *cp,
                                            void *usrptr ) ); 

boolean g2h_NLSplFillHoled ( GHoleDomaind *domain, const point3d *hole_cp,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const double *knu,
                                        int m, int lknv, const double *knv,
                                        const point3d *cp, void *usrptr ) );
boolean g2h_NLSplFillHoleConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                     int nconstr, const vector3d *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const double *knu,
                                        int m, int lknv, const double *knv,
                                        const point3d *cp, void *usrptr ) );
boolean g2h_NLSplFillHoleAltConstrd ( GHoleDomaind *domain, const point3d *hole_cp,
                     int nconstr, const double *constr,
                     double *acoeff, void *usrptr,
                     void (*outpatch) ( int n, int lknu, const double *knu,
                                        int m, int lknv, const double *knv,
                                        const point3d *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2h_GetFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                   CONST_ double *hole_cp, double *acoeff,
                                   void (*outcurve) ( int n, const double *cp ) );
boolean g2h_GetExtFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                      CONST_ double *hole_cp, double *acoeff,
                                      void (*outcurve) ( int n, const double *cp ) );
boolean g2h_GetSplFinalPatchCurvesd ( GHoleDomaind *domain, int spdimen,
                                      CONST_ double *hole_cp, double *acoeff,
                                      void (*outcurve) ( int n, int lkn,
                                               const double *kn, const double *cp ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* spline basis procedures */
boolean g2h_ComputeSplBasisd ( GHoleDomaind *domain, int nk, int m1, int m2 );
boolean g2h_ComputeSplFormMatrixd ( GHoleDomaind *domain );
boolean g2h_DecomposeSplMatrixd ( GHoleDomaind *domain );
boolean g2h_SplFillHoled ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,            
                                  const double *cp, void *usrptr ) );

int g2h_SplV0SpaceDimd ( GHoleDomaind *domain );
boolean g2h_SetSplConstraintMatrixd ( GHoleDomaind *domain,
                                      int nconstr, const double *cmat );
boolean g2h_SplFillHoleConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int nconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) );

boolean g2h_SetSplAltConstraintMatrixd ( GHoleDomaind *domain, int spdimen,
                                         int naconstr, const double *acmat );
boolean g2h_SplFillHoleAltConstrd ( GHoleDomaind *domain,
               int spdimen, CONST_ double *hole_cp,
               int naconstr, CONST_ double *constr,
               double *acoeff, void *usrptr,
               void (*outpatch) ( int n, int lknu, const double *knu,
                                  int m, int lknv, const double *knv,
                                  const double *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
/* drawing procedures */
void g2h_DrawDomAuxPatchesd ( GHoleDomaind *domain,
               void (*drawpatch) ( int n, int m, const point2d *cp ) );
void g2h_DrawBasAuxPatchesd ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int m, const double *cp ) );
void g2h_DrawJFunctiond ( GHoleDomaind *domain, int k, int l,
                          void (*drawpoly) ( int deg, const double *f ) );
void g2h_DrawDiPatchesd ( GHoleDomaind *domain,
                      void (*drawpatch) ( int n, int m, const point2d *cp ) );
void g2h_ExtractPartitiond ( GHoleDomaind *domain,
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
void g2h_ExtractCentralPointd ( GHoleDomaind *domain,
                                point2d *centp, vector2d *centder );
void g2h_DrawBasAFunctiond ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3d *cp ) );
void g2h_DrawBasBFunctiond ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int m, const point3d *cp ) );
void g2h_DrawBasCNetd ( GHoleDomaind *domain, int fn,
               void (*drawnet) ( int n, int m, const point3d *cp ) );
void g2h_DrawBFAomcd ( GHoleDomaind *domain, int fn, 
                       void (*drawpoly)(int degree, const double *coeff) );
void g2h_DrawBFBomcd ( GHoleDomaind *domain, int fn,
                       void (*drawpoly)(int degree, const double *coeff) );
void g2h_DrawFinalSurfBCd ( GHoleDomaind *domain,   
                            int spdimen, const double *hole_cp,
                            const double *acoeff, 
                            void (*drawcurve)(int degree, int spdimen,
                                              const double *cp) );
void g2h_ExtDrawFinalSurfBCd ( GHoleDomaind *domain,
                               int spdimen, const double *hole_cp,
                               const double *acoeff, 
                               void (*drawcurve)(int degree, int spdimen,
                                                 const double *cp) );
void g2h_DrawMatricesd ( GHoleDomaind *domain,
                         void (*drawmatrix)(int nfa, int nfb,
                                            double *amat, double *bmat) );
void g2h_DrawExtMatricesd ( GHoleDomaind *domain,
                            void (*drawmatrix)(int k, int r, int s,
                                               double *Aii, double *Bi) );
int g2h_DrawBFcpnd ( int hole_k, unsigned char *bfcpn );

boolean g2h_GetABasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                    double *bfpc );
boolean g2h_GetBBasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                    double *bfpc );

/* ///////////////////////////////////////////////////////////////////////// */
/* spline basis drawing procedures */
void g2h_DrawSplBasFuncNumd ( GHoleDomaind *domain,
                        int *nfunc_a, int *nfunc_b, int *nfunc_c, int *nfunc_d );

void g2h_DrawSplBasAuxPatchesd ( GHoleDomaind *domain, int fn,
               void (*drawpatch) ( int n, int lknu, const double *knu,
                                   int m, int lknv, const double *knv,
                                   const point3d *cp ) );

void g2h_DrawSplBasFunctiond ( GHoleDomaind *domain, int fn,
             void (*drawpatch) ( int n, int lknu, const double *knu,
                                 int m, int lknv, const double *knv,
                                 const point3d *cp ) );

void g2h_DrawSplBFAomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) );
void g2h_DrawSplBFBomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) );
void g2h_DrawSplBFDomcd ( GHoleDomaind *domain, int fn,
                          void (*drawpoly)(int degree, int lastknot,
                                           const double *knots,
                                           const double *coeff) );
void g2h_DrawSplFinalSurfBCd ( GHoleDomaind *domain,
                               int spdimen, const double *hole_cp,
                               const double *acoeff,
                               void (*drawcurve)(int degree, int spdimen,
                                             int lastknot, const double *knots,
                                             const double *cp) );

void g2h_DrawSplMatricesd ( GHoleDomaind *domain,
                            void (*drawmatrix)( int k, int r, int s, int t,
                                                double *A, double *B ) );

boolean g2h_GetSplABasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                       int *lkn, double *kn, double *bfpc );
boolean g2h_GetSplBBasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                       int *lkn, double *kn, double *bfpc );
boolean g2h_GetSplDBasisFPatchCurved ( GHoleDomaind *domain, int fn, int i,
                                       int *lkn, double *kn, double *bfpc );

/* ///////////////////////////////////////////////////////////////////////// */
int g2h_SymPatchMatrixSize ( int hole_k );
boolean g2h_GetSymPatchMatrixd ( GHoleDomaind *domain, double *patchmatrix );
boolean g2h_GetExtSymPatchMatrixd ( GHoleDomaind *domain, double *patchmatrix );
boolean g2h_MatrixFillSymHoled ( int hole_k, const double *patchmatrix,
                                 int spdimen, const double *hole_cp, void *usrptr,
                                 void (*outpatch) ( int n, int m, const double *cp,
                                                    void *usrptr ) );
/* ///////////////////////////////////////////////////////////////////////// */
int g2h_GetErrorCoded ( GHoleDomaind *domain, char **ErrorString );

#ifdef __cplusplus
}
#endif

#endif

