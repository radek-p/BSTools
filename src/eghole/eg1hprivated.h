
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* this header file is private for g1hole library functions; */
/* it is NOT intended to be #included in application source files */

#ifndef EG1HPRIVATED_H
#define EG1HPRIVATED_H

#include "eghprivated.h"
#include "eg1hprivate.h"

/* degrees of polynomial junction functions */
#define G1_BF01DEG 2
#define G1_CG01DEG 1
#define G1_BF11DEG 3
#define G1_CG11DEG 2

/* degrees of patch cross derivatives */
#define G1_CROSS00DEG 4
#define G1_CROSS01DEG 5
#define G1_CROSS10DEG 3
#define G1_CROSS11DEG 5

/* degree of the final patches - the greatest of the above four */
/* #define G1H_FINALDEG  5 - this is #defined in g1holef.h */

#define G1_CROSSDEGSUM \
        (G1_CROSS00DEG+G1_CROSS01DEG+G1_CROSS10DEG+G1_CROSS11DEG)

/* the following number is relevant for spline bases */
#define G1_AUXDEG0 (G1_CROSS01DEG-G1H_OMCDEG+1)

/* angle tolerances for partition analysis */
#define MIN_DELTA (PI/180.0)       /*   1 degree  */
#define MAX_DELTA (179.0*PI/180.0) /* 179 degrees */
#define ANGLETOL (0.001*PI/180)    /* 0.001 degree */

/* the #definitions below are for the construction with */
/* the extended basis */
#define G1_DBDIM 4    /* diagonal block dimension, 2*2 */
#define G1_DIAGBLSIZE (G1_DBDIM*(G1_DBDIM+1)/2)


/* number of quadrature knots */
#define G1_NQUAD 16
#define G1_NQUADSQ (G1_NQUAD*G1_NQUAD)
#define G1_QUAD_FACTOR 10  /* for integration of spline basis functions */


/* ////////////////////////////////////////////////////////////////////////// */
/* private data of the construction procedures */
typedef struct G1HolePrivateRecd {
    int      hole_m;
    int      nfunc_a,   /* number of "inner" basis functions */
             nfunc_b;   /* number of "outer" basis functions */
    char     opt_cp, opt_der1, opt_basis, opt_quad;  /* recorded options */
    vector2d *omcbc;    /* boundary conditions for auxiliary patches */
    vector2d *omc;      /* auxiliary patches representation */
    double   *jfunc;    /* coefficients of junction functions */
    vector2d *dicross;  /* curves describing cross derivatives */
                        /* of the domain patches */
    double   *AuxiMat;  /* matrices of transformations of derivatives */
                        /* for composing functions with domain patches */
    double   *basis_a,  /* Coons representations of "inner" basis functions */
             *basis_b;  /* Coons represetnationf of "outer" basis functions */
    unsigned char spdimen[2];   /* numbers of "inner" basis functions: */
                                /* "half polynomial" and "splines" of degree 2 */
    double   *partition;        /* partition of the full angle */
    GHoleSgnPartd *spartition;  /* sorted partition, having hole_k-hole_m */
                                /* elements */
    double   spart_alpha0;      /* sorted partition bisectrix */
    double   *Amat, *Bmat;      /* matrices of scalar products */
                                /* of basis functions */
    double   *Lmat;             /* decomposition factor L:  A=LL^T */
    double   *EAmat, *EBmat;    /* matrices for the extended basis */
    double   *ELmat;            /* decomposition factor for the */
                                /* extended matrix */

    int      nconstr;           /* number of constraint equations */
    double   *Cmat, *RCmat;     /* matrices of constraint equations */
    int      naconstr, acdim;
    double   *ACmat, *ARCmat;
    int      extnconstr;        /* number of constraints and */
    double   *ECmat, *ERCmat;   /* matrices of constraint equations for */
                                /* the extended matrix */
    int      extnaconstr, extacdim;
    double   *AECmat, *AERCmat;

    double   *BBmat;            /* matrix of scalar products, for */
                                /* evaluation of the functional */

                              /* private data for the quasi G2 constructions */
    double    C1b;                 /* constant in the definition of the form */
    double    *Q2AMat, *Q2BMat;    /* form matrices */
    double    *Q2LMat;
    double    *Q2RCMat, *Q2ARCMat; /* constraints */

    double    C1e;                 /* constant in the definition of the form */
    double    *Q2EAMat, *Q2EBMat;  /* form matrices for the extended basis */
    double    *Q2ELMat;
    double    *Q2ERCMat, *Q2AERCMat; /* constraints */
    int (*GetOption)( GHoleDomaind *domain, int query, int qn,
                      int *ndata, int **idata, double **fdata );
  } G1HolePrivateRecd;


/* private data for the spline basis construction */
typedef struct G1HoleSPrivateRecd {
    int    nk, m1, m2;   /* parameters, which determine the space dimension */
    int    nsfunc_c, nsfunc_d;
                       /* numbers of basis functions of various kinds */
    int    csize, dsize; /* sizes of the blocks in the form matrix */
    int    lastbezknot;
    double bezknots[2*(G1H_FINALDEG+1)];  /* knots for Bezier polynomials */
    int    lastomcknot;
    double *omcknots;  /* knots of basis function auxiliary patches */
    int    lastpvknot;
    double *pvknots;
                      /* knots of the basis functions of the D block */
    int    lastcknot;
    double *cknots;    /* knots of the basis functions of the C block */
    double *basis_d;   /* Coons representations of basis functions */   

    int    lastfpknot;
    double *fpknots;   /* knots of the final patches */
    double *SAMat, *SBMat;  /* form matrices for the spline basis */
    double *SLMat;          /* A=LL^T decomposition factor */

    int    splnconstr; /* number of constraints */
    double *SCmat, *SRCmat; /* matrices of constraint equations */
    int    splnaconstr, splacdim;
    double *ASCmat, *ASRCmat;

/* private data for the quasi G2 constructions */
    double C1s;                 /* constant in the definition of the form */
    double *Q2SAMat, *Q2SBMat;  /* form matrices for the spline basis */
    double *Q2SLMat;

    double *Q2SRCmat; /* matrices of constraint equations */
    double *Q2SARCmat;
  } G1HoleSPrivateRecd;


/* private data for the nonlinear constructions */
typedef struct G1HNLPrivated {
    int      auxc;    /* patch counter */
    vector3d *nldi;   /* NL patches */
    vector3d *acoeff;
    vector3d *rhole_cp;
    vector3d nlnv, reflv;
    vector2d *diu, *div, *diuu, *diuv, *divv;
    double   ddiam;   /* NL domain diameter */
    double   *jac;
    double   *psiu, *psiv, *psiuu, *psiuv, *psivv;
                   /* the following data are used only for Q2 */
    vector2d *diuuu, *diuuv, *diuvv, *divvv;
    double   *psiuuu, *psiuuv, *psiuvv, *psivvv;
    vector2d *ctang;
    double   *cpsiu, *cpsiv, *cpsiuu, *cpsiuv, *cpsivv;
  } G1HNLPrivated;  


typedef struct G1Q2HNLFuncd {
    double   pu, pv, jpuu, jpuv, jpvv;
    vector2d tang;
    double   psiu, psiv, jpsiuu, jpsiuv, jpsivv;
    double   psju, psjv, jpsjuu, jpsjuv, jpsjvv;  
    double   A, B, e0, e1, e2, e9, tGt, b1, b2, b3, b4, b5;
  } G1Q2HNLFuncd;


/* private data for spline nonlinear constructions */
typedef struct G1HNLSPrivated {
    G1HNLPrivated nlpr; /* must be the first field */
    int     nkn;          /* number of quadrature knots */
    double  *tkn;         /* quadrature knots */
    int     ftabsize;     /* function values table size */
    int     psize;        /* number of control points of each final patch  */
    int     *fkn, *lkn;   /* ranges of knots for the functions of block C  */
    double  *cb, *cbt, *cbtt, *cbttt; /* B-spline functions and their derivatives */
    int     *cfuncvi;     /* indexes of the first samples of C block functions */
    int     *dfuncvi;     /* indexes of the first samples of D block functions */
    int     jtabsize;     /* jump values table size */
    int     jcfs, jdfs;   /* indexes of the first samples of jump */
                          /* for block C and block D functions */
    short   njcurves;
    boolean jumpC, jumpD;
  } G1HNLSPrivated;


/* ////////////////////////////////////////////////////////////////////////// */
extern G1HNLPrivated *_g1h_nlprivd;

/* ////////////////////////////////////////////////////////////////////////// */
void _g1h_GetDiPatchCurvesd ( GHoleDomaind *domain, int i,
                point2d **c00, vector2d **c01, point2d **c10, vector2d **c11,
                point2d **d00, vector2d **d01, point2d **d10, vector2d **d11 );
boolean _g1h_GetABasisAuxpd ( GHoleDomaind *domain, int fn,
                              double *br0, double *br0cr1 );
boolean _g1h_GetBBasisAuxpd ( GHoleDomaind *domain, int fn,
                              double *bezfc, double *fcomc, double *fcomcd );
void _g1h_GetBFAPatchCurvesd ( GHoleDomaind *domain, int fn, int i,
                    double **c00, double **c01,
                    double **d00, double **d01 );
void _g1h_GetBFBPatchCurvesd ( GHoleDomaind *domain, int fn, int i,
                    double **c00, double **c01, double **c10, double **c11,
                    double **d00, double **d01, double **d10, double **d11 );

void g1h_TabCubicHFuncDer2d ( int nkn, const double *kn,
                              double *hfunc, double *dhfunc, double *ddhfunc );
boolean g1h_TabBicubicCoonsPatchDer2d (
      int spdimen, int nkn, const double *kn, const double *hfunc,
      const double *dhfunc, const double *ddhfunc,
      int degc00, const double *c00,
      int degc01, const double *c01,
      int degc10, const double *c10,
      int degc11, const double *c11,
      int degd00, const double *d00,
      int degd01, const double *d01,
      int degd10, const double *d10,
      int degd11, const double *d11,
      double *p, double *pu, double *pv, double *puu, double *puv, double *pvv );
void _g1h_DiJacobian2d ( const vector2d *du, const vector2d *dv,
                         const vector2d *duu, const vector2d *duv,
                         const vector2d *dvv,
                         double *jac, double *trd );
boolean _g1h_TabDiPatchJac2d ( int nkn, const double *kn, const double *hfunc,
                               const double *dhfunc, const double *ddhfunc,
                               const vector2d *c00, const vector2d *c01,
                               const vector2d *c10, const vector2d *c11,
                               const vector2d *d00, const vector2d *d01,
                               const vector2d *d10, const vector2d *d11,
                               double *jac, double *trd );
boolean _g1h_TabLaplaciand ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *fc00, const double *fc01,
        const double *fc10, const double *fc11,
        const double *fd00, const double *fd01,
        const double *fd10, const double *fd11,
        const double *trd,
        double *lap );
boolean _g1h_TabLaplacian0d ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *fc00, const double *fc01,
        const double *fd00, const double *fd01,
        const double *trd,
        double *lapgrad );
double _g1h_Integrald ( int hole_k, int nknsq, double *jac,
                       unsigned short supp1, double *func1,
                       unsigned short supp2, double *func2 );

boolean _g1h_VerifyJunctionFunctionsd ( GHoleDomaind *domain );
boolean _g1h_VerifyDomPatchesd ( GHoleDomaind *domain );
                                                    
/* ////////////////////////////////////////////////////////////////////////// */
void g1h_ReflectVectorsd ( int n, const vector3d *v, vector3d *w );
void g1h_nonlinoutpatchd ( int n, int m, const double *cp, void *usrptr );
boolean _g1h_StopItd ( int itn, double gn0, double gn,
                       double cn, double dcn, double scf );
boolean g1h_GetHoleSurrndPatchd ( GHoleDomaind *domain,
                                  const point3d *hole_cp,
                                  int i, int j, point3d *bcp );
boolean _g1h_ComputeNLNormald ( GHoleDomaind *domain,
                                G1HNLPrivated *nlprivate,
                                const point3d *hole_cp );
boolean _g1h_TabNLDer0d ( int nkn, const double *tkn,
             const double *hfunc, const double *dhfunc, const double *ddhfunc,
             vector2d *diu, vector2d *div,
             vector2d *diuu, vector2d *diuv, vector2d *divv,
             double *fc00, double *fc01, double *fd00, double *fd01,
             double *psiu, double *psiv,
             double *psiuu, double *psiuv, double *psivv );
boolean _g1h_TabNLDerd ( int nkn, double *tkn,
             const double *hfunc, const double *dhfunc, const double *ddhfunc,
             vector2d *diu, vector2d *div,
             vector2d *diuu, vector2d *diuv, vector2d *divv,
             double *fc00, double *fc01, double *fc10, double *fc11,
             double *fd00, double *fd01, double *fd10, double *fd11,
             double *psiu, double *psiv,
             double *psiuu, double *psiuv, double *psivv );
boolean _g1h_TabNLBasisFunctionsd ( GHoleDomaind *domain, int nkn,
                                    G1HNLPrivated *nlpr );
boolean _g1h_ReflectSplAltConstrMatrixd ( GHoleDomaind *domain,
                                          vector3d *reflv, double *RACmat );
void _g1h_IntFunc1d ( double pu, double pv,
                      double puu, double puv, double pvv, double jac,
                      double *c, double *cs, double *a, double *b,   
                      double *funct );
void _g1h_IntFunc2d ( double pu, double pv,
                      double puu, double puv, double pvv,
                      double psiu, double psiv,
                      double psiuu, double psiuv, double psivv,
                      double jac, double c, double A, double B, 
                      double *ai, double *bi, double *grad );
void _g1h_IntFunc3d ( double pu, double pv,       
                      double puu, double puv, double pvv,       
                      double psiu, double psiv,       
                      double psiuu, double psiuv, double psivv,       
                      double psju, double psjv,       
                      double psjuu, double psjuv, double psjvv,       
                      double jac, double c,       
                      double A, double B, double Ai, double Bi,       
                      double Aj, double Bj,       
                      double *hessian );

/* ////////////////////////////////////////////////////////////////////////// */
boolean _g1hq2_TabNLDer0d ( int nkn, const double *tkn,
             const double *hfunc, const double *dhfunc, const double *ddhfunc,
             const double *dddhfunc,
             vector2d *diu, vector2d *div,
             vector2d *diuu, vector2d *diuv, vector2d *divv,
             vector2d *diuuu, vector2d *diuuv, vector2d *diuvv, vector2d *divvv,
             double *fc00, double *fc01, double *fd00, double *fd01,
             double *psiu, double *psiv,
             double *psiuu, double *psiuv, double *psivv,
             double *psiuuu, double *psiuuv, double *psiuvv, double *psivvv );

boolean _g1hq2_TabNLDerd ( int nkn, double *tkn,
             const double *hfunc, const double *dhfunc, const double *ddhfunc,
             const double *dddhfunc,
             vector2d *diu, vector2d *div,
             vector2d *diuu, vector2d *diuv, vector2d *divv,
             vector2d *diuuu, vector2d *diuuv, vector2d *diuvv, vector2d *divvv,
             double *fc00, double *fc01, double *fc10, double *fc11,
             double *fd00, double *fd01, double *fd10, double *fd11,
             double *psiu, double *psiv,
             double *psiuu, double *psiuv, double *psivv,
             double *psiuuu, double *psiuuv, double *psiuvv, double *psivvv );

boolean _g1hq2_TabNLBasisFunctionsOmegad ( GHoleDomaind *domain, int nkn,
                                           G1HNLPrivated *nlpr, double *bc00 );
boolean _g1hq2_TabNLBasisFunctionsGammad ( GHoleDomaind *domain, int nkn,
                                           G1HNLPrivated *nlpr, double *ctrd,
                                           double *bc00 );
void _g1hq2_IntFunc1ad ( G1Q2HNLFuncd *f, double *funct );
void _g1hq2_IntFunc1bd ( G1Q2HNLFuncd *f, double *funct );
void _g1hq2_IntFunc1cd ( G1Q2HNLFuncd *f, double *funct );
void _g1hq2_IntFunc2bd ( G1Q2HNLFuncd *f, double *grad );
void _g1hq2_IntFunc2cd ( G1Q2HNLFuncd *f, double *Ai, double *Bi, double *grad );
void _g1hq2_IntFunc3cd ( G1Q2HNLFuncd *f, double Ai, double Bi, double Aj, double Bj,
                         double *hessian );

boolean _g1h_ReflectAltConstrMatrixd ( GHoleDomaind *domain,
                                       vector3d *reflv, double *RACmat );
boolean _g1h_ReflectExtAltConstrMatrixd ( GHoleDomaind *domain,
                                          vector3d *reflv, double *RACmat );

/* ////////////////////////////////////////////////////////////////////////// */
/* the following horrible macros are used a number of times in the code. */
/* They assume that some variables are declared and have proper initial */
/* values */
#define G1GetPolyAddr(jfunc,b01,c01,f01,g01,b11,c11,f11,g11)\
  b01 = jfunc; \
  c01 = &b01[hole_k*(G1_BF01DEG+1)];  f01 = &c01[hole_k*(G1_CG01DEG+1)]; \
  g01 = &f01[hole_k*(G1_BF01DEG+1)];  b11 = &g01[hole_k*(G1_CG01DEG+1)]; \
  c11 = &b11[hole_k*(G1_BF11DEG+1)];  f11 = &c11[hole_k*(G1_CG11DEG+1)]; \
  g11 = &f11[hole_k*(G1_BF11DEG+1)]

#define G1GetPolyAddr0(jfunc,b01,c01,f01,g01)\
  b01 = jfunc; \
  c01 = &b01[hole_k*(G1_BF01DEG+1)];  f01 = &c01[hole_k*(G1_CG01DEG+1)]; \
  g01 = &f01[hole_k*(G1_BF01DEG+1)]

#define G1GetDiCrossAddresses() \
  dir0cr1 = privateG1->dicross; \
  diq0cr1 = &dir0cr1[hole_k*(G1_CROSS01DEG+1)]; \
  dir1cr1 = &diq0cr1[hole_k*(G1_CROSS01DEG+1)]; \
  diq1cr1 = &dir1cr1[hole_k*(G1_CROSS11DEG+1)]

#define G1GetBFuncACrossAddresses() \
  bbr0 = privateG1->basis_a; \
  bbr0cr1 = &bbr0[nfunc_a*hole_k*(G1_CROSS00DEG+1)]; \
  bbq0cr1 = &bbr0cr1[nfunc_a*hole_k*(G1_CROSS01DEG+1)]

#define G1GetBFuncBCrossAddresses() \
  bbr0 = privateG1->basis_b; \
  bbr1 = &bbr0[nfunc_b*hole_k*(G1_CROSS00DEG+1)]; \
  bbq1 = &bbr1[nfunc_b*hole_k*(G1_CROSS10DEG+1)]; \
  bbr0cr1 = &bbq1[nfunc_b*hole_k*(G1_CROSS10DEG+1)]; \
  bbq0cr1 = &bbr0cr1[nfunc_b*hole_k*(G1_CROSS01DEG+1)]; \
  bbr1cr1 = &bbq0cr1[nfunc_b*hole_k*(G1_CROSS01DEG+1)]; \
  bbq1cr1 = &bbr1cr1[nfunc_b*hole_k*(G1_CROSS11DEG+1)]

#define G1GetFCAddresses() \
  fc01 = &fc00[hole_k*spdimen*(G1_CROSS00DEG+1)]; \
  fc10 = &fc01[hole_k*spdimen*(G1_CROSS01DEG+1)]; \
  fc11 = &fc10[hole_k*spdimen*(G1_CROSS10DEG+1)]; \
  fd00 = &fc11[hole_k*spdimen*(G1_CROSS11DEG+1)]; \
  fd01 = &fd00[hole_k*spdimen*(G1_CROSS00DEG+1)]; \
  fd10 = &fd01[hole_k*spdimen*(G1_CROSS01DEG+1)]; \
  fd11 = &fd10[hole_k*spdimen*(G1_CROSS10DEG+1)]

/* ///////////////////////////////////////////////////////////////////////// */
#define G1GetSFCAddresses() \
  sfc01 = &sfc00[hole_k*spdimen*(lastomcknot-G1_CROSS00DEG)]; \
  sfd00 = &sfc01[hole_k*spdimen*(lastpvknot-G1_CROSS01DEG)]; \
  sfd01 = &sfd00[hole_k*spdimen*(lastomcknot-G1_CROSS00DEG)]

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_GetExtBlockAddressesd ( GHoleDomaind *domain,
                                     double **Aii, double **Aki, double **Akk,
                                     double **Bi, double **Bk, double **Lii );

boolean _g1h_SetRightSided ( GHoleDomaind *domain, const double *Bmat,
                             int spdimen, CONST_ double *hole_cp,
                             double *fc00, double *b );
boolean _g1h_OutputPatchesd ( GHoleDomaind *domain, int spdimen,
                       CONST_ double *x, double *fc00, void *usrptr,
                       void (*outpatch) ( int n, int m, const double *cp,
                                          void *usrptr ) );

boolean _g1h_SetExtRightSided ( GHoleDomaind *domain,
                                const double *Bo, const double *Bk,
                                int spdimen, CONST_ double *hole_cp,
                                double *fc00, double *b );
boolean _g1h_OutputExtPatchesd ( GHoleDomaind *domain, int spdimen,
                           CONST_ double *x, double *fc00, void *usrptr,
                           void (*outpatch) ( int n, int m, const double *cp,
                                              void *usrptr ) );

boolean _g1h_TabTensBezPolyDer2d ( int nkn, const double *tkn,       
                     double *tbez, double *tbezu, double *tbezv,
                     double *tbezuu, double *tbezuv, double *tbezvv );

void _g1h_TensDer2d ( double p, double pu, double puu,
                      double q, double qv, double qvv,
                      double *pq );

/* ///////////////////////////////////////////////////////////////////////// */
void g1h_DestroySPrivateDatad ( GHoleDomaind *domain );
boolean _g1h_GetSplDBasisAuxpd ( GHoleDomaind *domain, int fn, int cn,
                                 int *nzc, double *fcomc, double *fcomcd );
boolean _g1h_GetSplDBasisCrossDerd ( GHoleDomaind *domain, int fn, int cn,
                                     double *fcomc, double *pv, double *pu );

boolean _g1h_TabBSFuncDer2d ( int deg, int lastknot, const double *knots,
                              int i0, int i1,
                              int n, const double *tkn, int *fkn, int *lkn,
                              double *b, double *bt, double *btt );

boolean _g1h_SetSplRightSided ( GHoleDomaind *domain,
                                int spdimen, CONST_ double *hole_cp,
                                const double *bmat,
                                double *fc00, double *b );
boolean _g1h_OutputSplPatchesd ( GHoleDomaind *domain,
                int spdimen, CONST_ double *x, double *fc00, void *usrptr,
                void (*outpatch) ( int n, int lknu, const double *knu,
                                   int m, int lknv, const double *knv,
                                   const double *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_Q2TabDiPatchJac3d ( int nkn, const double *kn, const double *hfunc,
             const double *dhfunc, const double *ddhfunc, const double *dddhfunc,
             const vector2d *c00, const vector2d *c01,
             const vector2d *c10, const vector2d *c11,
             const vector2d *d00, const vector2d *d01,
             const vector2d *d10, const vector2d *d11,
             double *jac, double *trd );
boolean _g1h_Q2TabLaplacianGradd ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *dddhfunc,
        const double *fc00, const double *fc01,
        const double *fc10, const double *fc11,
        const double *fd00, const double *fd01,
        const double *fd10, const double *fd11,
        const double *trd,
        vector2d *lapgrad );
boolean _g1h_Q2TabLaplacianGrad0d ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *dddhfunc,
        const double *fc00, const double *fc01,
        const double *fd00, const double *fd01,
        const double *trd,
        vector2d *lapgrad );
void _g1h_TabCurveJacobiand ( int deg, const point2d *cp,
                              int nkn, const double *kn, double *jac );
void _g1h_LapCoeffd ( const vector2d *du, const vector2d *dv,
                      const vector2d *duu, const vector2d *duv,
                      const vector2d *dvv, double *trd );
boolean _g1h_TabCurveLapCoeff0d ( const point2d *c00, const vector2d *c01,
                                  const point2d *c10, const vector2d *c11,
                                  const point2d *d00, const vector2d *d01,
                                  const point2d *d10, const vector2d *d11,
                                  int nkn, const double *tkn, const double *hfunc,
                                  const double *dhfunc, const double *ddhfunc,
                                  const double *atkn, const double *ahfunc,
                                  const double *adhfunc, const double *addhfunc,
                                  double *trdc00, double *trdc10,
                                  double *trdd00, double *trdd10 );
void _g1h_TabCurveLapCoeff1d ( const point2d *sicp, int nkn,
                               const double *tkn, double *trd );
boolean _g1h_TabTensBezPolyDer3d ( int nkn, const double *tkn,
             double *tbez, double *tbezu, double *tbezv,
             double *tbezuu, double *tbezuv, double *tbezvv,
             double *tbezuuu, double *tbezuuv, double *tbezuvv, double *tbezvvv );
boolean _g1h_Q2TabLaplacianJump0d ( int nkn, const double *tkn,
              const double *hfunc, const double *dhfunc, const double *ddhfunc,
              const double *atkn, const double *ahfunc, const double *adhfunc,
              const double *addhfunc,
              const double *ec00, const double *ec01,
              const double *ed00, const double *ed01, const double *etrdd00,
              const double *fc00, const double *fc01,
              const double *fd00, const double *fd01,
              const double *ftrdc00, const double *ftrdc10, const double *ftrdd10,
              double *lapjc00, double *lapjc10, double *lapjd10 );
boolean _g1h_Q2TabLaplacianJumpd ( int nkn, const double *tkn,
              const double *hfunc, const double *dhfunc, const double *ddhfunc,
              const double *atkn, const double *ahfunc, const double *adhfunc,
              const double *addhfunc,
              const double *ec00, const double *ec01,
              const double *ec10, const double *ec11,
              const double *ed00, const double *ed01,
              const double *ed10, const double *ed11, const double *etrdd00,
              const double *fc00, const double *fc01,
              const double *fc10, const double *fc11,
              const double *fd00, const double *fd01,
              const double *fd10, const double *fd11,
              const double *ftrdc00, const double *ftrdc10, const double *ftrdd10,
              const double *eicp1, const double *etrdc10,
              const double *eicp2, const double *etrdd10,
              double *lapjc00, double *lapjc10, double *lapjd10 );
unsigned short _g1h_ExtendSupport ( int hole_k, unsigned short supp );
boolean _g1h_TabLaplacianJump00d ( int nkn, const double *tkn, int fni,
                    const double *trdc00, const double *trdc10,
                    const double *trdd00, const double *trdd10,
                    double *lapc00, double *lapc10, double *lapd00, double *lapd10 );
double _g1h_Q2Integrald ( int hole_k, int nquad, double *jac,
                         unsigned short supp1, double *lapj1,
                         unsigned short supp2, double *lapj2 );
boolean _g1h_FuncDSuppd ( int hole_k, int nk, int m1, int fn, int i,
                          int *nzc, int *i0, int *i1, int *j0, int *j1 );

void g1h_splnloutpatchd ( int n, int lknu, const double *knu,
                          int m, int lknv, const double *knv,
                          const double *cp, void *usrptr );
boolean _g1h_TabBSFuncDer3d ( int deg, int lastknot, const double *knots,
                              int i0, int i1,
                              int n, const double *tkn, int *fkn, int *lkn,
                              double *b, double *bt, double *btt, double *bttt );
void _g1hq2_SetupCTrdd ( const vector2d *cdiu, const vector2d *cdiv,
           const vector2d *cdiuu, const vector2d *cdiuv, const vector2d *cdivv,
           double *ctrd );
boolean _g1hq2_FindDomSurrndPatchd ( GHoleDomaind *domain,
                                     G1HNLPrivated *nlpr,
                                     int i, int j, point2d *bezcp );
boolean _g1hq2_FindNLDomainDiameterd ( GHoleDomaind *domain,  
                                       G1HNLPrivated *nlpr );
boolean _g1h_TabBSFuncDer2Jd ( int rr, int deg, int nk, int m2,
             int lastcknot, const double *cknots,
             double *atbs, double *atbst, double *atbstt0, double *atbstt1 );

#endif

