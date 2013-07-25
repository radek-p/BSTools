
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* this header file is private for g1hole library functions; */
/* it is NOT intended to be #included in application source files */

#ifndef EG1HPRIVATEF_H
#define EG1HPRIVATEF_H

#ifndef EGHPRIVATEF_H
#include "eghprivatef.h"
#endif

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
typedef struct G1HolePrivateRecf {
    int      hole_m;
    int      nfunc_a,   /* number of "inner" basis functions */
             nfunc_b;   /* number of "outer" basis functions */
    char     opt_cp, opt_der1, opt_basis, opt_quad;  /* recorded options */
    vector2f *omcbc;    /* boundary conditions for auxiliary patches */
    vector2f *omc;      /* auxiliary patches representation */
    float    *jfunc;    /* coefficients of junction functions */
    vector2f *dicross;  /* curves describing cross derivatives */
                        /* of the domain patches */
    float    *AuxiMat;  /* matrices of transformations of derivatives */
                        /* for composing functions with domain patches */
    float    *basis_a,  /* Coons representations of "inner" basis functions */
             *basis_b;  /* Coons represetnationf of "outer" basis functions */
    unsigned char spdimen[2];   /* numbers of "inner" basis functions: */
                                /* "half polynomial" and "splines" of degree 2 */
    float    *partition;        /* partition of the full angle */
    GHoleSgnPartf *spartition;  /* sorted partition, having hole_k-hole_m */
                                /* elements */
    float    spart_alpha0;      /* sorted partition bisectrix */
    float    *Amat, *Bmat;      /* matrices of scalar products */
                                /* of basis functions */
    float    *Lmat;             /* decomposition factor L:  A=LL^T */
    float    *EAmat, *EBmat;    /* matrices for the extended basis */
    float    *ELmat;            /* decomposition factor for the */
                                /* extended matrix */

    int      nconstr;           /* number of constraint equations */
    float    *Cmat, *RCmat;     /* matrices of constraint equations */
    int      naconstr, acdim;
    float    *ACmat, *ARCmat;
    int      extnconstr;        /* number of constraints and */
    float    *ECmat, *ERCmat;   /* matrices of constraint equations for */
                                /* the extended matrix */
    int      extnaconstr, extacdim;
    float    *AECmat, *AERCmat;

    float    *BBmat;            /* matrix of scalar products, for */
                                /* evaluation of the functional */

                              /* private data for the quasi G2 constructions */
    float     C1b;                 /* constant in the definition of the form */
    float     *Q2AMat, *Q2BMat;    /* form matrices */
    float     *Q2LMat;
    float     *Q2RCMat, *Q2ARCMat; /* constraints */

    float     C1e;                 /* constant in the definition of the form */
    float     *Q2EAMat, *Q2EBMat;  /* form matrices for the extended basis */
    float     *Q2ELMat;
    float     *Q2ERCMat, *Q2AERCMat; /* constraints */
    int (*GetOption)( GHoleDomainf *domain, int query, int qn,
                      int *ndata, int **idata, float **fdata );
  } G1HolePrivateRecf;


/* private data for the spline basis construction */
typedef struct G1HoleSPrivateRecf {
    int   nk, m1, m2;   /* parameters, which determine the space dimension */
    int   nsfunc_c, nsfunc_d;
                      /* numbers of basis functions of various kinds */
    int   csize, dsize; /* sizes of the blocks in the form matrix */
    int   lastbezknot;
    float bezknots[2*(G1H_FINALDEG+1)];  /* knots for Bezier polynomials */
    int   lastomcknot;
    float *omcknots;  /* knots of basis function auxiliary patches */
    int   lastpvknot;
    float *pvknots;
                      /* knots of the basis functions of the D block */
    int   lastcknot;
    float *cknots;    /* knots of the basis functions of the C block */
    float *basis_d;   /* Coons representations of basis functions */   

    int   lastfpknot;
    float *fpknots;   /* knots of the final patches */
    float *SAMat, *SBMat;  /* form matrices for the spline basis */
    float *SLMat;          /* A=LL^T decomposition factor */

    int    splnconstr; /* number of constraints */
    float  *SCmat, *SRCmat; /* matrices of constraint equations */
    int    splnaconstr, splacdim;
    float  *ASCmat, *ASRCmat;

/* private data for the quasi G2 constructions */
    float C1s;                 /* constant in the definition of the form */
    float *Q2SAMat, *Q2SBMat;  /* form matrices for the spline basis */
    float *Q2SLMat;

    float  *Q2SRCmat; /* matrices of constraint equations */
    float  *Q2SARCmat;
  } G1HoleSPrivateRecf;


/* private data for the nonlinear constructions */
typedef struct G1HNLPrivatef {
    int      auxc;    /* patch counter */
    vector3f *nldi;   /* NL patches */
    vector3f *acoeff;
    vector3f *rhole_cp;
    vector3f nlnv, reflv;
    vector2f *diu, *div, *diuu, *diuv, *divv;
    float    ddiam;   /* NL domain diameter */
    float    *jac;
    float    *psiu, *psiv, *psiuu, *psiuv, *psivv;
                   /* the following data are used only for Q2 */
    vector2f *diuuu, *diuuv, *diuvv, *divvv;
    float    *psiuuu, *psiuuv, *psiuvv, *psivvv;
    vector2f *ctang;
    float    *cpsiu, *cpsiv, *cpsiuu, *cpsiuv, *cpsivv;
  } G1HNLPrivatef;  


typedef struct G1Q2HNLFuncf {
    float    pu, pv, jpuu, jpuv, jpvv;
    vector2f tang;
    float    psiu, psiv, jpsiuu, jpsiuv, jpsivv;
    float    psju, psjv, jpsjuu, jpsjuv, jpsjvv;
    float    A, B, e0, e1, e2, e9, tGt, b1, b2, b3, b4, b5;
  } G1Q2HNLFuncf;


/* private data for spline nonlinear constructions */
typedef struct G1HNLSPrivatef {
    G1HNLPrivatef nlpr; /* must be the first field */
    int     nkn;          /* number of quadrature knots */
    float   *tkn;         /* quadrature knots */
    int     ftabsize;     /* function values table size */
    int     psize;        /* number of control points of each final patch  */
    int     *fkn, *lkn;   /* ranges of knots for the functions of block C  */
    float   *cb, *cbt, *cbtt, *cbttt; /* B-spline functions and their derivatives */
    int     *cfuncvi;     /* indexes of the first samples of C block functions */
    int     *dfuncvi;     /* indexes of the first samples of D block functions */
    int     jtabsize;     /* jump values table size */
    int     jcfs, jdfs;   /* indexes of the first samples of jump */
                          /* for block C and block D functions */
    short   njcurves;
    boolean jumpC, jumpD;
  } G1HNLSPrivatef;


/* ////////////////////////////////////////////////////////////////////////// */
extern G1HNLPrivatef *_g1h_nlprivf;

/* ////////////////////////////////////////////////////////////////////////// */
void _g1h_GetDiPatchCurvesf ( GHoleDomainf *domain, int i,
                point2f **c00, vector2f **c01, point2f **c10, vector2f **c11,
                point2f **d00, vector2f **d01, point2f **d10, vector2f **d11 );
boolean _g1h_GetABasisAuxpf ( GHoleDomainf *domain, int fn,
                              float *br0, float *br0cr1 );
boolean _g1h_GetBBasisAuxpf ( GHoleDomainf *domain, int fn,
                              float *bezfc, float *fcomc, float *fcomcd );
void _g1h_GetBFAPatchCurvesf ( GHoleDomainf *domain, int fn, int i,
                    float **c00, float **c01,
                    float **d00, float **d01 );
void _g1h_GetBFBPatchCurvesf ( GHoleDomainf *domain, int fn, int i,
                    float **c00, float **c01, float **c10, float **c11,
                    float **d00, float **d01, float **d10, float **d11 );

void g1h_TabCubicHFuncDer2f ( int nkn, const float *kn,
                              float *hfunc, float *dhfunc, float *ddhfunc );
boolean g1h_TabBicubicCoonsPatchDer2f (
      int spdimen, int nkn, const float *kn, const float *hfunc,
      const float *dhfunc, const float *ddhfunc,
      int degc00, const float *c00,
      int degc01, const float *c01,
      int degc10, const float *c10,
      int degc11, const float *c11,
      int degd00, const float *d00,
      int degd01, const float *d01,
      int degd10, const float *d10,
      int degd11, const float *d11,
      float *p, float *pu, float *pv, float *puu, float *puv, float *pvv );
void _g1h_DiJacobian2f ( const vector2f *du, const vector2f *dv,
                         const vector2f *duu, const vector2f *duv,
                         const vector2f *dvv,
                         float *jac, float *trd );
boolean _g1h_TabDiPatchJac2f ( int nkn, const float *kn, const float *hfunc,
                               const float *dhfunc, const float *ddhfunc,
                               const vector2f *c00, const vector2f *c01,
                               const vector2f *c10, const vector2f *c11,
                               const vector2f *d00, const vector2f *d01,
                               const vector2f *d10, const vector2f *d11,
                               float *jac, float *trd );
boolean _g1h_TabLaplacianf ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *fc00, const float *fc01,
        const float *fc10, const float *fc11,
        const float *fd00, const float *fd01,
        const float *fd10, const float *fd11,
        const float *trd,
        float *lap );
boolean _g1h_TabLaplacian0f ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *fc00, const float *fc01,
        const float *fd00, const float *fd01,
        const float *trd,
        float *lapgrad );
float _g1h_Integralf ( int hole_k, int nknsq, float *jac,
                       unsigned short supp1, float *func1,
                       unsigned short supp2, float *func2 );

boolean _g1h_VerifyJunctionFunctionsf ( GHoleDomainf *domain );
boolean _g1h_VerifyDomPatchesf ( GHoleDomainf *domain );
                                                    
/* ////////////////////////////////////////////////////////////////////////// */
void g1h_ReflectVectorsf ( int n, const vector3f *v, vector3f *w );
void g1h_nonlinoutpatchf ( int n, int m, const float *cp, void *usrptr );
boolean _g1h_StopItf ( int itn, float gn0, float gn,
                       float cn, float dcn, float scf );
boolean g1h_GetHoleSurrndPatchf ( GHoleDomainf *domain,
                                  const point3f *hole_cp,
                                  int i, int j, point3f *bcp );
boolean _g1h_ComputeNLNormalf ( GHoleDomainf *domain,
                                G1HNLPrivatef *nlprivate,
                                const point3f *hole_cp );
boolean _g1h_TabNLDer0f ( int nkn, const float *tkn,
             const float *hfunc, const float *dhfunc, const float *ddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             float *fc00, float *fc01, float *fd00, float *fd01,
             float *psiu, float *psiv,
             float *psiuu, float *psiuv, float *psivv );
boolean _g1h_TabNLDerf ( int nkn, float *tkn,
             const float *hfunc, const float *dhfunc, const float *ddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             float *fc00, float *fc01, float *fc10, float *fc11,
             float *fd00, float *fd01, float *fd10, float *fd11,
             float *psiu, float *psiv,
             float *psiuu, float *psiuv, float *psivv );
boolean _g1h_TabNLBasisFunctionsf ( GHoleDomainf *domain, int nkn,
                                    G1HNLPrivatef *nlpr );
boolean _g1h_ReflectSplAltConstrMatrixf ( GHoleDomainf *domain,
                                          vector3f *reflv, float *RACmat );
void _g1h_IntFunc1f ( float pu, float pv,
                      float puu, float puv, float pvv, float jac,
                      float *c, float *cs, float *a, float *b,   
                      float *funct );
void _g1h_IntFunc2f ( float pu, float pv,
                      float puu, float puv, float pvv,
                      float psiu, float psiv,
                      float psiuu, float psiuv, float psivv,
                      float jac, float c, float A, float B, 
                      float *ai, float *bi, float *grad );
void _g1h_IntFunc3f ( float pu, float pv,       
                      float puu, float puv, float pvv,       
                      float psiu, float psiv,       
                      float psiuu, float psiuv, float psivv,       
                      float psju, float psjv,       
                      float psjuu, float psjuv, float psjvv,       
                      float jac, float c,       
                      float A, float B, float Ai, float Bi,       
                      float Aj, float Bj,       
                      float *hessian );
 
/* ////////////////////////////////////////////////////////////////////////// */
boolean _g1hq2_TabNLDer0f ( int nkn, const float *tkn,
             const float *hfunc, const float *dhfunc, const float *ddhfunc,
             const float *dddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             vector2f *diuuu, vector2f *diuuv, vector2f *diuvv, vector2f *divvv,
             float *fc00, float *fc01, float *fd00, float *fd01,
             float *psiu, float *psiv,
             float *psiuu, float *psiuv, float *psivv,
             float *psiuuu, float *psiuuv, float *psiuvv, float *psivvv );

boolean _g1hq2_TabNLDerf ( int nkn, float *tkn,
             const float *hfunc, const float *dhfunc, const float *ddhfunc,
             const float *dddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             vector2f *diuuu, vector2f *diuuv, vector2f *diuvv, vector2f *divvv,
             float *fc00, float *fc01, float *fc10, float *fc11,
             float *fd00, float *fd01, float *fd10, float *fd11,
             float *psiu, float *psiv,
             float *psiuu, float *psiuv, float *psivv,
             float *psiuuu, float *psiuuv, float *psiuvv, float *psivvv );

boolean _g1hq2_TabNLBasisFunctionsOmegaf ( GHoleDomainf *domain, int nkn,
                                           G1HNLPrivatef *nlpr, float *bc00 );
boolean _g1hq2_TabNLBasisFunctionsGammaf ( GHoleDomainf *domain, int nkn,
                                           G1HNLPrivatef *nlpr, float *ctrd,
                                           float *bc00 );
void _g1hq2_IntFunc1af ( G1Q2HNLFuncf *f, float *funct );
void _g1hq2_IntFunc1bf ( G1Q2HNLFuncf *f, float *funct );
void _g1hq2_IntFunc1cf ( G1Q2HNLFuncf *f, float *funct );
void _g1hq2_IntFunc2bf ( G1Q2HNLFuncf *f, float *grad );
void _g1hq2_IntFunc2cf ( G1Q2HNLFuncf *f, float *Ai, float *Bi, float *grad );
void _g1hq2_IntFunc3cf ( G1Q2HNLFuncf *f, float Ai, float Bi, float Aj, float Bj,
                         float *hessian );

boolean _g1h_ReflectAltConstrMatrixf ( GHoleDomainf *domain,
                                       vector3f *reflv, float *RACmat );
boolean _g1h_ReflectExtAltConstrMatrixf ( GHoleDomainf *domain,
                                          vector3f *reflv, float *RACmat );

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
boolean _g1h_GetExtBlockAddressesf ( GHoleDomainf *domain,
                                     float **Aii, float **Aki, float **Akk,
                                     float **Bi, float **Bk, float **Lii );

boolean _g1h_SetRightSidef ( GHoleDomainf *domain, const float *Bmat,
                             int spdimen, CONST_ float *hole_cp,
                             float *fc00, float *b );
boolean _g1h_OutputPatchesf ( GHoleDomainf *domain, int spdimen,
                       CONST_ float *x, float *fc00, void *usrptr,
                       void (*outpatch) ( int n, int m, const float *cp,
                                          void *usrptr ) );

boolean _g1h_SetExtRightSidef ( GHoleDomainf *domain,
                                const float *Bi, const float *Bk,
                                int spdimen, CONST_ float *hole_cp,
                                float *fc00, float *b );
boolean _g1h_OutputExtPatchesf ( GHoleDomainf *domain, int spdimen,
                           CONST_ float *x, float *fc00, void *usrptr,
                           void (*outpatch) ( int n, int m, const float *cp,
                                              void *usrptr ) );

boolean _g1h_TabTensBezPolyDer2f ( int nkn, const float *tkn,       
                     float *tbez, float *tbezu, float *tbezv,
                     float *tbezuu, float *tbezuv, float *tbezvv );

void _g1h_TensDer2f ( float p, float pu, float puu,
                      float q, float qv, float qvv,
                      float *pq );

/* ///////////////////////////////////////////////////////////////////////// */
void g1h_DestroySPrivateDataf ( GHoleDomainf *domain );
boolean _g1h_GetSplDBasisAuxpf ( GHoleDomainf *domain, int fn, int cn,
                                 int *nzc, float *fcomc, float *fcomcd );
boolean _g1h_GetSplDBasisCrossDerf ( GHoleDomainf *domain, int fn, int cn,
                                     float *fcomc, float *pv, float *pu );

boolean _g1h_TabBSFuncDer2f ( int deg, int lastknot, const float *knots,
                              int i0, int i1,
                              int n, const float *tkn, int *fkn, int *lkn,
                              float *b, float *bt, float *btt );

boolean _g1h_SetSplRightSidef ( GHoleDomainf *domain,
                                int spdimen, CONST_ float *hole_cp,
                                const float *bmat,
                                float *fc00, float *b );
boolean _g1h_OutputSplPatchesf ( GHoleDomainf *domain,
                int spdimen, CONST_ float *x, float *fc00, void *usrptr,
                void (*outpatch) ( int n, int lknu, const float *knu,
                                   int m, int lknv, const float *knv,
                                   const float *cp, void *usrptr ) );

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g1h_Q2TabDiPatchJac3f ( int nkn, const float *kn, const float *hfunc,
             const float *dhfunc, const float *ddhfunc, const float *dddhfunc,
             const vector2f *c00, const vector2f *c01,
             const vector2f *c10, const vector2f *c11,
             const vector2f *d00, const vector2f *d01,
             const vector2f *d10, const vector2f *d11,
             float *jac, float *trd );
boolean _g1h_Q2TabLaplacianGradf ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *dddhfunc,
        const float *fc00, const float *fc01,
        const float *fc10, const float *fc11,
        const float *fd00, const float *fd01,
        const float *fd10, const float *fd11,
        const float *trd,
        vector2f *lapgrad );
boolean _g1h_Q2TabLaplacianGrad0f ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *dddhfunc,
        const float *fc00, const float *fc01,
        const float *fd00, const float *fd01,
        const float *trd,
        vector2f *lapgrad );
void _g1h_TabCurveJacobianf ( int deg, const point2f *cp,
                              int nkn, const float *kn, float *jac );
void _g1h_LapCoefff ( const vector2f *du, const vector2f *dv,
                      const vector2f *duu, const vector2f *duv,
                      const vector2f *dvv, float *trd );
boolean _g1h_TabCurveLapCoeff0f ( const point2f *c00, const vector2f *c01,
                                  const point2f *c10, const vector2f *c11,
                                  const point2f *d00, const vector2f *d01,
                                  const point2f *d10, const vector2f *d11,
                                  int nkn, const float *tkn, const float *hfunc,
                                  const float *dhfunc, const float *ddhfunc,
                                  const float *atkn, const float *ahfunc,
                                  const float *adhfunc, const float *addhfunc,
                                  float *trdc00, float *trdc10,
                                  float *trdd00, float *trdd10 );
void _g1h_TabCurveLapCoeff1f ( const point2f *sicp, int nkn,
                               const float *tkn, float *trd );
boolean _g1h_TabTensBezPolyDer3f ( int nkn, const float *tkn,
             float *tbez, float *tbezu, float *tbezv,
             float *tbezuu, float *tbezuv, float *tbezvv,
             float *tbezuuu, float *tbezuuv, float *tbezuvv, float *tbezvvv );
boolean _g1h_Q2TabLaplacianJump0f ( int nkn, const float *tkn,
              const float *hfunc, const float *dhfunc, const float *ddhfunc,
              const float *atkn, const float *ahfunc, const float *adhfunc,
              const float *addhfunc,
              const float *ec00, const float *ec01,
              const float *ed00, const float *ed01, const float *etrdd00,
              const float *fc00, const float *fc01,
              const float *fd00, const float *fd01,
              const float *ftrdc00, const float *ftrdc10, const float *ftrdd10,
              float *lapjc00, float *lapjc10, float *lapjd10 );
boolean _g1h_Q2TabLaplacianJumpf ( int nkn, const float *tkn,
              const float *hfunc, const float *dhfunc, const float *ddhfunc,
              const float *atkn, const float *ahfunc, const float *adhfunc,
              const float *addhfunc,
              const float *ec00, const float *ec01,
              const float *ec10, const float *ec11,
              const float *ed00, const float *ed01,
              const float *ed10, const float *ed11, const float *etrdd00,
              const float *fc00, const float *fc01,
              const float *fc10, const float *fc11,
              const float *fd00, const float *fd01,
              const float *fd10, const float *fd11,
              const float *ftrdc00, const float *ftrdc10, const float *ftrdd10,
              const float *eicp1, const float *etrdc10,
              const float *eicp2, const float *etrdd10,
              float *lapjc00, float *lapjc10, float *lapjd10 );
unsigned short _g1h_ExtendSupport ( int hole_k, unsigned short supp );
boolean _g1h_TabLaplacianJump00f ( int nkn, const float *tkn, int fni,
                    const float *trdc00, const float *trdc10,
                    const float *trdd00, const float *trdd10,
                    float *lapc00, float *lapc10, float *lapd00, float *lapd10 );
float _g1h_Q2Integralf ( int hole_k, int nquad, float *jac,
                         unsigned short supp1, float *lapj1,
                         unsigned short supp2, float *lapj2 );
boolean _g1h_FuncDSuppf ( int hole_k, int nk, int m1, int fn, int i,
                          int *nzc, int *i0, int *i1, int *j0, int *j1 );

void g1h_splnloutpatchf ( int n, int lknu, const float *knu,
                          int m, int lknv, const float *knv,
                          const float *cp, void *usrptr );
boolean _g1h_TabBSFuncDer3f ( int deg, int lastknot, const float *knots,
                              int i0, int i1,
                              int n, const float *tkn, int *fkn, int *lkn,
                              float *b, float *bt, float *btt, float *bttt );
void _g1hq2_SetupCTrdf ( const vector2f *cdiu, const vector2f *cdiv,
           const vector2f *cdiuu, const vector2f *cdiuv, const vector2f *cdivv,
           float *ctrd );
boolean _g1hq2_FindDomSurrndPatchf ( GHoleDomainf *domain,
                                     G1HNLPrivatef *nlpr,
                                     int i, int j, point2f *bezcp );
boolean _g1hq2_FindNLDomainDiameterf ( GHoleDomainf *domain,
                                       G1HNLPrivatef *nlpr );
boolean _g1h_TabBSFuncDer2Jf ( int rr, int deg, int nk, int m2,
             int lastcknot, const float *cknots,
             float *atbs, float *atbst, float *atbstt0, float *atbstt1 );

#endif

