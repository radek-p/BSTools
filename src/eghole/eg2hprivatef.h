
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

/* this header file is private for g2hole library functions; */
/* it is NOT intended to be #included in application source files */

#ifndef EG2HPRIVATEF_H
#define EG2HPRIVATEF_H

#ifndef EGHPRIVATEF_H
#include "eghprivatef.h"
#endif

/* degrees of polynomial junction functions */
#define G2_BF01DEG 2
#define G2_CG01DEG 1
#define G2_BF11DEG 4
#define G2_CG11DEG 3
#define G2_BF02DEG 3
#define G2_CG02DEG 2
#define G2_BF12DEG 5
#define G2_CG12DEG 4

/* degrees of patch cross derivatives */
#define G2_CROSS00DEG 7  /* == G2H_OMCDEG */
#define G2_CROSS01DEG 8
#define G2_CROSS02DEG 9
#define G2_CROSS10DEG 3
#define G2_CROSS11DEG 6
#define G2_CROSS12DEG 9

#define G2_CROSSDEGSUM \
        (G2_CROSS00DEG+G2_CROSS01DEG+G2_CROSS02DEG+ \
         G2_CROSS10DEG+G2_CROSS11DEG+G2_CROSS12DEG)

/* the following two numbers are relevant for spline bases */
#define G2_AUXDEG0 (G2_CROSS01DEG-G2H_OMCDEG+1)
#define G2_AUXDEG1 (G2_CROSS02DEG-G2H_OMCDEG+2)

/* the #definitions below are for the construction with */
/* the extended basis */
#define G2_DBDIM 16    /* diagonal block dimension, 4*4 */
#define G2_DIAGBLSIZE (G2_DBDIM*(G2_DBDIM+1)/2)


/* number of quadrature knots */
#define G2_NQUAD 16
#define G2_NQUADSQ (G2_NQUAD*G2_NQUAD)
#define G2_QUAD_FACTOR 10  /* for integration of spline basis functions */


/* ////////////////////////////////////////////////////////////////////////// */
/* private data of the construction procedures */
typedef struct G2HolePrivateRecf {
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
    unsigned char spdimen[4];   /* numbers of "inner" basis functions: */
                                /* "half polynomial" and "splines" of degree */
                                /* 3 and 4 */
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

    int (*GetOption)( GHoleDomainf *domain, int query, int qn,
                      int *ndata, int **idata, float **fdata );
  } G2HolePrivateRecf;


/* private data for the nonlinear constructions */
typedef struct G2HNLPrivatef {
    int      auxc;    /* patch counter */
    vector3f *nldi;   /* NL patches */
    vector3f *acoeff;
    vector3f *rhole_cp;
    vector3f nlnv, reflv;
    vector2f *diu, *div, *diuu, *diuv, *divv,
             *diuuu, *diuuv, *diuvv, *divvv;
    float    *jac;
    float    *psiu, *psiv, *psiuu, *psiuv, *psivv,
             *psiuuu, *psiuuv, *psiuvv, *psivvv;
  } G2HNLPrivatef;

typedef struct G2HNLFuncf {
    float    pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv, jac;
    float    psiu, psiv, psiuu, psiuv, psivv, psiuuu, psiuuv, psiuvv, psivvv,
             psju, psjv, psjuu, psjuv, psjvv, psjuuu, psjuuv, psjuvv, psjvvv;
    vector2f L, Lpu, Lpv, Lpuu, Lpuv, Lpvv, Lpuuu, Lpuuv, Lpuvv, Lpvvv,
             Lpupu, Lpupv, Lpupuu, Lpupuv, Lpupvv, Lpupuuu, Lpupuuv, Lpupuvv,
             Lpupvvv, Lpvpv, Lpvpuu, Lpvpuv, Lpvpvv, Lpvpuuu, Lpvpuuv,
             Lpvpuvv, Lpvpvvv, Lpuupuu, Lpuupuv, Lpuupvv, Lpuvpuv, Lpuvpvv,
             Lpvvpvv, BLT;
    float    D, Dpu, Dpv, Dpupu, Dpupv, Dpvpv;
    float    B[3], Bpu[3], Bpv[3], Bpupu[3], Bpupv[3], Bpvpv[3];
    float    LBLT;
  } G2HNLFuncf;


/* private data for the spline basis construction */
typedef struct G2HoleSPrivateRecf {
    int nk, m1, m2;   /* parameters, which determine the space dimension */
    int nsfunc_c, nsfunc_d;
                      /* numbers of basis functions of various kinds */
    int csize, dsize; /* sizes of the blocks in the form matrix */
    int   lastbezknot;
    float bezknots[2*(G2H_FINALDEG+1)];  /* knots for Bezier polynomials */
    int   lastomcknot;
    float *omcknots;  /* knots of basis function auxiliary patches */
    int lastpvknot, lastpvvknot;
    float *pvknots, *pvvknots;
                      /* knots of the basis functions of the D block */
    int lastcknot;
    float *cknots;    /* knots of the basis functions of the C block */
    float *basis_d;   /* Coons representations of basis functions */

    int lastfpknot;
    float *fpknots;   /* knots of the final patches */
    float *SAMat, *SBMat;  /* form matrices for the spline basis */
    float *SLMat;          /* A=LL^T decomposition factor */

    int    splnconstr; /* number of constraints */
    float  *SCmat, *SRCmat; /* matrices of constraint equations */
    int    splnaconstr, splacdim;
    float  *ASCmat, *ASRCmat;
  } G2HoleSPrivateRecf;


/* private data for spline nonlinear constructions */
typedef struct G2HNLSPrivatef {
    G2HNLPrivatef nlpr; /* must be the first field */
    int   nkn;          /* number of quadrature knots */
    float *tkn;         /* quadrature knots */
    int   ftabsize;     /* function values table size */
    int   psize;        /* number of control points of each final patch  */
    int   *fkn, *lkn;   /* ranges of knots for the functions of block C  */
    float *cb, *cbt, *cbtt, *cbttt;  /* B-spline functions and their derivatives */
    int   *cfuncvi;     /* indexes of the first samples of C block functions */
    int   *dfuncvi;     /* indexes of the first samples of D block functions */
  } G2HNLSPrivatef;


/* ////////////////////////////////////////////////////////////////////////// */
extern G2HNLPrivatef *_g2h_nlprivf;

/* ////////////////////////////////////////////////////////////////////////// */
void _g2h_GetDiPatchCurvesf ( GHoleDomainf *domain, int i,
                    point2f **c00, vector2f **c01, vector2f **c02,
                    point2f **c10, vector2f **c11, vector2f **c12,
                    point2f **d00, vector2f **d01, vector2f **d02,
                    point2f **d10, vector2f **d11, vector2f **d12 );
boolean _g2h_GetABasisAuxpf ( GHoleDomainf *domain, int fn,
                              float *br0, float *br0cr1, float *br0cr2 );
boolean _g2h_GetBBasisAuxpf ( GHoleDomainf *domain, int fn,
                              float *bezfc, float *fcomc,
                              float *fcomcd, float *fcomcdd );
void _g2h_GetBFAPatchCurvesf ( GHoleDomainf *domain, int fn, int i,
                    float **c00, float **c01, float **c02,
                    float **d00, float **d01, float **d02 );
void _g2h_GetBFBPatchCurvesf ( GHoleDomainf *domain, int fn, int i,
                    float **c00, float **c01, float **c02,
                    float **c10, float **c11, float **c12,
                    float **d00, float **d01, float **d02,
                    float **d10, float **d11, float **d12 );

boolean _g2h_TabDiPatchJac3f ( int nkn, const float *kn, const float *hfunc,
        const float *dhfunc, const float *ddhfunc, const float *dddhfunc,   
        const vector2f *c00, const vector2f *c01, const vector2f *c02,      
        const vector2f *c10, const vector2f *c11, const vector2f *c12,      
        const vector2f *d00, const vector2f *d01, const vector2f *d02,      
        const vector2f *d10, const vector2f *d11, const vector2f *d12,      
        float *jac, float *trd );
boolean _g2h_TabLaplacianGradf ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *dddhfunc,
        const float *fc00, const float *fc01, const float *fc02,
        const float *fc10, const float *fc11, const float *fc12,
        const float *fd00, const float *fd01, const float *fd02,
        const float *fd10, const float *fd11, const float *fd12,
        const float *trd,
        vector2f *lapgrad );
boolean _g2h_TabLaplacianGrad0f ( int nkn, const float *tkn,
        const float *hfunc, const float *dhfunc, const float *ddhfunc,
        const float *dddhfunc,
        const float *fc00, const float *fc01, const float *fc02,
        const float *fd00, const float *fd01, const float *fd02,
        const float *trd,
        vector2f *lapgrad );

boolean _g2h_VerifyJunctionFunctionsf ( GHoleDomainf *domain );
boolean _g2h_VerifyDomPatchesf ( GHoleDomainf *domain );
                                                    
/* ////////////////////////////////////////////////////////////////////////// */
void g2h_ReflectVectorsf ( int n, const vector3f *v, vector3f *w );
void g2h_nonlinoutpatchf ( int n, int m, const float *cp, void *usrptr );
boolean _g2h_StopItf ( int itn, float gn0, float gn,
                       float cn, float dcn, float scf );
boolean g2h_GetHoleSurrndPatchf ( GHoleDomainf *domain,       
                                  const point3f *hole_cp,
                                  int i, int j, point3f *bcp );
boolean _g2h_ComputeNLNormalf ( GHoleDomainf *domain,       
                                G2HNLPrivatef *nlprivate,       
                                const point3f *hole_cp );
boolean _g2h_TabNLDer0f ( GHoleDomainf *domain,
             int nkn, const float *tkn,       
             const float *hfunc, const float *dhfunc,
             const float *ddhfunc, const float *dddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             vector2f *diuuu, vector2f *diuuv, vector2f *diuvv, vector2f *divvv,
             float *fc00, float *fc01, float *fc02,
             float *fd00, float *fd01, float *fd02,
             float *psiu, float *psiv, 
             float *psiuu, float *psiuv, float *psivv,
             float *psiuuu, float *psiuuv, float *psiuvv, float *psivvv );
boolean _g2h_TabNLDerf ( GHoleDomainf *domain,
             int nkn, float *tkn,       
             const float *hfunc, const float *dhfunc,
             const float *ddhfunc, const float *dddhfunc,
             vector2f *diu, vector2f *div,
             vector2f *diuu, vector2f *diuv, vector2f *divv,
             vector2f *diuuu, vector2f *diuuv, vector2f *diuvv, vector2f *divvv,
             float *fc00, float *fc01, float *fc02,
             float *fc10, float *fc11, float *fc12,
             float *fd00, float *fd01, float *fd02,
             float *fd10, float *fd11, float *fd12,
             float *psiu, float *psiv, 
             float *psiuu, float *psiuv, float *psivv,
             float *psiuuu, float *psiuuv, float *psiuvv, float *psivvv );
boolean _g2h_TabNLBasisFunctionsf ( GHoleDomainf *domain,       
                                    G2HNLPrivatef *nlpr );
void _g2h_IntFunc1af ( G2HNLFuncf *f, float *funct );
void _g2h_IntFunc1bf ( G2HNLFuncf *f, float *funct );
void _g2h_IntFunc1cf ( G2HNLFuncf *f, float *funct );
void _g2h_IntFunc2bf ( G2HNLFuncf *f, float *grad );
void _g2h_IntFunc2cf ( G2HNLFuncf *f, vector2f *Li, float *Bi, vector2f *BiLT,       
                       float *Di, float *grad );
void _g2h_IntFunc3cf ( G2HNLFuncf *f, vector2f *Li, vector2f *Lj,       
                       vector2f *BiLT, vector2f *BjLT,       
                       float Di, float Dj, float *hessian );
 
/* ////////////////////////////////////////////////////////////////////////// */
/* the following horrible macros are used a number of times in the code. */
/* They assume that some variables are declared and have proper initial */
/* values */
#define G2GetPolynomialAddresses(jfunc,b01,c01,f01,g01,b11,c11,f11,g11,\
    b02,c02,f02,g02,b12,c12,f12,g12,b01b01,twob01c01,c01c01,\
    f01f01,twof01g01,g01g01,b11b11,twob11c11,c11c11,\
    f11f11,twof11g11,g11g11)\
  b01 = jfunc; \
  c01 = &b01[hole_k*(G2_BF01DEG+1)];  f01 = &c01[hole_k*(G2_CG01DEG+1)]; \
  g01 = &f01[hole_k*(G2_BF01DEG+1)];  b11 = &g01[hole_k*(G2_CG01DEG+1)]; \
  c11 = &b11[hole_k*(G2_BF11DEG+1)];  f11 = &c11[hole_k*(G2_CG11DEG+1)]; \
  g11 = &f11[hole_k*(G2_BF11DEG+1)];  b02 = &g11[hole_k*(G2_CG11DEG+1)]; \
  c02 = &b02[hole_k*(G2_BF02DEG+1)];  f02 = &c02[hole_k*(G2_CG02DEG+1)]; \
  g02 = &f02[hole_k*(G2_BF02DEG+1)];  b12 = &g02[hole_k*(G2_CG02DEG+1)]; \
  c12 = &b12[hole_k*(G2_BF12DEG+1)];  f12 = &c12[hole_k*(G2_CG12DEG+1)]; \
  g12 = &f12[hole_k*(G2_BF12DEG+1)]; \
  b01b01    = &g12[hole_k*(G2_CG12DEG+1)]; \
  twob01c01 = &b01b01[hole_k*(2*G2_BF01DEG+1)]; \
  c01c01    = &twob01c01[hole_k*(G2_BF01DEG+G2_CG01DEG+1)]; \
  f01f01    = &c01c01[hole_k*(2*G2_CG01DEG+1)]; \
  twof01g01 = &f01f01[hole_k*(2*G2_BF01DEG+1)]; \
  g01g01    = &twof01g01[hole_k*(G2_BF01DEG+G2_CG01DEG+1)]; \
  b11b11    = &g01g01[hole_k*(2*G2_CG01DEG+1)]; \
  twob11c11 = &b11b11[hole_k*(2*G2_BF11DEG+1)]; \
  c11c11    = &twob11c11[hole_k*(G2_BF11DEG+G2_CG11DEG+1)]; \
  f11f11    = &c11c11[hole_k*(2*G2_CG11DEG+1)]; \
  twof11g11 = &f11f11[hole_k*(2*G2_BF11DEG+1)]; \
  g11g11    = &twof11g11[hole_k*(G2_BF11DEG+G2_CG11DEG+1)]

#define G2GetPolynomialAddresses0(jfunc,b01,c01,f01,g01,b11,c11,f11,g11,\
    b02,c02,f02,g02,b12,c12,f12,g12,b01b01,twob01c01,c01c01,\
    f01f01,twof01g01,g01g01)\
  b01 = jfunc; \
  c01 = &b01[hole_k*(G2_BF01DEG+1)];  f01 = &c01[hole_k*(G2_CG01DEG+1)]; \
  g01 = &f01[hole_k*(G2_BF01DEG+1)];  b11 = &g01[hole_k*(G2_CG01DEG+1)]; \
  c11 = &b11[hole_k*(G2_BF11DEG+1)];  f11 = &c11[hole_k*(G2_CG11DEG+1)]; \
  g11 = &f11[hole_k*(G2_BF11DEG+1)];  b02 = &g11[hole_k*(G2_CG11DEG+1)]; \
  c02 = &b02[hole_k*(G2_BF02DEG+1)];  f02 = &c02[hole_k*(G2_CG02DEG+1)]; \
  g02 = &f02[hole_k*(G2_BF02DEG+1)];  b12 = &g02[hole_k*(G2_CG02DEG+1)]; \
  c12 = &b12[hole_k*(G2_BF12DEG+1)];  f12 = &c12[hole_k*(G2_CG12DEG+1)]; \
  g12 = &f12[hole_k*(G2_BF12DEG+1)]; \
  b01b01    = &g12[hole_k*(G2_CG12DEG+1)]; \
  twob01c01 = &b01b01[hole_k*(2*G2_BF01DEG+1)]; \
  c01c01    = &twob01c01[hole_k*(G2_BF01DEG+G2_CG01DEG+1)]; \
  f01f01    = &c01c01[hole_k*(2*G2_CG01DEG+1)]; \
  twof01g01 = &f01f01[hole_k*(2*G2_BF01DEG+1)]; \
  g01g01    = &twof01g01[hole_k*(G2_BF01DEG+G2_CG01DEG+1)]

#define G2GetPolyAddr(jfunc,b01,c01,f01,g01,b11,c11,f11,g11)\
  b01 = jfunc; \
  c01 = &b01[hole_k*(G2_BF01DEG+1)];  f01 = &c01[hole_k*(G2_CG01DEG+1)]; \
  g01 = &f01[hole_k*(G2_BF01DEG+1)];  b11 = &g01[hole_k*(G2_CG01DEG+1)]; \
  c11 = &b11[hole_k*(G2_BF11DEG+1)];  f11 = &c11[hole_k*(G2_CG11DEG+1)]; \
  g11 = &f11[hole_k*(G2_BF11DEG+1)]

#define G2GetDiCrossAddresses() \
  dir0cr1 = privateG2->dicross; \
  diq0cr1 = &dir0cr1[hole_k*(G2_CROSS01DEG+1)]; \
  dir1cr1 = &diq0cr1[hole_k*(G2_CROSS01DEG+1)]; \
  diq1cr1 = &dir1cr1[hole_k*(G2_CROSS11DEG+1)]; \
  dir0cr2 = &diq1cr1[hole_k*(G2_CROSS11DEG+1)]; \
  diq0cr2 = &dir0cr2[hole_k*(G2_CROSS02DEG+1)]; \
  dir1cr2 = &diq0cr2[hole_k*(G2_CROSS02DEG+1)]; \
  diq1cr2 = &dir1cr2[hole_k*(G2_CROSS12DEG+1)]

#define G2GetBFuncACrossAddresses() \
  bbr0 = privateG2->basis_a; \
  bbr0cr1 = &bbr0[nfunc_a*hole_k*(G2_CROSS00DEG+1)]; \
  bbq0cr1 = &bbr0cr1[nfunc_a*hole_k*(G2_CROSS01DEG+1)]; \
  bbr0cr2 = &bbq0cr1[nfunc_a*hole_k*(G2_CROSS01DEG+1)]; \
  bbq0cr2 = &bbr0cr2[nfunc_a*hole_k*(G2_CROSS02DEG+1)]

#define G2GetBFuncBCrossAddresses() \
  bbr0 = privateG2->basis_b; \
  bbr1 = &bbr0[nfunc_b*hole_k*(G2_CROSS00DEG+1)]; \
  bbq1 = &bbr1[nfunc_b*hole_k*(G2_CROSS10DEG+1)]; \
  bbr0cr1 = &bbq1[nfunc_b*hole_k*(G2_CROSS10DEG+1)]; \
  bbq0cr1 = &bbr0cr1[nfunc_b*hole_k*(G2_CROSS01DEG+1)]; \
  bbr1cr1 = &bbq0cr1[nfunc_b*hole_k*(G2_CROSS01DEG+1)]; \
  bbq1cr1 = &bbr1cr1[nfunc_b*hole_k*(G2_CROSS11DEG+1)]; \
  bbr0cr2 = &bbq1cr1[nfunc_b*hole_k*(G2_CROSS11DEG+1)]; \
  bbq0cr2 = &bbr0cr2[nfunc_b*hole_k*(G2_CROSS02DEG+1)]; \
  bbr1cr2 = &bbq0cr2[nfunc_b*hole_k*(G2_CROSS02DEG+1)]; \
  bbq1cr2 = &bbr1cr2[nfunc_b*hole_k*(G2_CROSS12DEG+1)]

#define G2GetFCAddresses() \
  fc01 = &fc00[hole_k*spdimen*(G2_CROSS00DEG+1)]; \
  fc02 = &fc01[hole_k*spdimen*(G2_CROSS01DEG+1)]; \
  fc10 = &fc02[hole_k*spdimen*(G2_CROSS02DEG+1)]; \
  fc11 = &fc10[hole_k*spdimen*(G2_CROSS10DEG+1)]; \
  fc12 = &fc11[hole_k*spdimen*(G2_CROSS11DEG+1)]; \
  fd00 = &fc12[hole_k*spdimen*(G2_CROSS12DEG+1)]; \
  fd01 = &fd00[hole_k*spdimen*(G2_CROSS00DEG+1)]; \
  fd02 = &fd01[hole_k*spdimen*(G2_CROSS01DEG+1)]; \
  fd10 = &fd02[hole_k*spdimen*(G2_CROSS02DEG+1)]; \
  fd11 = &fd10[hole_k*spdimen*(G2_CROSS10DEG+1)]; \
  fd12 = &fd11[hole_k*spdimen*(G2_CROSS11DEG+1)]

/* ///////////////////////////////////////////////////////////////////////// */
#define G2GetSFCAddresses() \
  sfc01 = &sfc00[hole_k*spdimen*(lastomcknot-G2_CROSS00DEG)]; \
  sfc02 = &sfc01[hole_k*spdimen*(lastpvknot-G2_CROSS01DEG)]; \
  sfd00 = &sfc02[hole_k*spdimen*(lastpvvknot-G2_CROSS02DEG)]; \
  sfd01 = &sfd00[hole_k*spdimen*(lastomcknot-G2_CROSS00DEG)]; \
  sfd02 = &sfd01[hole_k*spdimen*(lastpvknot-G2_CROSS01DEG)]

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2h_GetExtBlockAddressesf ( GHoleDomainf *domain,
                                     float **Aii, float **Aki, float **Akk,
                                     float **Bi, float **Bk, float **Lii );

boolean _g2h_SetRightSidef ( GHoleDomainf *domain,
                             int spdimen, CONST_ float *hole_cp,
                             float *fc00, float *b );
boolean _g2h_OutputPatchesf ( GHoleDomainf *domain,
                       int spdimen, CONST_ float *x, float *fc00, void *usrptr,
                       void (*outpatch) ( int n, int m, const float *cp,
                                          void *usrptr ) );

boolean _g2h_SetExtRightSidef ( GHoleDomainf *domain,
                                const float *Bi, const float *Bk,
                                int spdimen, CONST_ float *hole_cp,
                                float *fc00, float *b );
boolean _g2h_OutputExtPatchesf ( GHoleDomainf *domain,
                           int spdimen, CONST_ float *x, float *fc00, void *usrptr,
                           void (*outpatch) ( int n, int m, const float *cp,
                                              void *usrptr ) );

boolean _g2h_TabTensBezPolyDer3f ( int nkn, const float *tkn,
             float *tbez, float *tbezu, float *tbezv,
             float *tbezuu, float *tbezuv, float *tbezvv,
             float *tbezuuu, float *tbezuuv, float *tbezuvv, float *tbezvvv );

/* ///////////////////////////////////////////////////////////////////////// */
void _g2h_IntFunc1af ( G2HNLFuncf *f, float *funct );
void _g2h_IntFunc1bf ( G2HNLFuncf *f, float *funct );
void _g2h_IntFunc1cf ( G2HNLFuncf *f, float *funct );
void _g2h_IntFunc2bf ( G2HNLFuncf *f, float *grad );
void _g2h_IntFunc2cf ( G2HNLFuncf *f, vector2f *Li, float *Bi, vector2f *BiLT,
                       float *Di, float *grad );
void _g2h_IntFunc3cf ( G2HNLFuncf *f, vector2f *Li, vector2f *Lj,
                       vector2f *BiLT, vector2f *BjLT,
                       float Di, float Dj, float *hessian );

/* ///////////////////////////////////////////////////////////////////////// */
void g2h_DestroySPrivateDataf ( GHoleDomainf *domain );
boolean _g2h_GetSplDBasisAuxpf ( GHoleDomainf *domain, int fn, int cn,
                   int *nzc, float *fcomc, float *fcomcd, float *fcomcdd );
boolean _g2h_GetSplDBasisCrossDerf ( GHoleDomainf *domain, int fn, int cn,
                 float *fcomc, float *pv, float *pvv, float *pu, float *puu );

boolean _g2h_TabBSFuncDer3f ( int deg, int lastknot, const float *knots,
                              int i0, int i1,
                              int n, const float *tkn, int *fkn, int *lkn,
                              float *b, float *bt, float *btt, float *bttt );
void _g2h_TensDer3f ( float p, float pu, float puu, float puuu,
                      float q, float qv, float qvv, float qvvv,
                      float *pq );

boolean _g2h_FuncDSuppf ( int hole_k, int nk, int m1, int fn, int i,
                          int *nzc, int *i0, int *i1, int *j0, int *j1 );

boolean _g2h_SetSplRightSidef ( GHoleDomainf *domain,
                                int spdimen, CONST_ float *hole_cp,
                                float *fc00, float *b );
boolean _g2h_OutputSplPatchesf ( GHoleDomainf *domain,
                int spdimen, CONST_ float *x, float *fc00, void *usrptr,
                void (*outpatch) ( int n, int lknu, const float *knu,
                                   int m, int lknv, const float *knv,
                                   const float *cp, void *usrptr ) );

#endif

