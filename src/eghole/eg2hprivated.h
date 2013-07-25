
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

#ifndef EG2HPRIVATED_H
#define EG2HPRIVATED_H

#ifndef EGHPRIVATED_H
#include "eghprivated.h"
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
typedef struct G2HolePrivateRecd {
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
    unsigned char spdimen[4];   /* numbers of "inner" basis functions: */
                                /* "half polynomial" and "splines" of degree */
                                /* 3 and 4 */
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

    int (*GetOption)( GHoleDomaind *domain, int query, int qn,
                      int *ndata, int **idata, double **fdata );
  } G2HolePrivateRecd;


/* private data for the nonlinear constructions */
typedef struct G2HNLPrivated {
    int      auxc;    /* patch counter */
    vector3d *nldi;   /* NL patches */
    vector3d *acoeff;
    vector3d *rhole_cp;
    vector3d nlnv, reflv;
    vector2d *diu, *div, *diuu, *diuv, *divv,
             *diuuu, *diuuv, *diuvv, *divvv;
    double   *jac;
    double   *psiu, *psiv, *psiuu, *psiuv, *psivv,
             *psiuuu, *psiuuv, *psiuvv, *psivvv;
  } G2HNLPrivated;

typedef struct G2HNLFuncd {
    double   pu, pv, puu, puv, pvv, puuu, puuv, puvv, pvvv, jac;
    double   psiu, psiv, psiuu, psiuv, psivv, psiuuu, psiuuv, psiuvv, psivvv,
             psju, psjv, psjuu, psjuv, psjvv, psjuuu, psjuuv, psjuvv, psjvvv;
    vector2d L, Lpu, Lpv, Lpuu, Lpuv, Lpvv, Lpuuu, Lpuuv, Lpuvv, Lpvvv,
             Lpupu, Lpupv, Lpupuu, Lpupuv, Lpupvv, Lpupuuu, Lpupuuv, Lpupuvv,
             Lpupvvv, Lpvpv, Lpvpuu, Lpvpuv, Lpvpvv, Lpvpuuu, Lpvpuuv,
             Lpvpuvv, Lpvpvvv, Lpuupuu, Lpuupuv, Lpuupvv, Lpuvpuv, Lpuvpvv,
             Lpvvpvv, BLT;
    double   D, Dpu, Dpv, Dpupu, Dpupv, Dpvpv;
    double   B[3], Bpu[3], Bpv[3], Bpupu[3], Bpupv[3], Bpvpv[3];
    double   LBLT;
  } G2HNLFuncd;


/* private data for the spline basis construction */
typedef struct G2HoleSPrivateRecd {
    int nk, m1, m2;   /* parameters, which determine the space dimension */
    int nsfunc_c, nsfunc_d;
                      /* numbers of basis functions of various kinds */
    int csize, dsize; /* sizes of the blocks in the form matrix */
    int    lastbezknot;
    double bezknots[2*(G2H_FINALDEG+1)];  /* knots for Bezier polynomials */
    int    lastomcknot;
    double *omcknots;  /* knots of basis function auxiliary patches */
    int    lastpvknot, lastpvvknot;
    double *pvknots, *pvvknots;
                      /* knots of the basis functions of the D block */
    int lastcknot;
    double *cknots;    /* knots of the basis functions of the C block */
    double *basis_d;   /* Coons representations of basis functions */

    int lastfpknot;
    double *fpknots;   /* knots of the final patches */
    double *SAMat, *SBMat;  /* form matrices for the spline basis */
    double *SLMat;          /* A=LL^T decomposition factor */

    int    splnconstr; /* number of constraints */
    double *SCmat, *SRCmat; /* matrices of constraint equations */
    int    splnaconstr, splacdim;
    double *ASCmat, *ASRCmat;
  } G2HoleSPrivateRecd;


/* private data for spline nonlinear constructions */
typedef struct G2HNLSPrivated {
    G2HNLPrivated nlpr; /* must be the first field */
    int    nkn;         /* number of quadrature knots */
    double *tkn;        /* quadrature knots */
    int    ftabsize;    /* function values table size */
    int    psize;       /* number of control points of each final patch  */
    int    *fkn, *lkn;  /* ranges of knots for the functions of block C  */
    double *cb, *cbt, *cbtt, *cbttt;  /* B-spline functions and their derivatives */
    int    *cfuncvi;    /* indexes of the first samples of C block functions */
    int    *dfuncvi;    /* indexes of the first samples of D block functions */
  } G2HNLSPrivated;

 
/* ////////////////////////////////////////////////////////////////////////// */
extern G2HNLPrivated *_g2h_nlprivd;

/* ////////////////////////////////////////////////////////////////////////// */
void _g2h_GetDiPatchCurvesd ( GHoleDomaind *domain, int i,
                    point2d **c00, vector2d **c01, vector2d **c02,
                    point2d **c10, vector2d **c11, vector2d **c12,
                    point2d **d00, vector2d **d01, vector2d **d02,
                    point2d **d10, vector2d **d11, vector2d **d12 );
boolean _g2h_GetABasisAuxpd ( GHoleDomaind *domain, int fn,
                              double *br0, double *br0cr1, double *br0cr2 );
boolean _g2h_GetBBasisAuxpd ( GHoleDomaind *domain, int fn,
                              double *bezfc, double *fcomc,
                              double *fcomcd, double *fcomcdd );
void _g2h_GetBFAPatchCurvesd ( GHoleDomaind *domain, int fn, int i,
                    double **c00, double **c01, double **c02,
                    double **d00, double **d01, double **d02 );
void _g2h_GetBFBPatchCurvesd ( GHoleDomaind *domain, int fn, int i,
                    double **c00, double **c01, double **c02,
                    double **c10, double **c11, double **c12,
                    double **d00, double **d01, double **d02,
                    double **d10, double **d11, double **d12 );

boolean _g2h_TabDiPatchJac3d ( int nkn, const double *kn, const double *hfunc,
        const double *dhfunc, const double *ddhfunc, const double *dddhfunc,   
        const vector2d *c00, const vector2d *c01, const vector2d *c02,      
        const vector2d *c10, const vector2d *c11, const vector2d *c12,      
        const vector2d *d00, const vector2d *d01, const vector2d *d02,      
        const vector2d *d10, const vector2d *d11, const vector2d *d12,      
        double *jac, double *trd );
boolean _g2h_TabLaplacianGradd ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *dddhfunc,
        const double *fc00, const double *fc01, const double *fc02,
        const double *fc10, const double *fc11, const double *fc12,
        const double *fd00, const double *fd01, const double *fd02,
        const double *fd10, const double *fd11, const double *fd12,
        const double *trd,
        vector2d *lapgrad );
boolean _g2h_TabLaplacianGrad0d ( int nkn, const double *tkn,
        const double *hfunc, const double *dhfunc, const double *ddhfunc,
        const double *dddhfunc,
        const double *fc00, const double *fc01, const double *fc02,
        const double *fd00, const double *fd01, const double *fd02,
        const double *trd,
        vector2d *lapgrad );

boolean _g2h_VerifyJunctionFunctionsd ( GHoleDomaind *domain );
boolean _g2h_VerifyDomPatchesd ( GHoleDomaind *domain );
                                                    
/* ////////////////////////////////////////////////////////////////////////// */
void g2h_ReflectVectorsd ( int n, const vector3d *v, vector3d *w );
void g2h_nonlinoutpatchd ( int n, int m, const double *cp, void *usrptr );
boolean _g2h_StopItd ( int itn, double gn0, double gn,
                       double cn, double dcn, double scf );
boolean g2h_GetHoleSurrndPatchd ( GHoleDomaind *domain,       
                                  const point3d *hole_cp,
                                  int i, int j, point3d *bcp );
boolean _g2h_ComputeNLNormald ( GHoleDomaind *domain,       
                                G2HNLPrivated *nlprivate,       
                                const point3d *hole_cp );
boolean _g2h_TabNLDer0d ( GHoleDomaind *domain,
             int nkn, const double *tkn,       
             const double *hfunc, const double *dhfunc,
             const double *ddhfunc, const double *dddhfunc,
             vector2d *diu, vector2d *div,
             vector2d *diuu, vector2d *diuv, vector2d *divv,
             vector2d *diuuu, vector2d *diuuv, vector2d *diuvv, vector2d *divvv,
             double *fc00, double *fc01, double *fc02,
             double *fd00, double *fd01, double *fd02,
             double *psiu, double *psiv, 
             double *psiuu, double *psiuv, double *psivv,
             double *psiuuu, double *psiuuv, double *psiuvv, double *psivvv );
boolean _g2h_TabNLDerd ( GHoleDomaind *domain,
             int nkn, double *tkn,       
             const double *hfunc, const double *dhfunc,
             const double *ddhfunc, const double *dddhfunc,
             vector2d *diu, vector2d *div,
             vector2d *diuu, vector2d *diuv, vector2d *divv,
             vector2d *diuuu, vector2d *diuuv, vector2d *diuvv, vector2d *divvv,
             double *fc00, double *fc01, double *fc02,
             double *fc10, double *fc11, double *fc12,
             double *fd00, double *fd01, double *fd02,
             double *fd10, double *fd11, double *fd12,
             double *psiu, double *psiv, 
             double *psiuu, double *psiuv, double *psivv,
             double *psiuuu, double *psiuuv, double *psiuvv, double *psivvv );
boolean _g2h_TabNLBasisFunctionsd ( GHoleDomaind *domain,       
                                    G2HNLPrivated *nlpr );
void _g2h_IntFunc1ad ( G2HNLFuncd *f, double *funct );
void _g2h_IntFunc1bd ( G2HNLFuncd *f, double *funct );
void _g2h_IntFunc1cd ( G2HNLFuncd *f, double *funct );
void _g2h_IntFunc2bd ( G2HNLFuncd *f, double *grad );
void _g2h_IntFunc2cd ( G2HNLFuncd *f, vector2d *Li, double *Bi, vector2d *BiLT,       
                       double *Di, double *grad );
void _g2h_IntFunc3cd ( G2HNLFuncd *f, vector2d *Li, vector2d *Lj,       
                       vector2d *BiLT, vector2d *BjLT,       
                       double Di, double Dj, double *hessian );
 
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
boolean _g2h_GetExtBlockAddressesd ( GHoleDomaind *domain,
                                     double **Aii, double **Aki, double **Akk,
                                     double **Bi, double **Bk, double **Lii );

boolean _g2h_SetRightSided ( GHoleDomaind *domain,
                             int spdimen, CONST_ double *hole_cp,
                             double *fc00, double *b );
boolean _g2h_OutputPatchesd ( GHoleDomaind *domain,
                       int spdimen, CONST_ double *x, double *fc00, void *usrptr,
                       void (*outpatch) ( int n, int m, const double *cp,
                                          void *usrptr ) );

boolean _g2h_SetExtRightSided ( GHoleDomaind *domain,
                                const double *Bi, const double *Bk,
                                int spdimen, CONST_ double *hole_cp,
                                double *fc00, double *b );
boolean _g2h_OutputExtPatchesd ( GHoleDomaind *domain,
                           int spdimen, CONST_ double *x, double *fc00, void *usrptr,
                           void (*outpatch) ( int n, int m, const double *cp,
                                              void *usrptr ) );

boolean _g2h_TabTensBezPolyDer3d ( int nkn, const double *tkn,
             double *tbez, double *tbezu, double *tbezv,
             double *tbezuu, double *tbezuv, double *tbezvv,
             double *tbezuuu, double *tbezuuv, double *tbezuvv, double *tbezvvv );

/* ///////////////////////////////////////////////////////////////////////// */
void _g2h_IntFunc1ad ( G2HNLFuncd *f, double *funct );
void _g2h_IntFunc1bd ( G2HNLFuncd *f, double *funct );
void _g2h_IntFunc1cd ( G2HNLFuncd *f, double *funct );
void _g2h_IntFunc2bd ( G2HNLFuncd *f, double *grad ); 
void _g2h_IntFunc2cd ( G2HNLFuncd *f, vector2d *Li, double *Bi, vector2d *BiLT,
                       double *Di, double *grad );
void _g2h_IntFunc3cd ( G2HNLFuncd *f, vector2d *Li, vector2d *Lj,
                       vector2d *BiLT, vector2d *BjLT,
                       double Di, double Dj, double *hessian );

/* ///////////////////////////////////////////////////////////////////////// */
void g2h_DestroySPrivateDatad ( GHoleDomaind *domain );
boolean _g2h_GetSplDBasisAuxpd ( GHoleDomaind *domain, int fn, int cn,
                   int *nzc, double *fcomc, double *fcomcd, double *fcomcdd );
boolean _g2h_GetSplDBasisCrossDerd ( GHoleDomaind *domain, int fn, int cn,
                 double *fcomc, double *pv, double *pvv, double *pu, double *puu );

boolean _g2h_TabBSFuncDer3d ( int deg, int lastknot, const double *knots,
                              int i0, int i1,
                              int n, const double *tkn, int *fkn, int *lkn,
                              double *b, double *bt, double *btt, double *bttt );
void _g2h_TensDer3d ( double p, double pu, double puu, double puuu,
                      double q, double qv, double qvv, double qvvv,
                      double *pq );

boolean _g2h_FuncDSuppd ( int hole_k, int nk, int m1, int fn, int i,
                          int *nzc, int *i0, int *i1, int *j0, int *j1 );

boolean _g2h_SetSplRightSided ( GHoleDomaind *domain,
                                int spdimen, CONST_ double *hole_cp,
                                double *fc00, double *b );
boolean _g2h_OutputSplPatchesd ( GHoleDomaind *domain,
                int spdimen, CONST_ double *x, double *fc00, void *usrptr,
                void (*outpatch) ( int n, int lknu, const double *knu,
                                   int m, int lknv, const double *knv,
                                   const double *cp, void *usrptr ) );

#endif

