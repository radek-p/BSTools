
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

#ifndef EGHPRIVATEF_H
#define EGHPRIVATEF_H

/* angle tolerances for partition analysis */
#define MIN_DELTA (PI/180.0)       /*   1 degree  */
#define MAX_DELTA (179.0*PI/180.0) /* 179 degrees */
#define ANGLETOL (0.001*PI/180)    /* 0.001 degree */

/* partition of the full angle description element */
typedef struct GHoleSgnPartf {
    float   alpha;  /* the first three fields of this structure */
    float   malpha; /* must be as they are, for the sake */
    float   salpha; /* of the sorting procedure */
    float   knot;
    boolean sgn, both;
  } GHoleSgnPartf;

typedef struct GHolePrivateRecf {
    point2f  *dombezpt;    /* Bezier patches surrounding the domain */
    vector2f *surrpc;      /* boundary curves and cross derivatives */
                           /* of the surrounding patches */
    vector2f *surrpcbc;    /* derivatives of the surrounding patches */
                           /* at relevant corners */
    float    *domsurrbf;   /* polynomials, which determine boundary */
                           /* conditions for the block B basis functions */
    float    diam;         /* domain "diameter" */
    unsigned char  *bfcpn;
    unsigned short *support_b;
  } GHolePrivateRecf;


int gh_GetDefaultOptionf ( GHoleDomainf *domain, int query, int qn,
                          int *ndata, int **idata, float **fdata );

float *_gh_GetKnotSequencef ( GHoleDomainf *domain, int i );
point2f *_gh_GetDomSurrndPatchf ( GHoleDomainf *domain, int i, int j );
boolean _gh_FindDomSurrndPatchf ( GHoleDomainf *domain,
                                  int i, int j, point2f *bezcp );
boolean _gh_FindBasisFuncbSupport ( int hole_k, int nfunc_b,  
                                    unsigned short *support_b,
                                    unsigned char *bfcpn );
boolean _gh_FindDomSurrndBFuncPatchesf ( GHoleDomainf *domain );
boolean _gh_AnalyzePartitionf ( GHoleDomainf *domain, int omcdeg,
                                vector2f *omcbc, int *_hole_m, int *_auxp_k,
                                float **_partition, GHoleSgnPartf **_spartition,
                                float *_spart_alpha0 );
void _gh_PrepareTabKnotsf ( int nquad, int opt, float *knots );
boolean _g2h_DiJacobian3f ( const vector2f *du, const vector2f *dv,
                            const vector2f *duu, const vector2f *duv,
                            const vector2f *dvv,
                            const vector2f *duuu, const vector2f *duuv,
                            const vector2f *duvv, const vector2f *dvvv,
                            float *jac, float *trd );
float _g2h_Integralf ( int hole_k, int nknsq, float *jac,
                       unsigned short supp1, vector2f *func1,
                       unsigned short supp2, vector2f *func2 );
float _g2h_Integral0f ( int nknsq, const float *jac,
                        const vector2f *func1, const vector2f *func2 );
float _g2h_SplIntegralf ( int nkn, float *jac,
                          int i00, int i01, int j00, int j01, vector2f *func0,
                          int i10, int i11, int j10, int j11, vector2f *func1 );
void _g2h_TabSplCLaplacianGradf ( int nkn,
            int fknp, int lknp, const float *tp, const float *tpu,
            const float *tpuu, const float *tpuuu,
            int fknq, int lknq, const float *tq, const float *tqv,
            const float *tqvv, const float *tqvvv,
            const float *trd, vector2f *lgr );

#endif

