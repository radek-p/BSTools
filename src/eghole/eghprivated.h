
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

#ifndef EGHPRIVATED_H
#define EGHPRIVATED_H

/* angle tolerances for partition analysis */
#define MIN_DELTA (PI/180.0)       /*   1 degree  */
#define MAX_DELTA (179.0*PI/180.0) /* 179 degrees */
#define ANGLETOL (0.001*PI/180)    /* 0.001 degree */

/* partition of the full angle description element */
typedef struct GHoleSgnPartd {
    double  alpha;  /* the first three fields of this structure */
    double  malpha; /* must be as they are, for the sake */
    double  salpha; /* of the sorting procedure */
    double  knot;
    boolean sgn, both;
  } GHoleSgnPartd;

typedef struct GHolePrivateRecd {
    point2d  *dombezpt;    /* Bezier patches surrounding the domain */
    vector2d *surrpc;      /* boundary curves and cross derivatives */
                           /* of the surrounding patches */
    vector2d *surrpcbc;    /* derivatives of the surrounding patches */
                           /* at relevant corners */
    double   *domsurrbf;   /* polynomials, which determine boundary */
                           /* conditions for the block B basis functions */
    double   diam;         /* domain "diameter" */
    unsigned char  *bfcpn;
    unsigned short *support_b;
  } GHolePrivateRecd;


int gh_GetDefaultOptiond ( GHoleDomaind *domain, int query, int qn,
                          int *ndata, int **idata, double **fdata );

double *_gh_GetKnotSequenced ( GHoleDomaind *domain, int i );
point2d *_gh_GetDomSurrndPatchd ( GHoleDomaind *domain, int i, int j );
boolean _gh_FindDomSurrndPatchd ( GHoleDomaind *domain,
                                  int i, int j, point2d *bezcp );
boolean _gh_FindDomSurrndBFuncPatchesd ( GHoleDomaind *domain );
boolean _gh_FindBasisFuncbSupport ( int hole_k, int nfunc_b,    
                                    unsigned short *support_b,
                                    unsigned char *bfcpn );
boolean _gh_AnalyzePartitiond ( GHoleDomaind *domain, int omcdeg,
                                vector2d *omcbc, int *_hole_m, int *_auxp_k,
                                double **_partition, GHoleSgnPartd **_spartition,
                                double *_spart_alpha0 );
void _gh_PrepareTabKnotsd ( int nquad, int opt, double *knots );
boolean _g2h_DiJacobian3d ( const vector2d *du, const vector2d *dv,
                            const vector2d *duu, const vector2d *duv,
                            const vector2d *dvv,
                            const vector2d *duuu, const vector2d *duuv,
                            const vector2d *duvv, const vector2d *dvvv,
                            double *jac, double *trd );
double _g2h_Integrald ( int hole_k, int nknsq, double *jac,
                        unsigned short supp1, vector2d *func1,
                        unsigned short supp2, vector2d *func2 );
double _g2h_Integral0d ( int nknsq, const double *jac,
                         const vector2d *func1, const vector2d *func2 );
double _g2h_SplIntegrald ( int nkn, double *jac,
                           int i00, int i01, int j00, int j01, vector2d *func0,
                           int i10, int i11, int j10, int j11, vector2d *func1 );
void _g2h_TabSplCLaplacianGradd ( int nkn,
            int fknp, int lknp, const double *tp, const double *tpu,
            const double *tpuu, const double *tpuuu,
            int fknq, int lknq, const double *tq, const double *tqv,
            const double *tqvv, const double *tqvvv,
            const double *trd, vector2d *lgr );
#endif

