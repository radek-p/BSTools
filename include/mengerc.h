
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#ifndef MENGERC_H
#define MENGERC_H

#ifndef PKVARIA_H
#include "pkvaria.h"
#endif
#ifndef PKNUM_H
#include "pknum.h"
#endif
#ifndef PKGEOM_H
#include "pkgeom.h"
#endif
#ifndef MULTIBS_H
#include "multibs.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* number of quadrature knots; 4 seems enough */
/* the limits below are related with the code in mengerc00.c */
#define MENGERC_MIN_NQKN  2
#define MENGERC_MAX_NQKN 10

/* number of penalty parameters */
#define MENGERC_NPPARAM 5

/* penalty parameters optimization options */
#define MENGERC_OPT_NONE  0
#define MENGERC_OPT_FULL1 1
#define MENGERC_OPT_FULL2 2
#define MENGERC_OPT_PART  3

/* ///////////////////////////////////////////////////////////////////////// */
typedef struct {
          /* quadrature data */
    int     nqkn;     /* number of quadrature knots */
    double  *qkn;     /* pointer to the array of knots */
    double  *qc;      /* quadrature coefficients */
          /* spline basis description */
    int     deg;      /* degree, should be 3 */
    double  *bsf;     /* basis function values at knots */
    double  *dbsf;    /* 1st order derivatives at knots  */
    double  *ddbsf;   /* 2nd order derivatives at knots */
    double  *bsf1;    /* basis function values of degree n-1 */
    double  *dbsf1;   /* 1st order derivatives of basis functions of degree n-1 */
          /* curve description */
    int     lkn;      /* last knot number */
    double  *knots;   /* knot sequence - consecutive integers from 0 to lkn */
    point3d *cpoints; /* control points */
    point3d sc;       /* gravity centre of the control points */
          /* Menger curvature exponent */
    double  w;
          /* desired length of the curve */
    double  L;
          /* penalty parameters */
    double  *penalty_param;
    int     mdi;  /* control point most distant from the gravity centre */
    boolean alt_scale;
          /* pre-transformation switch */
    boolean pretransf;
          /* debugging output */
    boolean k0, k1, k2, k3, k4;
    double  k1min, k1max, k2min, k2max, k3min, k3max;
          /* number of parallel threads */
    int     npthr;
          /* working space */
    int     nvars;
    double  *x, *g;
          /* terms of the function, gradient and Hessian */
            /* recorded by the procedure evaluating the function */
    double  ffkM, ffkMe, lgt, ffR[MENGERC_NPPARAM];
    double  *fx;
            /* recorded by the procedure evaluating the gradient */
    double  gfkMe, gfR1, gfR2, gfR3, gfR4, gfR5;
    double  *ggkMe, *ggR1, *ggR2, *ggR3, *ggR4, *ggR5;
    double  *gx;
            /* recorded by the procedure evaluating the Hessian */
    double  hfkMe, hfR1, hfR2, hfR3, hfR4, hfR5;
    double  *hgkMe, *hgR1, *hgR2, *hgR3, *hgR4, *hgR5;
    double  *hhkMe, *hhR1, *hhR2, *hhR3, *hhR4, *hhR5;
    double  *hx;
            /* extreme eigenvalues */
    boolean heigok;
    double  hmin, hmax, _emin, _emax, emin, emax, ef, f;
            /* other data for optimization of the penalty parameters */
    int     pp_opt;
    double  *mcpoints, fmin;
    boolean ppopt;
            /* Levenberg-Marquardt trajectory parameter */
    double  nu;
    int     fcnt, itc, itres;
    double  lastf, gn;
  } mengerc_data;

/* ///////////////////////////////////////////////////////////////////////// */
boolean mengerc_TabBasisFunctions ( int deg, int nqkn, mengerc_data *md );
boolean mengerc_BindACurve ( mengerc_data *md,
                             int deg, int lkn, double *knots, point3d *cpoints,
                             int nqkn, double w, double *penalty_param,
                             boolean alt_scale );
void mengerc_UntieTheCurve ( mengerc_data *md );

/* ///////////////////////////////////////////////////////////////////////// */
boolean mengerc_intF ( mengerc_data *md,
                       int lkn, double *knots, point3d *cpoints,
                       double *func );
boolean mengerc_gradIntF ( mengerc_data *md,
                           int lkn, double *knots, point3d *cpoints,
                           double *intf, double *grad );
boolean mengerc_hessIntF ( mengerc_data *md,
                           int lkn, double *knots, point3d *cpoints,
                           double *intf, double *grad, double *hess );

boolean _mengerc_intF ( mengerc_data *md, double *func );
boolean _mengerc_gradIntF ( mengerc_data *md, double *func, double *grad );
boolean _mengerc_hessIntF ( mengerc_data *md, double *func, double *grad,
                            double *hess );

/* ///////////////////////////////////////////////////////////////////////// */
boolean mengerc_intD ( mengerc_data *md,
                       int lkn, double *knots, point3d *cpoints,
                       double *dl, double *acp );
boolean mengerc_gradIntD ( mengerc_data *md,
                           int lkn, double *knots, point3d *cpoints,
                           double *dl, double *grdl, double *acp, double *gracp );
boolean mengerc_hessIntD ( mengerc_data *md,
                           int lkn, double *knots, point3d *cpoints,
                           double *dl, double *grdl, double *hesdl,
                           double *acp, double *gracp, double *hesacp );

/* ///////////////////////////////////////////////////////////////////////// */
boolean mengerc_IntegralMengerf ( int n, void *usrdata, double *x, double *f );
boolean mengerc_IntegralMengerfg ( int n, void *usrdata, double *x,
                                   double *f, double *g );
boolean mengerc_IntegralMengerfgh ( int n, void *usrdata, double *x,
                                    double *f, double *g, double *h );
boolean mengerc_IntegralMengerTransC ( int n, void *usrdata, double *x );
boolean mengerc_HomotopyTest ( int n, void *usrdata, double *x0, double *x1,
                               boolean *went_out );

/* ///////////////////////////////////////////////////////////////////////// */
int mengerc_FindRemotestPoint ( int np, point3d *cpoints, point3d *sc );
int mengerc_ModifyRemotestPoint ( int np, point3d *cpoints, point3d *sc, int mdi );

boolean mengerc_OptPenaltyParams1 ( mengerc_data *md, boolean wide );   
boolean mengerc_OptPenaltyParams2 ( mengerc_data *md );
boolean mengerc_OptPenaltyParams3 ( mengerc_data *md );

/* ///////////////////////////////////////////////////////////////////////// */
boolean mengerc_InitMCOptimization ( int deg, int lkn, double *knots,
                                     point3d *cpoints, double w,
                                     double penalty_param[MENGERC_NPPARAM],
                                     int nqkn, int npthr, int opt,
                                     mengerc_data *md );
boolean mengerc_IterMCOptimization ( mengerc_data *md, boolean *finished );

boolean mengerc_OptimizeMengerCurvature (
                      int deg, int lkn, double *knots, point3d *cpoints,
                      double w, double penalty_param[MENGERC_NPPARAM],
                      int nqkn, int npthr, int opt, int maxit,
                      void (*outiter)(void *usrdata,
                                      boolean ppopt, int mdi,
                                      int it, int itres, double f, double g),
                      void *usrdata );

#ifdef __cplusplus
}
#endif

#endif

