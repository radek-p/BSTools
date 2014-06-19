
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2014                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */ 

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "raybez.h"
#include "mengerc.h"

#include "mengercprivate.h"

/* ///////////////////////////////////////////////////////////////////////// */ 
void mengerc_GravityCentre ( int ncp, point3d *cpoints, point3d *sc )
{
  int i;

  *sc = cpoints[0];
  for ( i = 1; i < ncp; i++ )
    AddVector3d ( sc, &cpoints[i], sc );
  MultVector3d ( 1.0/(double)ncp, sc, sc );
} /*mengerc_GravityCentre*/

int mengerc_FindRemotestPoint ( int np, point3d *cpoints, point3d *sc )
{
  int      i, j;
  double   d, e;
  vector3d vi;

  j = 0;
  e = -1.0;
  for ( i = 0; i < np; i++ ) {
    SubtractPoints3d ( &cpoints[i], sc, &vi );
    d = DotProduct3d ( &vi, &vi );
    if ( d > e ) {
      e = d;
      j = i;
    }
  }
  return j;
} /*mengerc_FindRemotestPoint*/

int mengerc_ModifyRemotestPoint ( int np, point3d *cpoints, point3d *sc, int mdi )
{
  int      j, k;
  double   d, e;
  vector3d v;

  j = (mdi + 1) % np;
  k = (mdi + np - 1) % np;
  SubtractPoints3d ( &cpoints[mdi], sc, &v );
  d = DotProduct3d ( &v, &v );
  SubtractPoints3d ( &cpoints[j], sc, &v );
  e = DotProduct3d ( &v, &v );
  if ( e > d )
    return mengerc_FindRemotestPoint ( np, cpoints, sc );
  SubtractPoints3d ( &cpoints[k], sc, &v );
  e = DotProduct3d ( &v, &v );
  if ( e > d )
    return mengerc_FindRemotestPoint ( np, cpoints, sc );
  return mdi;
} /*mengerc_ModifyRemotestPoint*/

/* ///////////////////////////////////////////////////////////////////////// */ 
boolean mengerc_IntegralMengerf ( int n, void *usrdata, double *x, double *f )
{
  void     *sp;
  mengerc_data *md;
  int      deg, lkn, clcK, mdi;
  double   *knots, *penalty_param, clcT2;
  point3d  *cpoints, sc;
  vector3d vv;
  double   kM, kkM, L, L0, L02, R1, R2, R3, R4, R5, a, q;
  double   r1, r2, r3;

  sp = pkv_GetScratchMemTop ();
  md = (mengerc_data*)usrdata;
  deg  = md->n;
  lkn = md->lkn;
  knots = md->knots;
  cpoints = md->cpoints;
  clcK = lkn - 2*deg;  /* n == 3*clcK */
  clcT2 = (double)(clcK*clcK);
  if ( &cpoints[0].x != x )
    memcpy ( cpoints, x, n*sizeof(double) );
  memcpy ( &cpoints[clcK], x, deg*sizeof(point3d) );
  penalty_param = md->penalty_param;

  if ( !memcmp ( x, md->fx, n*sizeof(double) ) ) {
    kkM = md->ffkM;  R1 = md->ffR1;  R2 = md->ffR2;
    R3 = md->ffR3;   R4 = md->ffR4;  R5 = md->ffR5;
    goto sum_up;
  }
  if ( !memcmp ( x, md->gx, n*sizeof(double) ) ) {
    kkM = md->gfkM;  R1 = md->gfR1;  R2 = md->gfR2;
    R3 = md->gfR3;   R4 = md->gfR4;  R5 = md->gfR5;
    goto sum_up;
  }
  if ( !memcmp ( x, md->hx, n*sizeof(double) ) ) {
    kkM = md->hfkM;  R1 = md->hfR1;  R2 = md->hfR2;
    R3 = md->hfR3;   R4 = md->hfR4;  R5 = md->hfR5;
    goto sum_up;
  }

  L0 = md->L;
  L02 = L0*L0;
  if ( md->alt_scale )
    q = 1.0/(md->w - 3.0);
  else
    q = md->w - 3.0;

        /* obliczanie krzywizny calkowej Mengera */
  if ( !_mengerc_intF ( md, &kM ) )
    goto failure;
        /* obliczanie dlugosci krzywej i kary za zmienna predkosc parametryzacji */
  if ( !mengerc_intD ( md, lkn, knots, cpoints, &L, &R4 ) )
    goto failure;
  R4 /= L0;
        /* kompensowanie dlugosci */
  if ( md->alt_scale )
    kkM = L*pow ( kM, q );
  else
    kkM = kM;
        /* obliczanie kary za niewlasciwa dlugosc */
  a = L/L0 - 1.0;
  R1 = a*a;
        /* obliczanie kary za zle polozony srodek ciezkosci */
  mengerc_GravityCentre ( clcK, cpoints, &sc );
  md->sc = sc;
  R2 = DotProduct3d ( &sc, &sc )/L02;
        /* obliczanie kary za niewlasciwe zorientowanie */
  r1 = cpoints[1].y-sc.y;
  r2 = cpoints[1].z-sc.z;
  r3 = cpoints[2].z-cpoints[0].z;
  R3 = (r1*r1 + r2*r2 + 0.25*clcT2*r3*r3)/L02;
        /* obliczanie kary za przesuniecie parametryzacji */
  if ( (mdi = md->mdi) > 0 )
    SubtractPoints3d ( &cpoints[mdi+1], &cpoints[mdi-1], &vv );
  else
    SubtractPoints3d ( &cpoints[1], &cpoints[lkn-2*deg-1], &vv );
  R5 = DotProduct3d ( &cpoints[mdi], &vv );
  R5 *= clcT2*R5/(L02*L02);
        /* zapamietywanie skladnikow */
  memcpy ( md->fx, x, n*sizeof(double) );
  md->ffkM = kkM;  md->ffR1 = R1;  md->ffR2 = R2;
  md->ffR3 = R3;   md->ffR4 = R4;  md->ffR5 = R5;

sum_up:
  *f = kkM + penalty_param[0]*R1 + penalty_param[1]*R2 +
       penalty_param[2]*R3 + penalty_param[3]*R4 + penalty_param[4]*R5;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  memset ( md->fx, 0, n*sizeof(double) );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mengerc_IntegralMengerf*/

boolean mengerc_IntegralMengerfg ( int n, void *usrdata, double *x,
                                   double *f, double *g )
{
  void *sp;
  mengerc_data *md;
  int      deg, lkn, clcK;
  double   *knots, *penalty_param, clcT2;
  point3d  *cpoints, sc, sc1;
  double   kM, kkM, L, L0, L02, R1, R2, R3, R4, R5, a, q;
  double   KtoQ, LtoR;
  double   r1, r2, r3, r4;
  double   *grkM, *grL, *grR1, *grR2, *grR3, *grR4, *grR5;
  int      i, mdi, m3;
  vector3d vv;

  sp = pkv_GetScratchMemTop ();
  md = (mengerc_data*)usrdata;

  deg  = md->n;
  lkn = md->lkn;
  knots = md->knots;
  cpoints = md->cpoints;
  clcK = lkn - 2*deg;  /* n == 3*clcK */
  clcT2 = (double)(clcK*clcK);
  memcpy ( cpoints, x, n*sizeof(double) );
  memcpy ( &cpoints[clcK], x, deg*sizeof(point3d) );
  penalty_param = md->penalty_param;
  grL = pkv_GetScratchMemd ( n );
  if ( !grL )
    goto failure;

  if ( !memcmp ( x, md->hx, n*sizeof(double) ) ) {
    kkM = md->hfkM;  R1 = md->hfR1;  R2 = md->hfR2;
    R3 = md->hfR3;   R4 = md->hfR4;  R5 = md->hfR5;
    grkM = md->hgkM;  grR1 = md->hgR1;  grR2 = md->hgR2;
    grR3 = md->hgR3;  grR4 = md->hgR4;  grR5 = md->hgR5;
    goto sum_up;
  }
  grkM = md->ggkM;  grR1 = md->ggR1;  grR2 = md->ggR2;
  grR3 = md->ggR3;  grR4 = md->ggR4;  grR5 = md->ggR5;
  if ( !memcmp ( x, md->gx, n*sizeof(double) ) ) {
    kkM = md->gfkM;  R1 = md->gfR1;  R2 = md->gfR2;
    R3 = md->gfR3;   R4 = md->gfR4;  R5 = md->gfR5;
    goto sum_up;
  }

  L0 = md->L;
  L02 = L0*L0;
  if ( md->alt_scale )
    q = 1.0/(md->w - 3.0);
  else
    q = md->w - 3.0;

        /* obliczanie krzywizny calkowej Mengera */
  if ( !_mengerc_gradIntF ( md, &kM, grkM ) )
    goto failure;
        /* obliczanie dlugosci krzywej i kary za zmienna predkosc parametryzacji */
  if ( !mengerc_gradIntD ( md, lkn, knots, cpoints, &L, grL, &R4, grR4 ) )
    goto failure;
  R4 /= L0;
  for ( i = 0; i < n; i++ )
    grR4[i] /= L0;
        /* kompensowanie dlugosci */
  if ( md->alt_scale ) {
    KtoQ = pow ( kM, q );
    kkM = L*KtoQ;
    for ( i = 0; i < n; i++ )
      grkM[i] = (q*L*grkM[i]/kM+grL[i])*KtoQ;
  }
  else {
    LtoR = pow ( L, q );
    kkM = kM*LtoR;
    for ( i = 0; i < n; i++ )
      grkM[i] = (q*kM*grL[i]/L+grkM[i])*LtoR;
  }
        /* obliczanie kary za niewlasciwa dlugosc */
  a = L/L0 - 1.0;
  R1 = a*a;
  a = (a+a)/L0;
  for ( i=0 ; i<n; i++ )
    grR1[i] = a*grL[i];
        /* obliczanie kary za zle polozony srodek ciezkosci */
  mengerc_GravityCentre ( clcK, cpoints, &sc );
  md->sc = sc;
  R2 = DotProduct3d ( &sc, &sc )/L02;
  MultVector3d ( 6.0/((double)n*L02), &sc, &sc1 );
  for ( i = 0; i < n; i += 3 ) {
    grR2[i] = sc1.x;
    grR2[i+1] = sc1.y;
    grR2[i+2] = sc1.z;
  }
        /* obliczanie kary za niewlasciwe zorientowanie */
  r1 = cpoints[1].y-sc.y;
  r2 = cpoints[1].z-sc.z;
  r3 = cpoints[2].z-cpoints[0].z;
  R3 = (r1*r1 + r2*r2 + 0.25*clcT2*r3*r3)/L02;
  memset ( grR3, 0, n*sizeof(double) );
  r4 = 2.0/((double)clcK*L02);
  r1 *= r4;
  r2 *= r4;
  grR3[1] = -r1;
  grR3[2] = -r2;
  grR3[4] = (double)(clcK-1)*r1;
  grR3[5] = (double)(clcK-1)*r2;
  for ( i = 7; i < n-1; i += 3 ) {
    grR3[i] = -r1;
    grR3[i+1] = -r2;
  }
  grR3[8] += r3 = 0.5*clcT2*r3/L02;
  grR3[2] -= r3;
        /* obliczanie kary za przesuniecie parametryzacji */
  if ( (mdi = md->mdi) > 0 )
    SubtractPoints3d ( &cpoints[mdi+1], &cpoints[mdi-1], &vv );
  else
    SubtractPoints3d ( &cpoints[1], &cpoints[lkn-2*deg-1], &vv );
  memset ( grR5, 0, n*sizeof(double) );
  a = DotProduct3d ( &cpoints[mdi], &vv );
  R5 = clcT2*a*a/(L02*L02);
  a = clcT2*(a+a)/(L02*L02);
  m3 = 3*mdi;
  if ( mdi == 0 ) {
    grR5[n-3] = -(grR5[3] = a*cpoints[0].x);
    grR5[n-2] = -(grR5[4] = a*cpoints[0].y);
    grR5[n-1] = -(grR5[5] = a*cpoints[0].z);
  }
  else if ( mdi == clcK-1 ) {
    grR5[m3-3] = -(grR5[0] = a*cpoints[mdi].x);
    grR5[m3-2] = -(grR5[1] = a*cpoints[mdi].y);
    grR5[m3-1] = -(grR5[2] = a*cpoints[mdi].z);
  }
  else {
    grR5[m3-3] = -(grR5[m3+3] = a*cpoints[mdi].x);
    grR5[m3-2] = -(grR5[m3+4] = a*cpoints[mdi].y);
    grR5[m3-1] = -(grR5[m3+5] = a*cpoints[mdi].z);
  }
  grR5[m3]   = a*vv.x;
  grR5[m3+1] = a*vv.y;
  grR5[m3+2] = a*vv.z;

        /* zapamietywanie skladnikow */
  memcpy ( md->gx, x, n*sizeof(double) );
  md->gfkM = kkM;  md->gfR1 = R1;  md->gfR2 = R2;
  md->gfR3 = R3;   md->gfR4 = R4;  md->gfR5 = R5;

sum_up:
  *f = kkM + penalty_param[0]*R1 + penalty_param[1]*R2 +
       penalty_param[2]*R3 + penalty_param[3]*R4 + penalty_param[4]*R5;
  for ( i = 0; i < n; i++ )
    g[i] = grkM[i] + penalty_param[0]*grR1[i] + penalty_param[1]*grR2[i] +
           penalty_param[2]*grR3[i] + penalty_param[3]*grR4[i] +
           penalty_param[4]*grR5[i];
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  memset ( md->gx, 0, n*sizeof(double) );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mengerc_IntegralMengerfg*/

boolean mengerc_IntegralMengerfgh ( int n, void *usrdata, double *x,
                                    double *f, double *g, double *h )
{
  void *sp;
  mengerc_data *md;
  int      deg, lkn, clcK;
  double   *knots, *penalty_param, clcT2, twoT2;
  point3d  *cpoints, sc, sc1;
  vector3d vv;
  double   kM, kkM, L, L0, L02, R1, R2, R3, R4, R5, a, b, q, qm1;
  double   KtoQ, K2, KtoQm2, LtoR, L2, LtoRm2;
  double   r1, r2, r3, r4;
  double   *gkM, *gL, *gR1, *gR2, *gR3, *gR4, *gR5;
  double   *hkM, *hL, *hR1, *hR2, *hR3, *hR4, *hR5;
  int      i, j, nn, mdi, m3;

  sp = pkv_GetScratchMemTop ();
  md = (mengerc_data*)usrdata;

  penalty_param = md->penalty_param;
  nn = (n*(n+1))/2;
  gL = pkv_GetScratchMemd ( n );
  hL  = pkv_GetScratchMemd ( nn );
  if ( !gL || !hL )
    goto failure;
  gkM = md->hgkM;  gR1 = md->hgR1;  gR2 = md->hgR2;
  gR3 = md->hgR3;  gR4 = md->hgR4;  gR5 = md->hgR5;
  hkM = md->hhkM;  hR1 = md->hhR1;  hR2 = md->hhR2;
  hR3 = md->hhR3;  hR4 = md->hhR4;  hR5 = md->hhR5;

  if ( !memcmp ( x, md->hx, n*sizeof(double) ) ) {
    kkM = md->hfkM;  R1 = md->hfR1;  R2 = md->hfR2;
    R3 = md->hfR3;   R4 = md->hfR4;  R5 = md->hfR5;
    goto sum_up;
  }

  deg  = md->n;
  lkn = md->lkn;
  knots = md->knots;
  cpoints = md->cpoints;
  clcK = lkn - 2*deg;  /* n == 3*clcK */
  clcT2 = (double)(clcK*clcK);
  twoT2 = clcT2+clcT2;
  memcpy ( cpoints, x, n*sizeof(double) );
  memcpy ( &cpoints[clcK], x, deg*sizeof(point3d) );
  L0 = md->L;
  L02 = L0*L0;
  if ( md->alt_scale )
    q = 1.0/(md->w - 3.0);
  else
    q = md->w - 3.0;

        /* obliczanie krzywizny calkowej Mengera */
  if ( !_mengerc_hessIntF ( md, &kM, gkM, hkM ) )
    goto failure;
        /* obliczanie dlugosci krzywej i kary za zmienna predkosc parametryzacji */
  if ( !mengerc_hessIntD ( md, lkn, knots, cpoints, &L, gL, hL, &R4, gR4, hR4 ) )
    goto failure;
  R4 /= L0;
  for ( i = 0; i < n; i++ )
    gR4[i] /= L0;
  for ( i = 0; i < nn; i++ )
    hR4[i] /= L0;
        /* kompensowanie dlugosci */
  if ( md->alt_scale ) {
    qm1 = q-1.0;
    KtoQ = pow ( kM, q );
    K2 = kM*kM;
    KtoQm2 = KtoQ/K2;
    kkM = L*KtoQ;
    for ( i = 0; i < n; i++ )
      for ( j = 0; j <= i; j++ )
        hkM[pkn_LowerTrMatIndex(i,j)] = KtoQm2*(q*(qm1*gkM[i]*gkM[j]*L +
                   kM*(hkM[pkn_LowerTrMatIndex(i,j)]*L +
                          gkM[i]*gL[j] + gkM[j]*gL[i])) +
                   K2*hL[pkn_LowerTrMatIndex(i,j)]);
    for ( i = 0; i < n; i++ )
      gkM[i] = (q*L*gkM[i]/kM+gL[i])*KtoQ;
  }
  else {
    qm1 = q-1.0;
    LtoR = pow ( L, q );
    L2 = L*L;
    LtoRm2 = LtoR/L2;
    kkM = kM*LtoR;
    for ( i = 0; i < n; i++ )
      for ( j = 0; j <= i; j++ )
        hkM[pkn_LowerTrMatIndex(i,j)] = LtoRm2*(q*(qm1*gL[i]*gL[j]*kM +
                   L*(hL[pkn_LowerTrMatIndex(i,j)]*kM +
                      gL[i]*gkM[j] + gL[j]*gkM[i])) +
                   L2*hkM[pkn_LowerTrMatIndex(i,j)]);
    for ( i = 0; i < n; i++ )
      gkM[i] = (q*kM*gL[i]/L+gkM[i])*LtoR;
  }
        /* obliczanie kary za niewlasciwa dlugosc */
  a = L/L0 - 1.0;
  R1 = a*a;
  a = (a+a)/L0;
  for ( i = 0; i < n; i++ )
    gR1[i] = a*gL[i];
  for ( i = 0; i < n; i++ )
    for ( j = 0; j <= i; j++ )
      hR1[pkn_LowerTrMatIndex(i,j)] = a*hL[pkn_LowerTrMatIndex(i,j)] +
                                      2.0*gL[i]*gL[j]/L02;
        /* obliczanie kary za zle polozony srodek ciezkosci */
  mengerc_GravityCentre ( clcK, cpoints, &sc );
  md->sc = sc;
  R2 = DotProduct3d ( &sc, &sc )/L02;
  MultVector3d ( 6.0/((double)n*L02), &sc, &sc1 );
  for ( i = 0; i < n; i += 3 ) {
    gR2[i] = sc1.x;
    gR2[i+1] = sc1.y;
    gR2[i+2] = sc1.z;
  }
  a = 18.0/((double)(n*n)*L02);
  memset ( hR2, 0, nn*sizeof(double) );
  for ( i = 0; i < n; i += 3 )
    for ( j = 0; j <= i; j += 3 )
      hR2[pkn_LowerTrMatIndex(i,j)] = hR2[pkn_LowerTrMatIndex(i+1,j+1)] =
      hR2[pkn_LowerTrMatIndex(i+2,j+2)] = a;
        /* obliczanie kary za niewlasciwe zorientowanie */
  r1 = cpoints[1].y-sc.y;
  r2 = cpoints[1].z-sc.z;
  r3 = cpoints[2].z-cpoints[0].z;
  R3 = (r1*r1 + r2*r2 + 0.25*clcT2*r3*r3)/L02;
  memset ( gR3, 0, n*sizeof(double) );
  r4 = 2.0/((double)clcK*L02);
  r1 *= r4;
  r2 *= r4;
  gR3[1] = -r1;
  gR3[2] = -r2;
  gR3[4] = (double)(clcK-1)*r1;
  gR3[5] = (double)(clcK-1)*r2;
  for ( i = 7; i < n-1; i += 3 ) {
    gR3[i] = -r1;
    gR3[i+1] = -r2;
  }
  gR3[8] += r3 = 0.5*clcT2*r3/L02;
  gR3[2] -= r3;
  memset ( hR3, 0, nn*sizeof(double) );
  r3 = 2.0/(clcT2*L02);
  r2 = -r3*(double)(clcK-1);
  r1 = -r2*(double)(clcK-1);
  for ( i = 7; i < n-1; i += 3 )
    for ( j = 1; j <= i; j += 3 )
      hR3[pkn_LowerTrMatIndex(i,j)] = hR3[pkn_LowerTrMatIndex(i+1,j+1)] = r3;
  hR3[pkn_LowerTrMatIndex(1,1)] = hR3[pkn_LowerTrMatIndex(2,2)] = r3;
  hR3[pkn_LowerTrMatIndex(4,4)] = hR3[pkn_LowerTrMatIndex(5,5)] = r1;
  hR3[pkn_LowerTrMatIndex(4,1)] = hR3[pkn_LowerTrMatIndex(5,2)] = r2;
  for ( i = 7; i < n-1; i += 3 )
    hR3[pkn_LowerTrMatIndex(i,4)] = hR3[pkn_LowerTrMatIndex(i+1,5)] = r2;
  r4 = 0.5*clcT2/L02;
  hR3[pkn_LowerTrMatIndex(8,8)] += r4;
  hR3[pkn_LowerTrMatIndex(2,2)] += r4;
  hR3[pkn_LowerTrMatIndex(8,2)] -= r4;
        /* obliczanie kary za przesuniecie parametryzacji */
  if ( (mdi = md->mdi) > 0 )
    SubtractPoints3d ( &cpoints[mdi+1], &cpoints[mdi-1], &vv );
  else
    SubtractPoints3d ( &cpoints[1], &cpoints[lkn-2*deg-1], &vv );
  memset ( gR5, 0, n*sizeof(double) );
  memset ( hR5, 0, nn*sizeof(double) );
  a = DotProduct3d ( &cpoints[mdi], &vv );
  R5 = clcT2*a*a/(L02*L02);
  b = twoT2/(L02*L02);
  a *= b;
  m3 = 3*mdi;
  if ( mdi == 0 ) {
    gR5[n-3] = -(gR5[3] = a*cpoints[0].x);
    gR5[n-2] = -(gR5[4] = a*cpoints[0].y);
    gR5[n-1] = -(gR5[5] = a*cpoints[0].z);
    hR5[pkn_LowerTrMatIndex(m3+3,m3+3)] = hR5[pkn_LowerTrMatIndex(n-3,n-3)] = b*cpoints[mdi].x*cpoints[mdi].x;
    hR5[pkn_LowerTrMatIndex(m3+4,m3+3)] = hR5[pkn_LowerTrMatIndex(n-2,n-3)] = b*cpoints[mdi].y*cpoints[mdi].x;
    hR5[pkn_LowerTrMatIndex(m3+4,m3+4)] = hR5[pkn_LowerTrMatIndex(n-2,n-2)] = b*cpoints[mdi].y*cpoints[mdi].y;
    hR5[pkn_LowerTrMatIndex(m3+5,m3+3)] = hR5[pkn_LowerTrMatIndex(n-1,n-3)] = b*cpoints[mdi].z*cpoints[mdi].x;
    hR5[pkn_LowerTrMatIndex(m3+5,m3+4)] = hR5[pkn_LowerTrMatIndex(n-1,n-2)] = b*cpoints[mdi].z*cpoints[mdi].y;
    hR5[pkn_LowerTrMatIndex(m3+5,m3+5)] = hR5[pkn_LowerTrMatIndex(n-1,n-1)] = b*cpoints[mdi].z*cpoints[mdi].z;
    hR5[pkn_LowerTrMatIndex(n-3,m3+3)] = -hR5[pkn_LowerTrMatIndex(m3+3,m3+3)];
    hR5[pkn_LowerTrMatIndex(n-2,m3+3)] = hR5[pkn_LowerTrMatIndex(n-3,m3+4)] = -hR5[pkn_LowerTrMatIndex(m3+4,m3+3)];
    hR5[pkn_LowerTrMatIndex(n-2,m3+4)] = -hR5[pkn_LowerTrMatIndex(m3+4,m3+4)];
    hR5[pkn_LowerTrMatIndex(n-1,m3+3)] = hR5[pkn_LowerTrMatIndex(n-3,m3+5)] = -hR5[pkn_LowerTrMatIndex(m3+5,m3+3)];
    hR5[pkn_LowerTrMatIndex(n-1,m3+4)] = hR5[pkn_LowerTrMatIndex(n-2,m3+5)] = -hR5[pkn_LowerTrMatIndex(m3+5,m3+4)];
    hR5[pkn_LowerTrMatIndex(n-1,m3+5)] = -hR5[pkn_LowerTrMatIndex(m3+5,m3+5)];

    hR5[pkn_LowerTrMatIndex(n-3,m3)] = -(hR5[pkn_LowerTrMatIndex(m3+3,m3)] = b*vv.x*cpoints[mdi].x + a);
    hR5[pkn_LowerTrMatIndex(n-2,m3)] = -(hR5[pkn_LowerTrMatIndex(m3+4,m3)] = b*vv.x*cpoints[mdi].y);
    hR5[pkn_LowerTrMatIndex(n-1,m3)] = -(hR5[pkn_LowerTrMatIndex(m3+5,m3)] = b*vv.x*cpoints[mdi].z);
    hR5[pkn_LowerTrMatIndex(n-3,m3+1)] = -(hR5[pkn_LowerTrMatIndex(m3+3,m3+1)] = b*vv.y*cpoints[mdi].x);
    hR5[pkn_LowerTrMatIndex(n-2,m3+1)] = -(hR5[pkn_LowerTrMatIndex(m3+4,m3+1)] = b*vv.y*cpoints[mdi].y + a);
    hR5[pkn_LowerTrMatIndex(n-1,m3+1)] = -(hR5[pkn_LowerTrMatIndex(m3+5,m3+1)] = b*vv.y*cpoints[mdi].z);
    hR5[pkn_LowerTrMatIndex(n-3,m3+2)] = -(hR5[pkn_LowerTrMatIndex(m3+3,m3+2)] = b*vv.z*cpoints[mdi].x);
    hR5[pkn_LowerTrMatIndex(n-2,m3+2)] = -(hR5[pkn_LowerTrMatIndex(m3+4,m3+2)] = b*vv.z*cpoints[mdi].y);
    hR5[pkn_LowerTrMatIndex(n-1,m3+2)] = -(hR5[pkn_LowerTrMatIndex(m3+5,m3+2)] = b*vv.z*cpoints[mdi].z + a);

  }
  else if ( mdi == clcK-1 ) {
    gR5[m3-3] = -(gR5[0] = a*cpoints[mdi].x);
    gR5[m3-2] = -(gR5[1] = a*cpoints[mdi].y);
    gR5[m3-1] = -(gR5[2] = a*cpoints[mdi].z);
    hR5[pkn_LowerTrMatIndex(0,0)] = hR5[pkn_LowerTrMatIndex(m3-3,m3-3)] = b*cpoints[mdi].x*cpoints[mdi].x;
    hR5[pkn_LowerTrMatIndex(1,0)] = hR5[pkn_LowerTrMatIndex(m3-2,m3-3)] = b*cpoints[mdi].y*cpoints[mdi].x;
    hR5[pkn_LowerTrMatIndex(1,1)] = hR5[pkn_LowerTrMatIndex(m3-2,m3-2)] = b*cpoints[mdi].y*cpoints[mdi].y;
    hR5[pkn_LowerTrMatIndex(2,0)] = hR5[pkn_LowerTrMatIndex(m3-1,m3-3)] = b*cpoints[mdi].z*cpoints[mdi].x;
    hR5[pkn_LowerTrMatIndex(2,1)] = hR5[pkn_LowerTrMatIndex(m3-1,m3-2)] = b*cpoints[mdi].z*cpoints[mdi].y;
    hR5[pkn_LowerTrMatIndex(2,2)] = hR5[pkn_LowerTrMatIndex(m3-1,m3-1)] = b*cpoints[mdi].z*cpoints[mdi].z;
    hR5[pkn_LowerTrMatIndex(m3-3,0)] = -hR5[pkn_LowerTrMatIndex(0,0)];
    hR5[pkn_LowerTrMatIndex(m3-2,0)] = hR5[pkn_LowerTrMatIndex(m3-3,1)] = -hR5[pkn_LowerTrMatIndex(1,0)];
    hR5[pkn_LowerTrMatIndex(m3-2,1)] = -hR5[pkn_LowerTrMatIndex(1,1)];
    hR5[pkn_LowerTrMatIndex(m3-1,0)] = hR5[pkn_LowerTrMatIndex(m3-3,2)] = -hR5[pkn_LowerTrMatIndex(2,0)];
    hR5[pkn_LowerTrMatIndex(m3-1,1)] = hR5[pkn_LowerTrMatIndex(m3-2,2)] = -hR5[pkn_LowerTrMatIndex(2,1)];
    hR5[pkn_LowerTrMatIndex(m3-1,2)] = -hR5[pkn_LowerTrMatIndex(2,2)];
    hR5[pkn_LowerTrMatIndex(m3,m3-3)] = -(hR5[pkn_LowerTrMatIndex(m3,0)] = b*vv.x*cpoints[mdi].x + a);
    hR5[pkn_LowerTrMatIndex(m3,m3-2)] = -(hR5[pkn_LowerTrMatIndex(m3,1)] = b*vv.x*cpoints[mdi].y);
    hR5[pkn_LowerTrMatIndex(m3,m3-1)] = -(hR5[pkn_LowerTrMatIndex(m3,2)] = b*vv.x*cpoints[mdi].z);
    hR5[pkn_LowerTrMatIndex(m3+1,m3-3)] = -(hR5[pkn_LowerTrMatIndex(m3+1,0)] = b*vv.y*cpoints[mdi].x);
    hR5[pkn_LowerTrMatIndex(m3+1,m3-2)] = -(hR5[pkn_LowerTrMatIndex(m3+1,1)] = b*vv.y*cpoints[mdi].y + a);
    hR5[pkn_LowerTrMatIndex(m3+1,m3-1)] = -(hR5[pkn_LowerTrMatIndex(m3+1,2)] = b*vv.y*cpoints[mdi].z);
    hR5[pkn_LowerTrMatIndex(m3+2,m3-3)] = -(hR5[pkn_LowerTrMatIndex(m3+2,0)] = b*vv.z*cpoints[mdi].x);
    hR5[pkn_LowerTrMatIndex(m3+2,m3-2)] = -(hR5[pkn_LowerTrMatIndex(m3+2,1)] = b*vv.z*cpoints[mdi].y);
    hR5[pkn_LowerTrMatIndex(m3+2,m3-1)] = -(hR5[pkn_LowerTrMatIndex(m3+2,2)] = b*vv.z*cpoints[mdi].z + a);
  }
  else {
    gR5[m3-3] = -(gR5[m3+3] = a*cpoints[mdi].x);
    gR5[m3-2] = -(gR5[m3+4] = a*cpoints[mdi].y);
    gR5[m3-1] = -(gR5[m3+5] = a*cpoints[mdi].z);
    hR5[pkn_LowerTrMatIndex(m3+3,m3+3)] = hR5[pkn_LowerTrMatIndex(m3-3,m3-3)] = b*cpoints[mdi].x*cpoints[mdi].x;
    hR5[pkn_LowerTrMatIndex(m3+4,m3+3)] = hR5[pkn_LowerTrMatIndex(m3-2,m3-3)] = b*cpoints[mdi].y*cpoints[mdi].x;
    hR5[pkn_LowerTrMatIndex(m3+4,m3+4)] = hR5[pkn_LowerTrMatIndex(m3-2,m3-2)] = b*cpoints[mdi].y*cpoints[mdi].y;
    hR5[pkn_LowerTrMatIndex(m3+5,m3+3)] = hR5[pkn_LowerTrMatIndex(m3-1,m3-3)] = b*cpoints[mdi].z*cpoints[mdi].x;
    hR5[pkn_LowerTrMatIndex(m3+5,m3+4)] = hR5[pkn_LowerTrMatIndex(m3-1,m3-2)] = b*cpoints[mdi].z*cpoints[mdi].y;
    hR5[pkn_LowerTrMatIndex(m3+5,m3+5)] = hR5[pkn_LowerTrMatIndex(m3-1,m3-1)] = b*cpoints[mdi].z*cpoints[mdi].z;
    hR5[pkn_LowerTrMatIndex(m3+3,m3-3)] = -hR5[pkn_LowerTrMatIndex(m3+3,m3+3)];
    hR5[pkn_LowerTrMatIndex(m3+3,m3-2)] = hR5[pkn_LowerTrMatIndex(m3+4,m3-3)] = -hR5[pkn_LowerTrMatIndex(m3+4,m3+3)];
    hR5[pkn_LowerTrMatIndex(m3+4,m3-2)] = -hR5[pkn_LowerTrMatIndex(m3+4,m3+4)];
    hR5[pkn_LowerTrMatIndex(m3+3,m3-1)] = hR5[pkn_LowerTrMatIndex(m3+5,m3-3)] = -hR5[pkn_LowerTrMatIndex(m3+5,m3+3)];
    hR5[pkn_LowerTrMatIndex(m3+4,m3-1)] = hR5[pkn_LowerTrMatIndex(m3+5,m3-2)] = -hR5[pkn_LowerTrMatIndex(m3+5,m3+4)];
    hR5[pkn_LowerTrMatIndex(m3+5,m3-1)] = -hR5[pkn_LowerTrMatIndex(m3+5,m3+5)];
    hR5[pkn_LowerTrMatIndex(m3,m3-3)] = -(hR5[pkn_LowerTrMatIndex(m3+3,m3)] = b*vv.x*cpoints[mdi].x + a);
    hR5[pkn_LowerTrMatIndex(m3,m3-2)] = -(hR5[pkn_LowerTrMatIndex(m3+4,m3)] = b*vv.x*cpoints[mdi].y);
    hR5[pkn_LowerTrMatIndex(m3,m3-1)] = -(hR5[pkn_LowerTrMatIndex(m3+5,m3)] = b*vv.x*cpoints[mdi].z);
    hR5[pkn_LowerTrMatIndex(m3+1,m3-3)] = -(hR5[pkn_LowerTrMatIndex(m3+3,m3+1)] = b*vv.y*cpoints[mdi].x);
    hR5[pkn_LowerTrMatIndex(m3+1,m3-2)] = -(hR5[pkn_LowerTrMatIndex(m3+4,m3+1)] = b*vv.y*cpoints[mdi].y + a);
    hR5[pkn_LowerTrMatIndex(m3+1,m3-1)] = -(hR5[pkn_LowerTrMatIndex(m3+5,m3+1)] = b*vv.y*cpoints[mdi].z);
    hR5[pkn_LowerTrMatIndex(m3+2,m3-3)] = -(hR5[pkn_LowerTrMatIndex(m3+3,m3+2)] = b*vv.z*cpoints[mdi].x);
    hR5[pkn_LowerTrMatIndex(m3+2,m3-2)] = -(hR5[pkn_LowerTrMatIndex(m3+4,m3+2)] = b*vv.z*cpoints[mdi].y);
    hR5[pkn_LowerTrMatIndex(m3+2,m3-1)] = -(hR5[pkn_LowerTrMatIndex(m3+5,m3+2)] = b*vv.z*cpoints[mdi].z + a);
  }
  gR5[m3]   = a*vv.x;
  gR5[m3+1] = a*vv.y;
  gR5[m3+2] = a*vv.z;
  hR5[pkn_LowerTrMatIndex(m3,m3)]     = b*vv.x*vv.x;
  hR5[pkn_LowerTrMatIndex(m3+1,m3)]   = b*vv.y*vv.x;
  hR5[pkn_LowerTrMatIndex(m3+1,m3+1)] = b*vv.y*vv.y;
  hR5[pkn_LowerTrMatIndex(m3+2,m3)]   = b*vv.z*vv.x;
  hR5[pkn_LowerTrMatIndex(m3+2,m3+1)] = b*vv.z*vv.y;
  hR5[pkn_LowerTrMatIndex(m3+2,m3+2)] = b*vv.z*vv.z;

        /* zapamietywanie skladnikow */
  memcpy ( md->hx, x, n*sizeof(double) );
  md->hfkM = kkM;  md->hfR1 = R1;  md->hfR2 = R2;
  md->hfR3 = R3;   md->hfR4 = R4;  md->hfR5 = R5;

sum_up:
  *f = kkM + penalty_param[0]*R1 + penalty_param[1]*R2 +
       penalty_param[2]*R3 + penalty_param[3]*R4 + penalty_param[4]*R5;
  for ( i = 0; i < n; i++ )
    g[i] = gkM[i] + penalty_param[0]*gR1[i] + penalty_param[1]*gR2[i] +
           penalty_param[2]*gR3[i] + penalty_param[3]*gR4[i] +
           penalty_param[4]*gR5[i];
  for ( i = 0; i < nn; i++ )
    h[i] = hkM[i] + penalty_param[0]*hR1[i] + penalty_param[1]*hR2[i] +
           penalty_param[2]*hR3[i] + penalty_param[3]*hR4[i] +
           penalty_param[4]*hR5[i];
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  memset ( md->hx, 0, n*sizeof(double) );
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mengerc_IntegralMengerfgh*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean mengerc_IntegralMengerTransC ( int n, void *usrdata, double *x )
{
  void     *sp;
  mengerc_data *md;
  int      deg, lkn, clcK, ncp, i;
  point3d  *cpoints, sc;
  double   *knots, L, s, acp;
  vector3d v1, v2, v3;
  trans3d  tr;

  sp = pkv_GetScratchMemTop ();
  md = (mengerc_data*)usrdata;

  deg  = md->n;
  lkn = md->lkn;
  ncp = lkn-deg;
  knots = md->knots;
  cpoints = md->cpoints;
  clcK = lkn - 2*deg;
  memcpy ( cpoints, x, n*sizeof(double) );
  memcpy ( &cpoints[clcK], x, deg*sizeof(point3d) );
        /* skalowanie w celu otrzymania wlasciwej dlugosci */
  if ( !mengerc_intD ( md, lkn, knots, cpoints, &L, &acp ) )
    goto failure;
  s = md->L/L;
  for ( i = 0; i < ncp; i++ )
    MultVector3d ( s, &cpoints[i], &cpoints[i] );
        /* przesuwanie w celu otrzymania srodka ciezkosci w poczatku ukladu */
  mengerc_GravityCentre ( clcK, cpoints, &sc );
  for ( i = 0; i < ncp; i++ )
    SubtractPoints3d ( &cpoints[i], &sc, &cpoints[i] );
  if ( !md->pretransf )
    goto way_out;
        /* obracanie do zadanej pozycji */
  memcpy ( &v1, &cpoints[1], sizeof(vector3d) );
  NormalizeVector3d ( &v1 );
  SubtractPoints3d ( &cpoints[2], &cpoints[0], &v2 );
  OrtVector3d ( &v1, &v2, &v2 );
  NormalizeVector3d ( &v2 );
  CrossProduct3d ( &v1, &v2, &v3 );
  IdentTrans3d ( &tr );
  GeneralAffineTrans3d ( &tr, &v1, &v2, &v3 );
  for ( i = 0; i < ncp; i++ )
    TransContra3d ( &tr, &cpoints[i], &cpoints[i] );

  memcpy ( x, cpoints, n*sizeof(double) );
way_out:
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mengerc_IntegralMengerTransC*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean mengerc_HomotopyTest ( int n, void *usrdata, double *x0, double *x1,
                               boolean *went_out )
{
  void     *sp;
  mengerc_data *md;
  double   tfh[3];
  point3d  *cp0, *cp1;
  int      deg, lkn, ncp, clcK;
  boolean  homot, error;

  sp = pkv_GetScratchMemTop ();
  md = (mengerc_data*)usrdata;
  deg = md->n;
  lkn = md->lkn;
  ncp = lkn - deg;
  clcK = ncp - deg;
  cp0 = md->cpoints;
  cp1 = pkv_GetScratchMem ( ncp*sizeof(point3d) );
  if ( !cp1 )
    goto failure;
  memcpy ( cp0, x0, n*sizeof(double) );
  memcpy ( &cp0[clcK], x0, deg*sizeof(point3d) );
  memcpy ( cp1, x1, n*sizeof(double) );
  memcpy ( &cp1[clcK], x1, deg*sizeof(point3d) );
  homot = rbez_HomotopicClosedBSC3d ( deg, lkn, md->knots, cp0, cp1, tfh, &error );
  if ( error )
    goto failure;
  *went_out = !homot;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*mengerc_HomotopyTest*/

