
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2011                            */
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
#include "bsmesh.h"
#include "egholed.h"
#include "g2blendingd.h"
#include "g2mblendingd.h"

#include "g2blprivated.h"
#include "g2mblprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
void _g2mbl_UCompRDerd ( int nkn2, double *Nitab, int knot,
                         int *cpind, point3d *mvcp, vector3d pder[11] )
{
  int      fi;
  vector3d *_cp;
  double   *Ni;

  memset ( pder, 0, 9*sizeof(vector3d) );
  for ( fi = 0; fi < 16; fi++ ) {
    _cp = &mvcp[cpind[fi]];
    Ni = &Nitab[(fi*nkn2+knot)*9];
    pkn_MultMatrixAddd ( 9, 1, 1, Ni, 3, 0, &_cp->x, 3, &pder[0].x );
  }
  AddVector3d ( &puuu, &puvv, &pder[9] );
  AddVector3d ( &puuv, &pvvv, &pder[10] );
} /*_g2mbl_UCompRDerd*/

void g2mbl_UFuncRSQd ( int nkn, double *qcoeff, double *Nitab,
                       int *cpind, point3d *mvcp,
                       double C, double *ftab )
{
  int      nkn2, i, j, knot;
  double   sum0, sum1, f, g;
  vector3d pder[11];

  nkn2 = nkn*nkn;
  sum0 = 0.0;
  for ( i = knot = 0;  i < nkn;  i++ ) {
    sum1 = 0.0;
    for ( j = 0;  j < nkn;  j++, knot++ ) {
        /* compute the parameterization derivatives */
      _g2mbl_UCompRDerd ( nkn2, Nitab, knot, cpind, mvcp, pder );
        /* compute the integrand terms */
      _g2bl_UFuncSQIntegrandd ( pder, &f, &g );
        /* add the quadrature term */
      if ( C > 0.0 )
        sum1 += (f + C*g)*qcoeff[j];
      else
        sum1 += f*qcoeff[j];
    }
    sum0 += sum1*qcoeff[i];
  }
  *ftab = sum0;
} /*g2mbl_UFuncRSQd*/

void _g2mbl_UCompSDerd ( int nkn2, double *Nitab, int knot, int k, int l,
                         int *cpind, point3d *mvcp, vector3d pder[11] )
{
  int     m, n, p, q, r, s;
  double *Ni, *d;

  r = nkn2*9;
  d = &pder[0].x;
  Ni = &Nitab[knot*9];
  pkn_MultMatrixd ( 9, 1, 1, Ni, 3, 3, &mvcp[cpind[0]].x, 3, d );
  for ( m = 0; m < k; m++ ) {
    n = (l-m+k) % k;
    for ( s = 0, p = 6*n+1, q = (6*m+1)*r;  s < 6;  s++, q += r )
      pkn_MultMatrixAddd ( 9, 1, 1, &Ni[q], 3, 3, &mvcp[cpind[p+s]].x, 3, d );
  }
  AddVector3d ( &puuu, &puvv, &pder[9] );
  AddVector3d ( &puuv, &pvvv, &pder[10] );
} /*_g2mbl_UCompSDerd*/

void g2mbl_UFuncSSQd ( int nkn, double *qcoeff, int k, double *Nitab, double *Jac,
                       int *cpind, point3d *mvcp, double C,
                       double *ftab )
{
  int      nkn2, i, j, l, knot;
  double   sum0, sum1, f, g;
  vector3d pder[11];

  nkn2 = nkn*nkn;
  sum0 = 0.0;
  for ( l = 0; l < k; l++ ) {
    for ( i = knot = 0;  i < nkn;  i++ ) {
      sum1 = 0.0;
      for ( j = 0;  j < nkn;  j++, knot++ ) {
          /* compute the parameterization derivatives */
        _g2mbl_UCompSDerd ( nkn2, Nitab, knot, k, l, cpind, mvcp, pder );
          /* compute the integrand terms */
        _g2bl_UFuncSQIntegrandd ( pder, &f, &g );
          /* add the quadrature term */
        if ( C > 0.0 )
          sum1 += (f + C*g)*Jac[knot]*qcoeff[j];
        else
          sum1 += f*Jac[knot]*qcoeff[j];
      }
      sum0 += sum1*qcoeff[i];
    }
  }
  *ftab = sum0;
} /*g2mbl_UFuncSSQd*/

/* ///////////////////////////////////////////////////////////////////////// */
void _g2mbl_UCompDGStard ( int nkn2, int nf, double *Nitab, int knot, int k, int l,
                           const vector3d *pder, int *cpind, int *nncpi,
                           double *DGstar )
{
  int       fi, p, q;
  vector3d *Dgstar;
  double   *Ni;

  for ( fi = 0; fi < nf; fi++ )
    if ( nncpi[cpind[fi]] >= 0 ) {
      Dgstar = (vector3d*)&DGstar[fi*3*9];
      if ( k == 4 || fi == 0 )
        Ni = &Nitab[(fi*nkn2+knot)*9];
      else {
        p = (fi-1) / 6;
        p = (l - p + k) % k;
        q = (fi-1) % 6;
        Ni = &Nitab[((6*p+q+1)*nkn2+knot)*9];
      }
      MultVector3d ( 2.0*Ni10, &pu, &dg11 );
      MultVector3d ( Ni10, &pv, &dg12 );
      AddVector3Md ( &dg12, &pu, Ni01, &dg12 );
      MultVector3d ( 2.0*Ni01, &pv, &dg22 );

      MultVector3d ( 2.0*Ni20, &pu, &dg11u );
      AddVector3Md ( &dg11u, &puu, 2.0*Ni10, &dg11u );
      MultVector3d ( Ni20, &pv, &dg12u );
      AddVector3Md ( &dg12u, &puu, Ni01, &dg12u );
      AddVector3Md ( &dg12u, &puv, Ni10, &dg12u );
      AddVector3Md ( &dg12u, &pu, Ni11, &dg12u );
      MultVector3d ( 2.0*Ni11, &pv, &dg22u );
      AddVector3Md ( &dg22u, &puv, 2.0*Ni01, &dg22u );

      MultVector3d ( 2.0*Ni11, &pu, &dg11v );
      AddVector3Md ( &dg11v, &puv, 2.0*Ni10, &dg11v );
      MultVector3d ( Ni11, &pv, &dg12v );
      AddVector3Md ( &dg12v, &puv, Ni01, &dg12v );
      AddVector3Md ( &dg12v, &pvv, Ni10, &dg12v );
      AddVector3Md ( &dg12v, &pu, Ni02, &dg12v );
      MultVector3d ( 2.0*Ni02, &pv, &dg22v );
      AddVector3Md ( &dg22v, &pvv, 2.0*Ni01, &dg22v );
    }
} /*_g2mbl_UCompDGStard*/

void _g2mbl_UCompDBStard ( int nkn2, int nf, double *Nitab, int knot, int k, int l,
                           const vector3d *pder, int *cpind, int *nncpi,
                           double *DBstar )
{
  int       fi, p, q;
  vector3d *Dbstar, v, w, z, y;
  double   *Ni;

  for ( fi = 0; fi < nf; fi++ )
    if ( nncpi[cpind[fi]] >= 0 ) {
      Dbstar = (vector3d*)&DBstar[fi*3*11];
      if ( k == 4 || fi == 0 )
        Ni = &Nitab[(fi*nkn2+knot)*9];
      else {
        p = (fi-1) / 6;
        p = (l-p+k) % k;
        q = (fi-1) % 6;
        Ni = &Nitab[((6*p+q+1)*nkn2+knot)*9];
      }

      CrossProduct3d ( &pu, &pv, &v );
      MultVector3d ( Ni20, &v, &dtb11 );
      MultVector3d ( Ni11, &v, &dtb12 );
      MultVector3d ( Ni02, &v, &dtb22 );
      MultVector3d ( Ni30, &v, &dtb11u );
      MultVector3d ( Ni21, &v, &dtb12u );
      MultVector3d ( Ni12, &v, &dtb22u );
      MultVector3d ( Ni03, &v, &dtb22v );
      AddVector3d ( &dtb11u, &dtb22u, &Dbstar[9] );
      AddVector3d ( &dtb12u, &dtb22v, &Dbstar[10] );

      CrossProduct3d ( &puu, &pu, &v );
      AddVector3Md ( &dtb11, &v, Ni01, &dtb11 );
      AddVector3Md ( &dtb11u, &v, Ni11, &dtb11u );
      MultVector3d ( Ni02, &v, &dtb11v );

      CrossProduct3d ( &pv, &puu, &v );
      AddVector3Md ( &dtb11, &v, Ni10, &dtb11 );
      MultVector3d ( Ni11, &v, &z );
      MultVector3d ( -Ni02, &v, &y );

      CrossProduct3d ( &pv, &puv, &v );
      AddVector3Md ( &dtb12, &v, Ni10, &dtb12 );
      AddVector3Md ( &z, &v, -Ni20, &z );
      AddVector3Md ( &dtb22v, &v, -Ni02, &dtb22v );

      CrossProduct3d ( &puv, &pu, &v );
      AddVector3Md ( &dtb12, &v, Ni01, &dtb12 );
      AddVector3Md ( &dtb11u, &v, -Ni20, &dtb11u );
      MultVector3d ( Ni02, &v, &w );

      CrossProduct3d ( &pv, &pvv, &v );
      AddVector3Md ( &dtb22, &v, Ni10, &dtb22 );
      AddVector3Md ( &dtb22v, &v, Ni11, &dtb22v );
      AddVector3Md ( &y, &v, Ni20, &y );

      CrossProduct3d ( &pvv, &pu, &v );
      AddVector3Md ( &dtb22, &v, Ni01, &dtb22 );
      AddVector3Md ( &dtb11v, &v, -Ni20, &dtb11v );
      AddVector3Md ( &w, &v, -Ni11, &w );

      CrossProduct3d ( &pv, &puuu, &v );
      AddVector3Md ( &dtb11u, &v, Ni10, &dtb11u );
      CrossProduct3d ( &puuu, &pu, &v );
      AddVector3Md ( &dtb11u, &v, Ni01, &dtb11u );

      CrossProduct3d ( &puv, &puu, &v );
      AddVector3Md ( &dtb11u, &v, Ni10, &dtb11u );
      AddVector3Md ( &z, &v, -Ni01, &z );

      CrossProduct3d ( &pv, &puuv, &v );
      AddVector3Md ( &dtb12u, &v, Ni10, &dtb12u );

      CrossProduct3d ( &puuv, &pu, &v );
      AddVector3Md ( &dtb12u, &v, Ni01, &dtb12u );

      CrossProduct3d ( &pvv, &puv, &v );
      AddVector3Md ( &w, &v, Ni10, &w );
      AddVector3Md ( &dtb22v, &v, Ni01, &dtb22v );

      CrossProduct3d ( &pv, &puvv, &v );
      AddVector3Md ( &dtb22u, &v, Ni10, &dtb22u );

      CrossProduct3d ( &puvv, &pu, &v );
      AddVector3Md ( &dtb22u, &v, Ni01, &dtb22u );
      AddVector3d ( &dtb22u, &w, &dtb12v );
      SubtractPoints3d ( &dtb22u, &w, &dtb22u );

      AddVector3d ( &dtb11v, &z, &dtb11v );
      AddVector3d ( &dtb12u, &dtb11v, &dtb11v );
      SubtractPoints3d ( &dtb12u, &z, &dtb12u );

      CrossProduct3d ( &pv, &pvvv, &v );
      AddVector3Md ( &dtb22v, &v, Ni10, &dtb22v );

      CrossProduct3d ( &pvvv, &pu, &v );
      AddVector3Md ( &dtb22v, &v, Ni01, &dtb22v );

      CrossProduct3d ( &pvv, &puu, &v );
      AddVector3Md ( &dtb22u, &v, Ni01, &dtb22u );
      AddVector3Md ( &dtb11v, &v, Ni10, &dtb11v );

      AddVector3d ( &dtb22u, &y, &dtb22u );

      AddVector3d ( &puuu, &puvv, &w );
      CrossProduct3d ( &pv, &w, &v );
      AddVector3Md ( &Dbstar[9], &v, Ni10, &Dbstar[9] );
      CrossProduct3d ( &w, &pu, &v );
      AddVector3Md ( &Dbstar[9], &v, Ni01, &Dbstar[9] );

      AddVector3d ( &puuv, &pvvv, &w );
      CrossProduct3d ( &pv, &w, &v );
      AddVector3Md ( &Dbstar[10], &v, Ni10, &Dbstar[10] );
      CrossProduct3d ( &w, &pu, &v );
      AddVector3Md ( &Dbstar[10], &v, Ni01, &Dbstar[10] );
    }
} /*_g2mbl_UCompDBStard*/

void _g2mbl_UFuncGradSQIntegrandd ( int nkn2, int nf, double *Nitab,
                                    int knot, int k, int l, vector3d pder[11],
                                    int *cpind, int *nncpi,
                                    double Gstar[9], double *DGstar,
                                    double Bstar[11], double *DBstar,
                                    double C,
                                    double *F, double *DF )
{
  double detG, detGu, detGv, tH, tHu, tHv;
  double DdetG[3], DdetGu[3], DdetGv[3], DtH[3], DtHu[3], DtHv[3];
  double L[2], M[3], N, DL[6], DM[9], DN[3];
  double MLT[2], LMLT, DMLT[6], DLMLT[3], LDMLT[3];
  double g, aa, bb, DA[3], BB, DBB[3];
  int    fn, fn3, p, q;
  double first, *DGf, *DBf, *Ni, Df[3], Dg[3];

        /* the first functional */
          /* compute L, M, N */
  M[2] = g11;
  M[1] = -g12;
  M[0] = g22; 
  detG = g11*g22 - g12*g12;
  N = detG*detG, N = 1.0/(N*N*detG*sqrt(detG));
  detGu = g11u*g22 + g11*g22u - 2.0*g12*g12u;  
  detGv = g11v*g22 + g11*g22v - 2.0*g12*g12v;  
  tH = g11*tb22 + tb11*g22 - 2.0*g12*tb12;     
  tHu = g11u*tb22 + g11*tb22u + tb11u*g22 + tb11*g22u - 2.0*(g12u*tb12 + g12*tb12u);
  tHv = g11v*tb22 + g11*tb22v + tb11v*g22 + tb11*g22v - 2.0*(g12v*tb12 + g12*tb12v);
  L[0] = detG*tHu - 1.5*detGu*tH;
  L[1] = detG*tHv - 1.5*detGv*tH;
  MLT[0] = M[0]*L[0] + M[1]*L[1];
  MLT[1] = M[1]*L[0] + M[2]*L[1];
  LMLT = L[0]*MLT[0] + L[1]*MLT[1];   
  first = 0.25*LMLT*N;
        /* the second functional */
          /* Frobenius norm of the Laplacian gradient */
  if ( C > 0.0 ) {
    g = DotProduct3d ( &pder[9], &pder[9] ) +
        DotProduct3d ( &pder[10], &pder[10] );
          /* Frobenius norm of the projection of the Laplacian gradient */
          /* on the normal vector */
    BB = Bstar[9]*Bstar[9] + Bstar[10]*Bstar[10];
    g -= BB/detG;
    *F = first + C*g;
  }
  else
    *F = first;
        /* compute the functional derivatives */
  for ( fn = fn3 = 0;  fn < nf;  fn++, fn3 += 3 )
    if ( nncpi[cpind[fn]] >= 0 ) {
      DGf = &DGstar[9*fn3];
      DBf = &DBstar[11*fn3];
          /* derivatives of det G */
      pkn_MVectorLinCombd ( 3, 3, DdetG,
                &DGf[0], g22, &DGf[6], g11, &DGf[3], -2.0*g12 );
          /* derivatives of (det G)_u */
      pkn_MVectorLinCombd ( 6, 3, DdetGu,
                &DGf[9], g22, &DGf[6], g11u, &DGf[0], g22u,
                &DGf[15], g11, &DGf[3], -2.0*g12u, &DGf[12], -2.0*g12 );
          /* derivatives of (det G)_u */
      pkn_MVectorLinCombd ( 6, 3, DdetGv,
                &DGf[18], g22, &DGf[6], g11v, &DGf[0], g22v,
                &DGf[24], g11, &DGf[3], -2.0*g12v, &DGf[21], -2.0*g12 );
          /* derivatives of tH */
      pkn_MVectorLinCombd ( 6, 3, DtH,
                &DGf[0], tb22, &DBf[6], g11, &DBf[0], g22,
                &DGf[6], tb11, &DGf[3], -2.0*tb12, &DBf[3], -2.0*g12 );
          /* derivatives of tH_u */
      pkn_MVectorLinCombd ( 12, 3, DtHu,
                &DGf[9], tb22, &DBf[6], g11u, &DGf[0], tb22u,
                &DBf[15], g11, &DBf[9], g22, &DGf[6], tb11u, 
                &DBf[0], g22u, &DGf[15], tb11, &DGf[12], -2.0*tb12,
                &DBf[3], -2.0*g12u, &DGf[3], -2.0*tb12u, &DBf[12], -2.0*g12 );
          /* derivatives of tH_v */
      pkn_MVectorLinCombd ( 12, 3, DtHv,
                &DGf[18], tb22, &DBf[6], g11v, &DGf[0], tb22v,
                &DBf[24], g11, &DBf[18], g22, &DGf[6], tb11v, 
                &DBf[0], g22v, &DGf[24], tb11, &DGf[21], -2.0*tb12,
                &DBf[3], -2.0*g12v, &DGf[3], -2.0*tb12v, &DBf[21], -2.0*g12 );
          /* derivatives of N */
      MultVector3d ( -5.5*N/detG, (vector3d*)DdetG, (vector3d*)&DN[0] );
          /* derivatives of M */
      memcpy ( &DM[0], &DGf[6], sizeof(vector3d) );
      SetVector3d ( (vector3d*)&DM[3], -DGf[3], -DGf[4], -DGf[5] );
      memcpy ( &DM[6], DGf, sizeof(vector3d) );
          /* derivatives of L */
      pkn_MVectorLinCombd ( 4, 3, &DL[0],
                DdetG, tHu, DtHu, detG, DdetGu, -1.5*tH, DtH, -1.5*detGu );
      pkn_MVectorLinCombd ( 4, 3, &DL[3],
                DdetG, tHv, DtHv, detG, DdetGv, -1.5*tH, DtH, -1.5*detGv );
          /* compute derivatives of f */
      pkn_MVectorLinCombd ( 2, 3, &DLMLT[0], &DL[0], MLT[0], &DL[3], MLT[1] );
      pkn_MVectorLinCombd ( 2, 3, &DMLT[0], &DM[0], L[0], &DM[3], L[1] );
      pkn_MVectorLinCombd ( 2, 3, &DMLT[3], &DM[3], L[0], &DM[6], L[1] );
      pkn_MVectorLinCombd ( 2, 3, LDMLT, &DMLT[0], L[0], &DMLT[3], L[1] );
      pkn_MVectorLinCombd ( 3, 3, Df, DLMLT, 2.0*N, LDMLT, N, &DN[0], LMLT );
      if ( C > 0.0 ) {
          /* compute derivatives of g */
            /* derivative of the Frobenius norm of Laplacian gradient */
        if ( k == 4 || fn == 0 )
          Ni = &Nitab[(fn*nkn2+knot)*9];
        else {
          p = (fn-1) / 6;
          p = (l - p + k) % k;
          q = (fn-1) % 6;
          Ni = &Nitab[((6*p+q+1)*nkn2+knot)*9];
        }
        aa = Ni[5] + Ni[7];
        bb = Ni[6] + Ni[8];
        Dg[0] = 2.0*(aa*pder[9].x + bb*pder[10].x);
        Dg[1] = 2.0*(aa*pder[9].y + bb*pder[10].y);
        Dg[2] = 2.0*(aa*pder[9].z + bb*pder[10].z);
            /* derivative of the projection on the surface normal */
        MultVector3d ( -1.0/(detG*detG), (vector3d*)DdetG, (vector3d*)DA );
        AddVector3Md ( (vector3d*)Dg, (vector3d*)DA, -BB, (vector3d*)Dg );
        MultVector3d ( Bstar[9], (vector3d*)&DBf[9*3], (vector3d*)DBB );
        AddVector3Md ( (vector3d*)DBB, (vector3d*)&DBf[10*3], Bstar[10],
                       (vector3d*)DBB );
        AddVector3Md ( (vector3d*)Dg, (vector3d*)DBB, -2.0/detG, (vector3d*)Dg );
            /* add the derivatives of the two terms */
        MultVector3d ( 0.25, (vector3d*)Df, (vector3d*)Df );
        AddVector3Md ( (vector3d*)Df, (vector3d*)Dg, C, (vector3d*)&DF[fn3] );
      }
      else
        MultVector3d ( 0.25, (vector3d*)Df, (vector3d*)&DF[fn3] );
    }
} /*_g2mbl_UFuncGradSQIntegrandd*/

void g2mbl_UFuncGradRSQd ( int nkn, double *qcoeff, double *Nitab,
                           int *cpind, int *nncpi, point3d *mvcp,
                           double C, double *ftab, double *gtab )
{
  int      nkn2, i, j, knot;
  vector3d pder[11];
  double   Gstar[9], DGstar[9*3*16], Bstar[11], DBstar[11*3*16];
  double   sum0, sum1, f;
  double   Df[SQUAREGDS], sumD1[SQUAREGDS];

  nkn2 = nkn*nkn;
  memset ( Df, 0, SQUAREGDS*sizeof(double) );
  sum0 = 0.0;
  memset ( gtab, 0, SQUAREGDS*sizeof(double) );
  for ( i = knot = 0;  i < nkn;  i++ ) {
    sum1 = 0.0;
    memset ( sumD1, 0, SQUAREGDS*sizeof(double) );
    for ( j = 0;  j < nkn;  j++, knot++ ) {
        /* compute the parameterization derivatives */
      _g2mbl_UCompRDerd ( nkn2, Nitab, knot, cpind, mvcp, pder );
        /* compute the first fundamental form and its derivatives */
      _g2bl_UCompGStard ( pder, Gstar );
      _g2mbl_UCompDGStard ( nkn2, 16, Nitab, knot, 4, 0, pder, cpind, nncpi, DGstar );
        /* compute the second fundamental form and its derivatives */
      _g2bl_UCompBStard ( pder, Bstar );
      _g2mbl_UCompDBStard ( nkn2, 16, Nitab, knot, 4, 0, pder, cpind, nncpi, DBstar );
        /* compute the integrand terms */
      _g2mbl_UFuncGradSQIntegrandd ( nkn2, 16, Nitab, knot, 4, 0, pder,
                        cpind, nncpi, Gstar, DGstar, Bstar, DBstar, C, &f, Df );
        /* add the quadrature term */
      sum1 += f*qcoeff[j];
      pkn_AddMatrixMd ( 1, 3*16, 0, sumD1, 0, Df, qcoeff[j], 0, sumD1 );
    }
    sum0 += sum1*qcoeff[i];
    pkn_AddMatrixMd ( 1, 3*16, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
  }
  *ftab = sum0;
} /*g2mbl_UFuncGradRSQd*/

boolean g2mbl_UFuncGradSSQd ( int nkn, double *qcoeff, int k,
                              double *Nitab, double *Jac,
                              int *cpind, int *nncpi, point3d *mvcp,
                              double C, double *ftab, double *gtab )
{
  void     *sp;
  int      nkn2, i, j, l, knot, nf, nf3;
  vector3d *pder;
  double   *Gstar, *DGstar, *Bstar, *DBstar;
  double   sum0, sum1, f;
  double   *Df, *sumD1;

  sp = pkv_GetScratchMemTop ();
  nkn2 = nkn*nkn;
  nf = 6*k+1;
  nf3 = 3*nf;
  pder = (point3d*)pkv_GetScratchMemd ( 53 + 68*nf );
  if ( !pder )
    goto failure;
  Gstar  = &pder[11].x;
  DGstar = &Gstar[9];
  Bstar  = &DGstar[27*nf];
  DBstar = &Bstar[11];
  Df     = &DBstar[33*nf];
  sumD1  = &Df[nf3];

  memset ( Df, 0, nf3*sizeof(double) );
  sum0 = 0.0;
  memset ( gtab, 0, nf3*sizeof(double) );
  for ( l = 0; l < k; l++ ) {
    for ( i = knot = 0;  i < nkn;  i++ ) {
      sum1 = 0.0;
      memset ( sumD1, 0, nf3*sizeof(double) );
      for ( j = 0;  j < nkn;  j++, knot++ ) {
          /* compute the parameterization derivatives */
        _g2mbl_UCompSDerd ( nkn2, Nitab, knot, k, l, cpind, mvcp, pder );
          /* compute the first fundamental form and its derivatives */
        _g2bl_UCompGStard ( pder, Gstar );
        _g2mbl_UCompDGStard ( nkn2, nf, Nitab, knot, k, l, pder, cpind, nncpi, DGstar );
          /* compute the second fundamental form and its derivatives */
        _g2bl_UCompBStard ( pder, Bstar );
        _g2mbl_UCompDBStard ( nkn2, nf, Nitab, knot, k, l, pder, cpind, nncpi, DBstar );
          /* compute the integrand terms */
        _g2mbl_UFuncGradSQIntegrandd ( nkn2, nf, Nitab, knot, k, l, pder,
                          cpind, nncpi, Gstar, DGstar, Bstar, DBstar, C, &f, Df );
          /* add the quadrature term */
        sum1 += f*Jac[knot]*qcoeff[j];
        pkn_AddMatrixMd ( 1, nf3, 0, sumD1, 0, Df, Jac[knot]*qcoeff[j], 0, sumD1 );
      }
      sum0 += sum1*qcoeff[i];
      pkn_AddMatrixMd ( 1, nf3, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
    }
  }
  *ftab = sum0;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_UFuncGradSSQd*/

/* ///////////////////////////////////////////////////////////////////////// */
void _g2mbl_UCompDDGstard ( int nkn2, int nf, double *Nijtab,
                            int knot, int k, int l,
                            int *cpind, int *nncpi, double *DDGstar )
{
  int    fi, fj, fii, fjj, p, q;
  double *Nij;

  for ( fi = 0; fi < nf; fi++ )
    if ( nncpi[cpind[fi]] >= 0 ) {
      if ( k != 4 && fi > 0 ) {
        p = (fi-1) / 6;
        p = (l-p+k) % k;
        q = (fi-1) % 6;
        fii = 6*p+q+1;
      }
      else fii = fi;
      for ( fj = 0; fj <= fi; fj++ )
        if ( nncpi[cpind[fj]] >= 0 ) {
          if ( k != 4 && fj > 0 ) {
            p = (fj-1) / 6;
            p = (l-p+k) % k;
            q = (fj-1) % 6;
            fjj = 6*p+q+1;
          }
          else fjj = fj;
          Nij = &Nijtab[(pkn_SymMatIndex(fii,fjj)*nkn2+knot)*9];
          memcpy ( &DDGstar[9*pkn_SymMatIndex(fi,fj)], Nij, 9*sizeof(double) );
        }
    }
} /*_g2mbl_UCompDDGstard*/

void _g2mbl_UCompDDBstard ( int nkn2, int nf, double *Mijtab, int knot, int k, int l,
                            vector3d *pder, int *cpind, int *nncpi,
                            double *DDBstar )
{
#define ddtb11  DDB[0]
#define ddtb12  DDB[1]
#define ddtb22  DDB[2]
#define ddtb11u DDB[3]
#define ddtb12u DDB[4]
#define ddtb22u DDB[5]
#define ddtb11v DDB[6]
#define ddtb12v DDB[7]
#define ddtb22v DDB[8]
#define ddtlg1  DDB[9]
#define ddtlg2  DDB[10]
  int      fi, fj, fii, fjj, p, q;
  double   *Mij, negMij[18];
  vector3d *DDB;
  vector3d DDtb301001, DDtb211001, DDtb121001, DDtb031001, DDtb201011,
           DDtb112001, DDtb022001, DDtb201002, DDtb111002, DDtb021101;

  for ( fi = 0; fi < nf; fi++ )
    if ( nncpi[cpind[fi]] >= 0 ) {
      if ( k == 4 || fi == 0 )
        fii = fi;
      else {
        p = (fi-1) / 6;
        p = (l-p+k) % k;
        q = (fi-1) % 6;
        fii = 6*p+q+1;
      }
      DDB = (vector3d*)&DDBstar[11*3*pkn_SymMatIndex(fi,fi)];
      memset ( DDB, 0, 11*3*sizeof(double) );
      for ( fj = fi+1; fj < nf; fj++ )
        if ( nncpi[cpind[fj]] >= 0 ) {
          if ( k == 4 )
            fjj = fj;
          else {
            p = (fj-1) / 6;
            p = (l-p+k) % k;
            q = (fj-1) % 6;
            fjj = 6*p+q+1;
          }
          DDB = (vector3d*)&DDBstar[11*3*pkn_SymMatIndex(fi,fj)];
          if ( fii < fjj )
            Mij = &Mijtab[(pkn_SymMatIndex(fii,fjj-1)*nkn2+knot)*18];
          else {
            Mij = &Mijtab[(pkn_SymMatIndex(fii-1,fjj)*nkn2+knot)*18];
            for ( p = 0; p < 18; p++ )
              negMij[p] = -Mij[p];
            Mij = &negMij[0];
          }
            /* DDtb301001 */
          MultVector3d ( Mij1001, &puuu, &DDtb301001 );
          AddVector3Md ( &DDtb301001, &pu, Mij0130, &DDtb301001 );
          AddVector3Md ( &DDtb301001, &pv, Mij3010, &DDtb301001 );
            /* DDtb211001 */
          MultVector3d ( Mij1001, &puuv, &DDtb211001 );
          AddVector3Md ( &DDtb211001, &pu, Mij0121, &DDtb211001 );
          AddVector3Md ( &DDtb211001, &pv, Mij2110, &DDtb211001 );
            /* DDtb121001 */
          MultVector3d ( Mij1001, &puvv, &DDtb121001 );
          AddVector3Md ( &DDtb121001, &pu, Mij0112, &DDtb121001 );
          AddVector3Md ( &DDtb121001, &pv, Mij1210, &DDtb121001 );
            /* DDtb031001 */
          MultVector3d ( Mij1001, &pvvv, &DDtb031001 );
          AddVector3Md ( &DDtb031001, &pu, Mij0103, &DDtb031001 );
          AddVector3Md ( &DDtb031001, &pv, Mij0310, &DDtb031001 );
            /* DDtb201011 */
          MultVector3d ( Mij1011, &puu, &DDtb201011 );
          AddVector3Md ( &DDtb201011, &pu, Mij1120, &DDtb201011 );
          AddVector3Md ( &DDtb201011, &puv, Mij2010, &DDtb201011 );
            /* DDtb112001 */
          MultVector3d ( Mij2001, &puv, &DDtb112001 );
          AddVector3Md ( &DDtb112001, &puu, Mij0111, &DDtb112001 );
          AddVector3Md ( &DDtb112001, &pv, Mij1120, &DDtb112001 );
            /* DDtb022001 */
          MultVector3d ( Mij2001, &pvv, &DDtb022001 );
          AddVector3Md ( &DDtb022001, &puu, Mij0102, &DDtb022001 );
          AddVector3Md ( &DDtb022001, &pv, Mij0220, &DDtb022001 );
            /* DDtb201002 */
          MultVector3d ( Mij1002, &puu, &DDtb201002 );
          AddVector3Md ( &DDtb201002, &pu, Mij0220, &DDtb201002 );
          AddVector3Md ( &DDtb201002, &pvv, Mij2010, &DDtb201002 );
            /* DDtb111002 */ 
          MultVector3d ( Mij1002, &puv, &DDtb111002 );
          AddVector3Md ( &DDtb111002, &pu, Mij0211, &DDtb111002 );
          AddVector3Md ( &DDtb111002, &pvv, Mij1110, &DDtb111002 );
            /* DDtb021101 */
          MultVector3d ( Mij1101, &pvv, &DDtb021101 );
          AddVector3Md ( &DDtb021101, &puv, Mij0102, &DDtb021101 );
          AddVector3Md ( &DDtb021101, &pv, Mij0211, &DDtb021101 );
            /* DDtb11 */
          MultVector3d ( Mij1001, &puu, &ddtb11 );
          AddVector3Md ( &ddtb11, &pu, Mij0120, &ddtb11 );
          AddVector3Md ( &ddtb11, &pv, Mij2010, &ddtb11 );
            /* DDtb12 */
          MultVector3d ( Mij1001, &puv, &ddtb12 );
          AddVector3Md ( &ddtb12, &pu, Mij0111, &ddtb12 );
          AddVector3Md ( &ddtb12, &pv, Mij1110, &ddtb12 );
            /* DDtb22 */
          MultVector3d ( Mij1001, &pvv, &ddtb22 );
          AddVector3Md ( &ddtb22, &pu, Mij0102, &ddtb22 );
          AddVector3Md ( &ddtb22, &pv, Mij0210, &ddtb22 );
            /* DDtb11u */
          AddVector3d ( &DDtb301001, &DDtb201011, &ddtb11u );
            /* DDtb12u */
          AddVector3d ( &DDtb211001, &DDtb112001, &ddtb12u );
            /* DDtb22u */
          AddVector3d ( &DDtb121001, &DDtb022001, &ddtb22u );
          SubtractPoints3d ( &ddtb22u, &DDtb111002, &ddtb22u );
            /* DDtb11v */
          SubtractPoints3d ( &DDtb211001, &DDtb112001, &ddtb11v );
          AddVector3d ( &ddtb11v, &DDtb201002, &ddtb11v );
            /* DDtb12v */
          AddVector3d ( &DDtb121001, &DDtb111002, &ddtb12v );
            /* DDtb22v */
          AddVector3d ( &DDtb031001, &DDtb021101, &ddtb22v );
            /* DD(h100130+h100112) */
          AddVector3d ( &DDtb301001, &DDtb121001, &ddtlg1 );
            /* DD(h100121+h100103) */
          AddVector3d ( &DDtb211001, &DDtb031001, &ddtlg2 );
        }
    }
#undef ddtb11
#undef ddtb12
#undef ddtb22
#undef ddtb11u
#undef ddtb12u
#undef ddtb22u
#undef ddtb11v
#undef ddtb12v
#undef ddtb22v
#undef ddtlg1
#undef ddtlg2
} /*_g2mbl_UCompDDBstard*/

boolean _g2mbl_UFuncGradHessSQIntegrandd ( int nkn2, int nf, double *Nitab,
                                int knot, int k, int l, vector3d *pder,
                                int *cpind, int *nncpi,
                                double *Gstar, double *DGstar, double *DDGstar,
                                double *Bstar, double *DBstar, double *DDBstar,
                                double C,
                                double *F, double *DF, double *DDF )
{
  void   *sp;
  double detG, detGu, detGv, tH, tHu, tHv;
  double *DdetG, *DdetGu, *DdetGv, *DtH, *DtHu, *DtHv;
  double L[2], M[3], N, *DL, *DM, *DN;
  double DDdetG[9], DDdetGu[9], DDdetGv[9], DDtH[9], DDtHu[9], DDtHv[9];
  double DDL[18], DDM[3], DDN[9], DDaux[9];
  double MLT[2], LMLT, DMLT[6], DLMLT[3], LDMLT[3];
  double a, g, aa, bb, cc, dd, *DA, DDA[9], BB, *DBB, aux[3];
  int    fi, fi3, fj, fj3, p, q, ll, mm, nn;
  double *Na, *Nb, *DGa, *DGb, *DBa, *DBb, *gstar, *bstar;
  double first, Df[3], Dg[3], DDf[9], DDg[9];

  sp = pkv_GetScratchMemTop ();
  DdetG = pkv_GetScratchMemd ( 14*3*nf );
  if ( !DdetG )
    goto failure;
  DdetGu = &DdetG[3*nf];
  DdetGv = &DdetGu[3*nf];
  DtH    = &DdetGv[3*nf];
  DtHu   = &DtH[3*nf];
  DtHv   = &DtHu[3*nf];
  DL     = &DtHv[3*nf];
  DM     = &DL[6*nf];
  DN     = &DM[9*nf];
  DA     = &DN[3*nf];
  DBB    = &DA[3*nf];
        /* the first functional */
          /* compute L, M, N */
  M[2] = g11;
  M[1] = -g12;
  M[0] = g22; 
  detG = g11*g22 - g12*g12;
  N = detG*detG, N = 1.0/(N*N*detG*sqrt(detG));
  detGu = g11u*g22 + g11*g22u - 2.0*g12*g12u;  
  detGv = g11v*g22 + g11*g22v - 2.0*g12*g12v;  
  tH = g11*tb22 + tb11*g22 - 2.0*g12*tb12;     
  tHu = g11u*tb22 + g11*tb22u + tb11u*g22 + tb11*g22u - 2.0*(g12u*tb12 + g12*tb12u);
  tHv = g11v*tb22 + g11*tb22v + tb11v*g22 + tb11*g22v - 2.0*(g12v*tb12 + g12*tb12v);
  L[0] = detG*tHu - 1.5*detGu*tH;
  L[1] = detG*tHv - 1.5*detGv*tH;
  MLT[0] = M[0]*L[0] + M[1]*L[1];
  MLT[1] = M[1]*L[0] + M[2]*L[1];
  LMLT = L[0]*MLT[0] + L[1]*MLT[1];   
  first = 0.25*LMLT*N;
  if ( C > 0.0 ) {
        /* the second functional */
          /* Frobenius norm of the Laplacian gradient */
    g = DotProduct3d ( &pder[9], &pder[9] ) + DotProduct3d ( &pder[10], &pder[10] );
          /* Frobenius norm of the projection of the Laplacian gradient */
          /* on the normal vector */
    BB = Bstar[9]*Bstar[9] + Bstar[10]*Bstar[10];
    *F = first + C*(g - BB/detG);
  }
  else
    *F = first;
        /* compute the functional derivatives */
  for ( fi = fi3 = 0;  fi < nf;  fi++, fi3 += 3 )
    if ( nncpi[cpind[fi]] >= 0 ) {
      DGa = &DGstar[9*fi3];
      DBa = &DBstar[11*fi3];
          /* derivatives of det G */
      pkn_MVectorLinCombd ( 3, 3, &DdetG[fi3],
                &DGa[0], g22, &DGa[6], g11, &DGa[3], -2.0*g12 );
          /* derivatives of (det G)_u */
      pkn_MVectorLinCombd ( 6, 3, &DdetGu[fi3],
                &DGa[9], g22, &DGa[6], g11u, &DGa[0], g22u,
                &DGa[15], g11, &DGa[3], -2.0*g12u, &DGa[12], -2.0*g12 );
          /* derivatives of (det G)_u */
      pkn_MVectorLinCombd ( 6, 3, &DdetGv[fi3],
                &DGa[18], g22, &DGa[6], g11v, &DGa[0], g22v,
                &DGa[24], g11, &DGa[3], -2.0*g12v, &DGa[21], -2.0*g12 );
          /* derivatives of tH */
      pkn_MVectorLinCombd ( 6, 3, &DtH[fi3],
                &DGa[0], tb22, &DBa[6], g11, &DBa[0], g22,
                &DGa[6], tb11, &DGa[3], -2.0*tb12, &DBa[3], -2.0*g12 );
          /* derivatives of tH_u */
      pkn_MVectorLinCombd ( 12, 3, &DtHu[fi3],
                &DGa[9], tb22, &DBa[6], g11u, &DGa[0], tb22u,
                &DBa[15], g11, &DBa[9], g22, &DGa[6], tb11u, 
                &DBa[0], g22u, &DGa[15], tb11, &DGa[12], -2.0*tb12,
                &DBa[3], -2.0*g12u, &DGa[3], -2.0*tb12u, &DBa[12], -2.0*g12 );
          /* derivatives of tH_v */
      pkn_MVectorLinCombd ( 12, 3, &DtHv[fi3],
                &DGa[18], tb22, &DBa[6], g11v, &DGa[0], tb22v,
                &DBa[24], g11, &DBa[18], g22, &DGa[6], tb11v, 
                &DBa[0], g22v, &DGa[24], tb11, &DGa[21], -2.0*tb12,
                &DBa[3], -2.0*g12v, &DGa[3], -2.0*tb12v, &DBa[21], -2.0*g12 );
          /* derivatives of N */
      MultVector3d ( -5.5*N/detG, (vector3d*)&DdetG[fi3], (vector3d*)&DN[fi3] );
          /* derivatives of M */
      memcpy ( &DM[3*fi3], &DGa[6], sizeof(vector3d) );
      SetVector3d ( (vector3d*)&DM[3*fi3+3], -DGa[3], -DGa[4], -DGa[5] );
      memcpy ( &DM[3*fi3+6], &DGstar[9*fi3], sizeof(vector3d) );
          /* derivatives of L */
      pkn_MVectorLinCombd ( 4, 3, &DL[2*fi3],
                &DdetG[fi3], tHu, &DtHu[fi3], detG,
                &DdetGu[fi3], -1.5*tH, &DtH[fi3], -1.5*detGu );
      pkn_MVectorLinCombd ( 4, 3, &DL[2*fi3+3],
                &DdetG[fi3], tHv, &DtHv[fi3], detG,
                &DdetGv[fi3], -1.5*tH, &DtH[fi3], -1.5*detGv );
          /* compute derivatives of f */
      pkn_MVectorLinCombd ( 2, 3, &DLMLT[0], &DL[2*fi3], MLT[0], &DL[2*fi3+3], MLT[1] );
      pkn_MVectorLinCombd ( 2, 3, &DMLT[0], &DM[3*fi3], L[0], &DM[3*fi3+3], L[1] );
      pkn_MVectorLinCombd ( 2, 3, &DMLT[3], &DM[3*fi3+3], L[0], &DM[3*fi3+6], L[1] );
      pkn_MVectorLinCombd ( 2, 3, LDMLT, &DMLT[0], L[0], &DMLT[3], L[1] );
      pkn_MVectorLinCombd ( 3, 3, Df, DLMLT, 2.0*N, LDMLT, N, &DN[fi3], LMLT );
      if ( C > 0.0 ) {
          /* compute derivatives of g */
            /* derivative of the Frobenius norm of Laplacian gradient */
        if ( k == 4 || fi == 0 )
          Na = &Nitab[(fi*nkn2+knot)*9];
        else {
          p = (fi-1) / 6;
          p = (l - p + k) % k;
          q = (fi-1) % 6;
          Na = &Nitab[((6*p+q+1)*nkn2+knot)*9];
        }
        aa = Na[5] + Na[7];
        bb = Na[6] + Na[8];
        Dg[0] = 2.0*(aa*pder[9].x + bb*pder[10].x);
        Dg[1] = 2.0*(aa*pder[9].y + bb*pder[10].y);
        Dg[2] = 2.0*(aa*pder[9].z + bb*pder[10].z);
            /* derivative of the projection on the surface normal */
        MultVector3d ( -1.0/(detG*detG), (vector3d*)&DdetG[fi3], (vector3d*)&DA[fi3] );
        AddVector3Md ( (vector3d*)Dg, (vector3d*)&DA[fi3], -BB,
                       (vector3d*)Dg );
        MultVector3d ( Bstar[9], (vector3d*)&DBa[9*3], (vector3d*)&DBB[fi3] );
        AddVector3Md ( (vector3d*)&DBB[fi3], (vector3d*)&DBa[10*3], Bstar[10],
                       (vector3d*)&DBB[fi3] );
        AddVector3Md ( (vector3d*)Dg, (vector3d*)&DBB[fi3], -2.0/detG,
                       (vector3d*)Dg );
            /* add the derivatives of the two terms */
        MultVector3d ( 0.25, (vector3d*)Df, (vector3d*)Df );
        AddVector3Md ( (vector3d*)Df, (vector3d*)Dg, C, (vector3d*)&DF[fi3] );
      }
      else
        MultVector3d ( 0.25, (vector3d*)Df, (vector3d*)&DF[fi3] );
    }

        /* compute the second order derivatives */
  for ( fi = fi3 = 0;  fi < nf;  fi++, fi3 += 3 )
    if ( nncpi[cpind[fi]] >= 0 ) {
      if ( k == 4 || fi == 0 )
        Na = &Nitab[(fi*nkn2+knot)*9];
      else {
        p = (fi-1) / 6;
        p = (l - p + k) % k;
        q = (fi-1) % 6;
        Na = &Nitab[((6*p+q+1)*nkn2+knot)*9];
      }
      DGa = &DGstar[9*fi3];
      DBa = &DBstar[11*fi3];
      for ( fj = fj3 = 0;  fj <= fi;  fj++, fj3 += 3 )
        if ( nncpi[cpind[fj]] >= 0 ) {
          if ( k == 4 || fj == 0 )
            Nb = &Nitab[(fj*nkn2+knot)*9];
          else {
            p = (fj-1) / 6;
            p = (l - p + k) % k;
            q = (fj-1) % 6;
            Nb = &Nitab[((6*p+q+1)*nkn2+knot)*9];
          }
          DGb = &DGstar[9*fj3];
          DBb = &DBstar[11*fj3];
          gstar = &DDGstar[9*pkn_SymMatIndex(fi,fj)];
          bstar = &DDBstar[3*11*pkn_SymMatIndex(fi,fj)];

          DDdetG[0] = DDdetG[4] = DDdetG[8] =
            g22*gstar[0] + g11*gstar[2] - 2.0*g12*gstar[1];
          DDdetG[0] += DGb[0]*DGa[6] + DGb[6]*DGa[0] - 2.0*DGb[3]*DGa[3];
          DDdetG[1] = DGb[0]*DGa[7] + DGb[6]*DGa[1] - 2.0*DGb[3]*DGa[4];
          DDdetG[2] = DGb[0]*DGa[8] + DGb[6]*DGa[2] - 2.0*DGb[3]*DGa[5];
          DDdetG[3] = DGb[1]*DGa[6] + DGb[7]*DGa[0] - 2.0*DGb[4]*DGa[3];
          DDdetG[4] += DGb[1]*DGa[7] + DGb[7]*DGa[1] - 2.0*DGb[4]*DGa[4];
          DDdetG[5] = DGb[1]*DGa[8] + DGb[7]*DGa[2] - 2.0*DGb[4]*DGa[5];
          DDdetG[6] = DGb[2]*DGa[6] + DGb[8]*DGa[0] - 2.0*DGb[5]*DGa[3];
          DDdetG[7] = DGb[2]*DGa[7] + DGb[8]*DGa[1] - 2.0*DGb[5]*DGa[4];
          DDdetG[8] += DGb[2]*DGa[8] + DGb[8]*DGa[2] - 2.0*DGb[5]*DGa[5];

          DDdetGu[0] = DDdetGu[4] = DDdetGu[8] =
                       g22*gstar[3] + g11u*gstar[2] + g22u*gstar[0] +
                       g11*gstar[5] - 2.0*(g12*gstar[4] + g12u*gstar[1]);
          DDdetGu[0] += DGb[9]*DGa[6] + DGb[6]*DGa[9] + DGb[0]*DGa[15] +
                        DGb[15]*DGa[0] - 2.0*(DGb[3]*DGa[12] + DGb[12]*DGa[3]);
          DDdetGu[1] = DGb[9]*DGa[7] + DGb[6]*DGa[10] + DGb[0]*DGa[16] +
                       DGb[15]*DGa[1] - 2.0*(DGb[3]*DGa[13] + DGb[12]*DGa[4]);
          DDdetGu[2] = DGb[9]*DGa[8] + DGb[6]*DGa[11] + DGb[0]*DGa[17] +
                       DGb[15]*DGa[2] - 2.0*(DGb[3]*DGa[14] + DGb[12]*DGa[5]);
          DDdetGu[3] = DGb[10]*DGa[6] + DGb[7]*DGa[9] + DGb[1]*DGa[15] +
                       DGb[16]*DGa[0] - 2.0*(DGb[4]*DGa[12] + DGb[13]*DGa[3]);
          DDdetGu[4] += DGb[10]*DGa[7] + DGb[7]*DGa[10] + DGb[1]*DGa[16] +
                        DGb[16]*DGa[1] - 2.0*(DGb[4]*DGa[13] + DGb[13]*DGa[4]);
          DDdetGu[5] = DGb[10]*DGa[8] + DGb[7]*DGa[11] + DGb[1]*DGa[17] +
                       DGb[16]*DGa[2] - 2.0*(DGb[4]*DGa[14] + DGb[13]*DGa[5]);
          DDdetGu[6] = DGb[11]*DGa[6] + DGb[8]*DGa[9] + DGb[2]*DGa[15] +
                       DGb[17]*DGa[0] - 2.0*(DGb[5]*DGa[12] + DGb[14]*DGa[3]);
          DDdetGu[7] = DGb[11]*DGa[7] + DGb[8]*DGa[10] + DGb[2]*DGa[16] +
                       DGb[17]*DGa[1] - 2.0*(DGb[5]*DGa[13] + DGb[14]*DGa[4]);
          DDdetGu[8] += DGb[11]*DGa[8] + DGb[8]*DGa[11] + DGb[2]*DGa[17] +
                        DGb[17]*DGa[2] - 2.0*(DGb[5]*DGa[14] + DGb[14]*DGa[5]);

          DDdetGv[0] = DDdetGv[4] = DDdetGv[8] =
                       g22*gstar[6] + g11v*gstar[2] + g22v*gstar[0] +
                       g11*gstar[8] - 2.0*(g12*gstar[7] + g12v*gstar[1]);
          DDdetGv[0] += DGb[18]*DGa[6] + DGb[6]*DGa[18] + DGb[0]*DGa[24] +
                        DGb[24]*DGa[0] - 2.0*(DGb[3]*DGa[21] + DGb[21]*DGa[3]);
          DDdetGv[1] = DGb[18]*DGa[7] + DGb[6]*DGa[19] + DGb[0]*DGa[25] +
                       DGb[24]*DGa[1] - 2.0*(DGb[3]*DGa[22] + DGb[21]*DGa[4]);
          DDdetGv[2] = DGb[18]*DGa[8] + DGb[6]*DGa[20] + DGb[0]*DGa[26] +
                       DGb[24]*DGa[2] - 2.0*(DGb[3]*DGa[23] + DGb[21]*DGa[5]);
          DDdetGv[3] = DGb[19]*DGa[6] + DGb[7]*DGa[18] + DGb[1]*DGa[24] +
                       DGb[25]*DGa[0] - 2.0*(DGb[4]*DGa[21] + DGb[22]*DGa[3]);
          DDdetGv[4] += DGb[19]*DGa[7] + DGb[7]*DGa[19] + DGb[1]*DGa[25] +
                        DGb[25]*DGa[1] - 2.0*(DGb[4]*DGa[22] + DGb[22]*DGa[4]);
          DDdetGv[5] = DGb[19]*DGa[8] + DGb[7]*DGa[20] + DGb[1]*DGa[26] +
                       DGb[25]*DGa[2] - 2.0*(DGb[4]*DGa[23] + DGb[22]*DGa[5]);
          DDdetGv[6] = DGb[20]*DGa[6] + DGb[8]*DGa[18] + DGb[2]*DGa[24] +
                       DGb[26]*DGa[0] - 2.0*(DGb[5]*DGa[21] + DGb[23]*DGa[3]);
          DDdetGv[7] = DGb[20]*DGa[7] + DGb[8]*DGa[19] + DGb[2]*DGa[25] +
                       DGb[26]*DGa[1] - 2.0*(DGb[5]*DGa[22] + DGb[23]*DGa[4]);
          DDdetGv[8] += DGb[20]*DGa[8] + DGb[8]*DGa[20] + DGb[2]*DGa[26] +
                        DGb[26]*DGa[2] - 2.0*(DGb[5]*DGa[23] + DGb[23]*DGa[5]);

          DDtH[0] = DDtH[4] = DDtH[8] =
                        tb22*gstar[0] + tb11*gstar[2] - 2.0*tb12*gstar[1];
          a = g11*bstar[8] + g22*bstar[2] - 2.0*g12*bstar[5]; 
          DDtH[1] = a;   DDtH[3] = -a;
          a = g11*bstar[7] + g22*bstar[1] - 2.0*g12*bstar[4];
          DDtH[2] = -a;  DDtH[6] = a;
          a = g11*bstar[6] + g22*bstar[0] - 2.0*g12*bstar[3];
          DDtH[5] = a;   DDtH[7] = -a;
          DDtH[0] += DGb[0]*DBa[6] + DGb[6]*DBa[0] + DBb[0]*DGa[6] + DBb[6]*DGa[0] -
                     2.0*(DGb[3]*DBa[3] + DBb[3]*DGa[3]);
          DDtH[1] += DGb[0]*DBa[7] + DGb[6]*DBa[1] + DBb[0]*DGa[7] + DBb[6]*DGa[1] -
                     2.0*(DGb[3]*DBa[4] + DBb[3]*DGa[4]);
          DDtH[2] += DGb[0]*DBa[8] + DGb[6]*DBa[2] + DBb[0]*DGa[8] + DBb[6]*DGa[2] -
                     2.0*(DGb[3]*DBa[5] + DBb[3]*DGa[5]);
          DDtH[3] += DGb[1]*DBa[6] + DGb[7]*DBa[0] + DBb[1]*DGa[6] + DBb[7]*DGa[0] -
                     2.0*(DGb[4]*DBa[3] + DBb[4]*DGa[3]);
          DDtH[4] += DGb[1]*DBa[7] + DGb[7]*DBa[1] + DBb[1]*DGa[7] + DBb[7]*DGa[1] -
                     2.0*(DGb[4]*DBa[4] + DBb[4]*DGa[4]);
          DDtH[5] += DGb[1]*DBa[8] + DGb[7]*DBa[2] + DBb[1]*DGa[8] + DBb[7]*DGa[2] -
                     2.0*(DGb[4]*DBa[5] + DBb[4]*DGa[5]);
          DDtH[6] += DGb[2]*DBa[6] + DGb[8]*DBa[0] + DBb[2]*DGa[6] + DBb[8]*DGa[0] -
                     2.0*(DGb[5]*DBa[3] + DBb[5]*DGa[3]);
          DDtH[7] += DGb[2]*DBa[7] + DGb[8]*DBa[1] + DBb[2]*DGa[7] + DBb[8]*DGa[1] -
                     2.0*(DGb[5]*DBa[4] + DBb[5]*DGa[4]);
          DDtH[8] += DGb[2]*DBa[8] + DGb[8]*DBa[2] + DBb[2]*DGa[8] + DBb[8]*DGa[2] -
                     2.0*(DGb[5]*DBa[5] + DBb[5]*DGa[5]);

          DDtHu[0] = DDtHu[4] = DDtHu[8] =
                tb22*gstar[3] + tb22u*gstar[0] + tb11*gstar[5] + tb11u*gstar[2] -
                2.0*(tb12*gstar[4] + tb12u*gstar[1]);
          a = g11u*bstar[8] + g11*bstar[17] + g22*bstar[11] + g22u*bstar[2] -
              2.0*(g12u*bstar[5] + g12*bstar[14]);
          DDtHu[1] = a;   DDtHu[3] = -a;
          a = g11u*bstar[7] + g11*bstar[16] + g22*bstar[10] + g22u*bstar[1] -
              2.0*(g12u*bstar[4] + g12*bstar[13]);
          DDtHu[2] = -a;  DDtHu[6] = a;
          a = g11u*bstar[6] + g11*bstar[15] + g22*bstar[9] + g22u*bstar[0] -
              2.0*(g12u*bstar[3] + g12*bstar[12]);
          DDtHu[5] = a;   DDtHu[7] = -a;
          DDtHu[0] += DGb[9]*DBa[6] + DGb[0]*DBa[15] + DGb[6]*DBa[9] + DGb[15]*DBa[0] +
                      DBb[9]*DGa[6] + DBb[0]*DGa[15] + DBb[6]*DGa[9] + DBb[15]*DGa[0] -
                      2.0*(DGb[12]*DBa[3] + DGb[3]*DBa[12] +
                           DBb[12]*DGa[3] + DBb[3]*DGa[12]);
          DDtHu[1] += DGb[9]*DBa[7] + DGb[0]*DBa[16] + DGb[6]*DBa[10] + DGb[15]*DBa[1] +
                      DBb[9]*DGa[7] + DBb[0]*DGa[16] + DBb[6]*DGa[10] + DBb[15]*DGa[1] -
                      2.0*(DGb[12]*DBa[4] + DGb[3]*DBa[13] +
                           DBb[12]*DGa[4] + DBb[3]*DGa[13]);
          DDtHu[2] += DGb[9]*DBa[8] + DGb[0]*DBa[17] + DGb[6]*DBa[11] + DGb[15]*DBa[2] +
                      DBb[9]*DGa[8] + DBb[0]*DGa[17] + DBb[6]*DGa[11] + DBb[15]*DGa[2] -
                      2.0*(DGb[12]*DBa[5] + DGb[3]*DBa[14] +
                           DBb[12]*DGa[5] + DBb[3]*DGa[14]);
          DDtHu[3] += DGb[10]*DBa[6] + DGb[1]*DBa[15] + DGb[7]*DBa[9] + DGb[16]*DBa[0] +
                      DBb[10]*DGa[6] + DBb[1]*DGa[15] + DBb[7]*DGa[9] + DBb[16]*DGa[0] -
                      2.0*(DGb[13]*DBa[3] + DGb[4]*DBa[12] +
                           DBb[13]*DGa[3] + DBb[4]*DGa[12]);
          DDtHu[4] += DGb[10]*DBa[7] + DGb[1]*DBa[16] + DGb[7]*DBa[10] + DGb[16]*DBa[1] +
                      DBb[10]*DGa[7] + DBb[1]*DGa[16] + DBb[7]*DGa[10] + DBb[16]*DGa[1] -
                      2.0*(DGb[13]*DBa[4] + DGb[4]*DBa[13] +
                           DBb[13]*DGa[4] + DBb[4]*DGa[13]);
          DDtHu[5] += DGb[10]*DBa[8] + DGb[1]*DBa[17] + DGb[7]*DBa[11] + DGb[16]*DBa[2] +
                      DBb[10]*DGa[8] + DBb[1]*DGa[17] + DBb[7]*DGa[11] + DBb[16]*DGa[2] -
                      2.0*(DGb[13]*DBa[5] + DGb[4]*DBa[14] +
                           DBb[13]*DGa[5] + DBb[4]*DGa[14]);
          DDtHu[6] += DGb[11]*DBa[6] + DGb[2]*DBa[15] + DGb[8]*DBa[9] + DGb[17]*DBa[0] +
                      DBb[11]*DGa[6] + DBb[2]*DGa[15] + DBb[8]*DGa[9] + DBb[17]*DGa[0] -
                      2.0*(DGb[14]*DBa[3] + DGb[5]*DBa[12] +
                           DBb[14]*DGa[3] + DBb[5]*DGa[12]);
          DDtHu[7] += DGb[11]*DBa[7] + DGb[2]*DBa[16] + DGb[8]*DBa[10] + DGb[17]*DBa[1] +
                      DBb[11]*DGa[7] + DBb[2]*DGa[16] + DBb[8]*DGa[10] + DBb[17]*DGa[1] -
                      2.0*(DGb[14]*DBa[4] + DGb[5]*DBa[13] +
                           DBb[14]*DGa[4] + DBb[5]*DGa[13]);
          DDtHu[8] += DGb[11]*DBa[8] + DGb[2]*DBa[17] + DGb[8]*DBa[11] + DGb[17]*DBa[2] +
                      DBb[11]*DGa[8] + DBb[2]*DGa[17] + DBb[8]*DGa[11] + DBb[17]*DGa[2] -
                      2.0*(DGb[14]*DBa[5] + DGb[5]*DBa[14] +
                           DBb[14]*DGa[5] + DBb[5]*DGa[14]);

          DDtHv[0] = DDtHv[4] = DDtHv[8] =
                tb22*gstar[6] + tb22v*gstar[0] + tb11*gstar[8] + tb11v*gstar[2] -
                2.0*(tb12*gstar[7] + tb12v*gstar[1]);
          a = g11v*bstar[8] + g11*bstar[26] + g22*bstar[20] + g22v*bstar[2] -
              2.0*(g12v*bstar[5] + g12*bstar[23]);
          DDtHv[1] = a;   DDtHv[3] = -a;
          a = g11v*bstar[7] + g11*bstar[25] + g22*bstar[19] + g22v*bstar[1] -
              2.0*(g12v*bstar[4] + g12*bstar[22]);
          DDtHv[2] = -a;  DDtHv[6] = a;
          a = g11v*bstar[6] + g11*bstar[24] + g22*bstar[18] + g22v*bstar[0] -
              2.0*(g12v*bstar[3] + g12*bstar[21]);
          DDtHv[5] = a;   DDtHv[7] = -a;
          DDtHv[0] += DGb[18]*DBa[6] + DGb[0]*DBa[24] + DGb[6]*DBa[18] + DGb[24]*DBa[0] +
                      DBb[18]*DGa[6] + DBb[0]*DGa[24] + DBb[6]*DGa[18] + DBb[24]*DGa[0] -
                      2.0*(DGb[21]*DBa[3] + DGb[3]*DBa[21] +
                           DBb[21]*DGa[3] + DBb[3]*DGa[21]);
          DDtHv[1] += DGb[18]*DBa[7] + DGb[0]*DBa[25] + DGb[6]*DBa[19] + DGb[24]*DBa[1] +
                      DBb[18]*DGa[7] + DBb[0]*DGa[25] + DBb[6]*DGa[19] + DBb[24]*DGa[1] -
                      2.0*(DGb[21]*DBa[4] + DGb[3]*DBa[22] +
                           DBb[21]*DGa[4] + DBb[3]*DGa[22]);
          DDtHv[2] += DGb[18]*DBa[8] + DGb[0]*DBa[26] + DGb[6]*DBa[20] + DGb[24]*DBa[2] +
                      DBb[18]*DGa[8] + DBb[0]*DGa[26] + DBb[6]*DGa[20] + DBb[24]*DGa[2] -
                      2.0*(DGb[21]*DBa[5] + DGb[3]*DBa[23] +
                           DBb[21]*DGa[5] + DBb[3]*DGa[23]);
          DDtHv[3] += DGb[19]*DBa[6] + DGb[1]*DBa[24] + DGb[7]*DBa[18] + DGb[25]*DBa[0] +
                      DBb[19]*DGa[6] + DBb[1]*DGa[24] + DBb[7]*DGa[18] + DBb[25]*DGa[0] -
                      2.0*(DGb[22]*DBa[3] + DGb[4]*DBa[21] +
                           DBb[22]*DGa[3] + DBb[4]*DGa[21]);
          DDtHv[4] += DGb[19]*DBa[7] + DGb[1]*DBa[25] + DGb[7]*DBa[19] + DGb[25]*DBa[1] +
                      DBb[19]*DGa[7] + DBb[1]*DGa[25] + DBb[7]*DGa[19] + DBb[25]*DGa[1] -
                      2.0*(DGb[22]*DBa[4] + DGb[4]*DBa[22] +
                           DBb[22]*DGa[4] + DBb[4]*DGa[22]);
          DDtHv[5] += DGb[19]*DBa[8] + DGb[1]*DBa[26] + DGb[7]*DBa[20] + DGb[25]*DBa[2] +
                      DBb[19]*DGa[8] + DBb[1]*DGa[26] + DBb[7]*DGa[20] + DBb[25]*DGa[2] -
                      2.0*(DGb[22]*DBa[5] + DGb[4]*DBa[23] +
                           DBb[22]*DGa[5] + DBb[4]*DGa[23]);
          DDtHv[6] += DGb[20]*DBa[6] + DGb[2]*DBa[24] + DGb[8]*DBa[18] + DGb[26]*DBa[0] +
                      DBb[20]*DGa[6] + DBb[2]*DGa[24] + DBb[8]*DGa[18] + DBb[26]*DGa[0] -
                      2.0*(DGb[23]*DBa[3] + DGb[5]*DBa[21] +
                           DBb[23]*DGa[3] + DBb[5]*DGa[21]);
          DDtHv[7] += DGb[20]*DBa[7] + DGb[2]*DBa[25] + DGb[8]*DBa[19] + DGb[26]*DBa[1] +
                      DBb[20]*DGa[7] + DBb[2]*DGa[25] + DBb[8]*DGa[19] + DBb[26]*DGa[1] -
                      2.0*(DGb[23]*DBa[4] + DGb[5]*DBa[22] +
                           DBb[23]*DGa[4] + DBb[5]*DGa[22]);
          DDtHv[8] += DGb[20]*DBa[8] + DGb[2]*DBa[26] + DGb[8]*DBa[20] + DGb[26]*DBa[2] +
                      DBb[20]*DGa[8] + DBb[2]*DGa[26] + DBb[8]*DGa[20] + DBb[26]*DGa[2] -
                      2.0*(DGb[23]*DBa[5] + DGb[5]*DBa[23] +
                           DBb[23]*DGa[5] + DBb[5]*DGa[23]);

          /* the second order derivatives of N */
          pkn_MultMatrixd ( 3, 1, 1, &DdetG[fj3], 3, 0, &DdetG[fi3], 3, DDN );
          pkn_MatrixLinCombd ( 1, 9, 0, DDN, 13.0, 0, DDdetG, -2.0*detG, 0, DDN );
          pkn_MultMatrixNumd ( 1, 9, 0, DDN, 11.0*N/(4.0*detG*detG), 0, DDN );
          /* the second order derivatives of M */
            /* the 3 x 3 matrices of the second order derivatives of the */
            /* three coefficients of M are diagonal, with all diagonal */
            /* coefficients equal, so only one number is needed to represent */
            /* each of them */
          DDM[0] = gstar[2];
          DDM[1] = -gstar[1];
          DDM[2] = gstar[0];
          /* the second order derivatives of L */
          pkn_MatrixLinCombd ( 1, 9, 0, DDdetG, tHu, 0, DDtHu, detG, 0, &DDL[0] );
          pkn_MultMatrixAddd ( 3, 1, 1, &DdetG[fj3], 3, 0, &DtHu[fi3], 3, &DDL[0] );
          pkn_MultMatrixAddd ( 3, 1, 1, &DtHu[fj3], 3, 0, &DdetG[fi3], 3, &DDL[0] );
          pkn_MatrixLinCombd ( 1, 9, 0, DDdetGu, tH, 0, DDtH, detGu, 0, DDaux );
          pkn_MultMatrixAddd ( 3, 1, 1, &DdetGu[fj3], 3, 0, &DtH[fi3], 3, DDaux );
          pkn_MultMatrixAddd ( 3, 1, 1, &DtH[fj3], 3, 0, &DdetGu[fi3], 3, DDaux );
          pkn_AddMatrixMd ( 1, 9, 0, &DDL[0], 0, DDaux, -1.5, 0, &DDL[0] );

          pkn_MatrixLinCombd ( 1, 9, 0, DDdetG, tHv, 0, DDtHv, detG, 0, &DDL[9] );
          pkn_MultMatrixAddd ( 3, 1, 1, &DdetG[fj3], 3, 0, &DtHv[fi3], 3, &DDL[9] );
          pkn_MultMatrixAddd ( 3, 1, 1, &DtHv[fj3], 3, 0, &DdetG[fi3], 3, &DDL[9] );
          pkn_MatrixLinCombd ( 1, 9, 0, DDdetGv, tH, 0, DDtH, detGv, 0, DDaux );
          pkn_MultMatrixAddd ( 3, 1, 1, &DdetGv[fj3], 3, 0, &DtH[fi3], 3, DDaux );
          pkn_MultMatrixAddd ( 3, 1, 1, &DtH[fj3], 3, 0, &DdetGv[fi3], 3, DDaux );
          pkn_AddMatrixMd ( 1, 9, 0, &DDL[9], 0, DDaux, -1.5, 0, &DDL[9] );

          /* the second order derivatives of f */
          for ( ll = nn = 0;  ll < 3;  ll++ )
            for ( mm = 0;  mm < 3;  mm++, nn++ ) {
              if ( ll == mm )
                a = L[0]*(DDM[0]*L[0]+DDM[1]*L[1]) + L[1]*(DDM[1]*L[0]+DDM[2]*L[1]);
              else
                a = 0.0;
              DDf[3*mm+ll] =
                  (2.0*(DDL[nn]*MLT[0]+DDL[9+nn]*MLT[1] +
                        DL[6*fj+ll]*(DM[9*fi+mm]*L[0]+DM[9*fi+3+mm]*L[1]) +
                        DL[6*fj+3+ll]*(DM[9*fi+3+mm]*L[0]+DM[9*fi+6+mm]*L[1]) +
                        DL[6*fi+mm]*(DM[9*fj+ll]*L[0]+DM[9*fj+3+ll]*L[1]) +
                        DL[6*fi+3+mm]*(DM[9*fj+3+ll]*L[0]+DM[9*fj+6+ll]*L[1]) +
                        DL[6*fj+ll]*(M[0]*DL[6*fi+mm]+M[1]*DL[6*fi+3+mm]) +
                        DL[6*fj+3+ll]*(M[1]*DL[6*fi+mm]+M[2]*DL[6*fi+3+mm])) + a)*N +
                  (2.0*(DL[6*fj+ll]*MLT[0]+DL[6*fj+3+ll]*MLT[1]) +
                        L[0]*(DM[9*fj+ll]*L[0]+DM[9*fj+3+ll]*L[1])+
                        L[1]*(DM[9*fj+3+ll]*L[0]+DM[9*fj+6+ll]*L[1]))*DN[fi3+mm] +
                  (2.0*(DL[6*fi+mm]*MLT[0]+DL[6*fi+3+mm]*MLT[1])+
                        L[0]*(DM[9*fi+mm]*L[0]+DM[9*fi+3+mm]*L[1])+
                        L[1]*(DM[9*fi+3+mm]*L[0]+DM[9*fi+6+mm]*L[1]))*DN[fj3+ll] +
                  LMLT*DDN[nn];
            }
          if ( C > 0.0 ) {
          /* compute the second order derivatives of g */
            /* Frobenius norm of the Laplacian gradient */
            memset ( DDg, 0, 9*sizeof(double) );
            aa = Na[5] + Na[7];
            bb = Na[6] + Na[8];
            cc = Nb[5] + Nb[7];
            dd = Nb[6] + Nb[8];
            DDg[0] = DDg[4] = DDg[8] = 2.0*(aa*cc + bb*dd);
            /* projection on the surface normal */
            pkn_MultMatrixd ( 3, 1, 1, &DdetG[fj3], 3, 3, &DdetG[fi3], 3, DDA );
            pkn_MatrixLinCombd ( 1, 9, 0, DDA, 2.0/(detG*detG*detG),
                                 0, DDdetG, -1.0/(detG*detG), 0, DDA );
            pkn_AddMatrixMd ( 1, 9, 0, DDg, 0, DDA, -BB, 0, DDg );

            pkn_MultMatrixd ( 3, 1, 1, (double*)&DBB[fj3],
                              3, 3, (double*)&DA[fi3], 3, DDaux );
            pkn_AddMatrixMd ( 1, 9, 0, DDg, 0, DDaux, -2.0, 0, DDg );
            pkn_MultMatrixd ( 3, 1, 1, (double*)&DA[fj3],
                              3, 3, (double*)&DBB[fi3], 3, DDaux );
            pkn_AddMatrixMd ( 1, 9, 0, DDg, 0, DDaux, -2.0, 0, DDg );

            pkn_MatrixLinCombd ( 1, 3, 0, &bstar[9*3], Bstar[9],
                                 0, &bstar[10*3], Bstar[10], 0, aux );
            pkn_MultMatrixd ( 3, 1, 1, &DBb[9*3], 3, 3, &DBa[9*3], 3, DDaux );
            pkn_MultMatrixAddd ( 3, 1, 1, &DBb[10*3], 3, 3, &DBa[10*3], 3, DDaux );
            DDaux[1] += aux[2];  DDaux[3] -= aux[2];
            DDaux[2] -= aux[1];  DDaux[6] += aux[1];
            DDaux[5] += aux[0];  DDaux[7] -= aux[0];
            pkn_AddMatrixMd ( 1, 9, 0, DDg, 0, DDaux, -2.0/detG, 0, DDaux );
              /* for historic reasons this is the transposition of the */
              /* necessary matrix, to be fixed later, now just transpose it */
            pkv_TransposeMatrixd ( 3, 3, 3, DDaux, 3, DDg );
              /* add the derivatives of the two terms */
            pkn_MatrixLinCombd ( 1, 9, 0, DDf, 0.25, 0, DDg, C, 0, DDaux );
          }
          else
            pkn_MultMatrixNumd ( 1, 9, 0, DDf, 0.25, 0, DDaux );
          if ( nncpi[cpind[fi]] < nncpi[cpind[fj]] )
            pkv_TransposeMatrixd ( 3, 3, 3, DDaux, 3, &DDF[9*pkn_SymMatIndex(fi,fj)] );
          else
            memcpy ( &DDF[9*pkn_SymMatIndex(fi,fj)], DDaux, 9*sizeof(double) );
/* DEBUG - to remove */
/*
memset ( &DDF[9*pkn_SymMatIndex(fi,fj)], 0, 9*sizeof(double) );
DDF[9*pkn_SymMatIndex(fi,fj)] = 
DDF[9*pkn_SymMatIndex(fi,fj)+4] = 
DDF[9*pkn_SymMatIndex(fi,fj)+8] = DDM[2];
*/
/*
DDF[9*pkn_SymMatIndex(fi,fj)] = 
DDF[9*pkn_SymMatIndex(fi,fj)+4] = 
DDF[9*pkn_SymMatIndex(fi,fj)+8] = 0.0;
if ( nncpi[cpind[fi]] < nncpi[cpind[fj]] ) {
  DDF[9*pkn_SymMatIndex(fi,fj)+1] =  bstar[24+2];
  DDF[9*pkn_SymMatIndex(fi,fj)+3] = -bstar[24+2];
  DDF[9*pkn_SymMatIndex(fi,fj)+2] = -bstar[24+1];
  DDF[9*pkn_SymMatIndex(fi,fj)+6] =  bstar[24+1];
  DDF[9*pkn_SymMatIndex(fi,fj)+5] =  bstar[24+0];
  DDF[9*pkn_SymMatIndex(fi,fj)+7] = -bstar[24+0];
}
else {
  DDF[9*pkn_SymMatIndex(fi,fj)+1] = -bstar[24+2];
  DDF[9*pkn_SymMatIndex(fi,fj)+3] =  bstar[24+2];
  DDF[9*pkn_SymMatIndex(fi,fj)+2] =  bstar[24+1];
  DDF[9*pkn_SymMatIndex(fi,fj)+6] = -bstar[24+1];
  DDF[9*pkn_SymMatIndex(fi,fj)+5] = -bstar[24+0];
  DDF[9*pkn_SymMatIndex(fi,fj)+7] =  bstar[24+0];
}
*/
/*
if ( nncpi[cpind[fi]] > nncpi[cpind[fj]] )
  pkv_TransposeMatrixd ( 3, 3, 3, &DDL[9], 3, &DDF[9*pkn_SymMatIndex(fi,fj)] );
else
  memcpy ( &DDF[9*pkn_SymMatIndex(fi,fj)], &DDL[9], 9*sizeof(double) );
*/
      }
    }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_UFuncGradHessSQIntegrandd*/

boolean g2mbl_UFuncGradHessRSQd ( int nkn, double *qcoeff,
                                  double *Nitab, double *Nijtab, double *Mijtab,
                                  int *cpind, int *nncpi, point3d *mvcp,
                                  double C,
                                  double *ftab, double *gtab, double *htab )
{
  int      nkn2, i, j, knot;
  vector3d pder[11];
  double   Gstar[9], DGstar[9*3*16], DDGstar[9*136],
           Bstar[11], DBstar[11*3*16], DDBstar[11*3*136];
  double   sum0, sum1, f;
  double   Df[SQUAREGDS], sumD1[SQUAREGDS], DDf[SQUAREHDS], sumDD1[SQUAREHDS];

  nkn2 = nkn*nkn;
  memset ( Df, 0, SQUAREGDS*sizeof(double) );
  memset ( DDf, 0, SQUAREHDS*sizeof(double) );
  sum0 = 0.0;
  memset ( gtab, 0, SQUAREGDS*sizeof(double) );
  memset ( htab, 0, SQUAREHDS*sizeof(double) );
  for ( i = knot = 0;  i < nkn;  i++ ) {
    sum1 = 0.0;
    memset ( sumD1, 0, SQUAREGDS*sizeof(double) );
    memset ( sumDD1, 0, SQUAREHDS*sizeof(double) );
    for ( j = 0;  j < nkn;  j++, knot++ ) {
        /* compute the parameterization derivatives */
      _g2mbl_UCompRDerd ( nkn2, Nitab, knot, cpind, mvcp, pder );
        /* compute the first fundamental form and its derivatives */
      _g2bl_UCompGStard ( pder, Gstar );
      _g2mbl_UCompDGStard ( nkn2, 16, Nitab, knot, 4, 0, pder, cpind, nncpi, DGstar );
      _g2mbl_UCompDDGstard ( nkn2, 16, Nijtab, knot, 4, 0, cpind, nncpi, DDGstar );
        /* compute the second fundamental form and its derivatives */
      _g2bl_UCompBStard ( pder, Bstar );
      _g2mbl_UCompDBStard ( nkn2, 16, Nitab, knot, 4, 0, pder, cpind, nncpi, DBstar );
      _g2mbl_UCompDDBstard ( nkn2, 16, Mijtab, knot, 4, 0,
                             pder, cpind, nncpi, DDBstar );

        /* compute the integrand terms */
      if ( !_g2mbl_UFuncGradHessSQIntegrandd ( nkn2, 16, Nitab, knot, 4, 0, pder,
                                         cpind, nncpi, Gstar, DGstar, DDGstar,
                                         Bstar, DBstar, DDBstar, C,
                                         &f, Df, DDf ) )
        return false;
        /* add the quadrature term */
      sum1 += f*qcoeff[j];
      pkn_AddMatrixMd ( 1, 3*16, 0, sumD1, 0, Df, qcoeff[j], 0, sumD1 );
      pkn_AddMatrixMd ( 1, 9*8*17, 0, sumDD1, 0, DDf, qcoeff[j], 0, sumDD1 );
    }
    sum0 += sum1*qcoeff[i];
    pkn_AddMatrixMd ( 1, 3*16, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
    pkn_AddMatrixMd ( 1, 9*8*17, 0, htab, 0, sumDD1, qcoeff[i], 0, htab );
  }
  *ftab = sum0;
  return true;
} /*g2mbl_UFuncGradHessRSQd*/

boolean g2mbl_UFuncGradHessSSQd ( int nkn, double *qcoeff, int k,
                                  double *Nitab, double *Nijtab, double *Mijtab,
                                  double *Jac,
                                  int *cpind, int *nncpi, point3d *mvcp,
                                  double C,
                                  double *ftab, double *gtab, double *htab )
{
  void     *sp;
  int      nkn2, i, j, l, knot, nf, nf3, nh, nh9;
  vector3d *pder;
  double   *Gstar, *DGstar, *DDGstar, *Bstar, *DBstar, *DDBstar;
  double   sum0, sum1, f;
  double   *Df, *sumD1, *DDf, *sumDD1;

  sp = pkv_GetScratchMemTop ();
  nkn2 = nkn*nkn;
  nf = 6*k+1;
  nf3 = 3*nf;
  nh = (nf*(nf+1))/2;
  nh9 = 9*nh;
  pder = (point3d*)pkv_GetScratchMemd ( 53 + 12*nf3 + 69*nh );
  if ( !pder )
    goto failure;
  Gstar = &pder[11].x;
  DGstar = &Gstar[9];
  DDGstar = &DGstar[9*nf3];
  Bstar = &DDGstar[nh9];
  DBstar = &Bstar[11];
  DDBstar = &DBstar[11*nf3];
  Df  = &DDBstar[11*3*nh];
  DDf = &Df[nf3];
  sumD1 = &DDf[nh9];
  sumDD1 = &sumD1[nf3];

  memset ( Df, 0, nf3*sizeof(double) );
  memset ( DDf, 0, nh9*sizeof(double) );
  sum0 = 0.0;
  memset ( gtab, 0, nf3*sizeof(double) );
  memset ( htab, 0, nh9*sizeof(double) );
  for ( l = 0; l < k; l++ ) {
    for ( i = knot = 0;  i < nkn;  i++ ) {
      sum1 = 0.0;
      memset ( sumD1, 0, nf3*sizeof(double) );
      memset ( sumDD1, 0, nh9*sizeof(double) );
      for ( j = 0;  j < nkn;  j++, knot++ ) {
        /* compute the parameterization derivatives */
        _g2mbl_UCompSDerd ( nkn2, Nitab, knot, k, l, cpind, mvcp, pder );
        /* compute the first fundamental form and its derivatives */
        _g2bl_UCompGStard ( pder, Gstar );
        _g2mbl_UCompDGStard ( nkn2, nf, Nitab, knot, k, l, pder, cpind, nncpi, DGstar );
        _g2mbl_UCompDDGstard ( nkn2, nf, Nijtab, knot, k, l, cpind, nncpi, DDGstar );
        /* compute the second fundamental form and its derivatives */
        _g2bl_UCompBStard ( pder, Bstar );
        _g2mbl_UCompDBStard ( nkn2, nf, Nitab, knot, k, l, pder, cpind, nncpi, DBstar );
        _g2mbl_UCompDDBstard ( nkn2, nf, Mijtab, knot, k, l,
                               pder, cpind, nncpi, DDBstar );
        /* compute the integrand terms */
        if ( !_g2mbl_UFuncGradHessSQIntegrandd ( nkn2, nf, Nitab, knot, k, l, pder,
                                           cpind, nncpi, Gstar, DGstar, DDGstar,
                                           Bstar, DBstar, DDBstar, C,
                                           &f, Df, DDf ) )
          goto failure;
        /* add the quadrature term */
        sum1 += f*Jac[knot]*qcoeff[j];
        pkn_AddMatrixMd ( 1, nf3, 0, sumD1, 0, Df, Jac[knot]*qcoeff[j], 0, sumD1 );
        pkn_AddMatrixMd ( 1, nh9, 0, sumDD1, 0, DDf, Jac[knot]*qcoeff[j], 0, sumDD1 );
      }
      sum0 += sum1*qcoeff[i];
      pkn_AddMatrixMd ( 1, nf3, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
      pkn_AddMatrixMd ( 1, nh9, 0, htab, 0, sumDD1, qcoeff[i], 0, htab );
    }
  }
  *ftab = sum0;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_UFuncGradHessSSQd*/

