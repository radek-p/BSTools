
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "g2blendingd.h"

#include "g2blprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
void _g2bl_UFuncSQIntegrandd ( vector3d pder[11], double *first, double *second )
{
  double   Gstar[9], detG, detGu, detGv;
  double   Bstar[11], tH, tHu, tHv;
  double   L[2], M[3], N, g;

        /* compute the first fundamental form and its derivatives */
  _g2bl_UCompGStard ( pder, Gstar );
        /* compute the second fundamental form and its derivatives */
  _g2bl_UCompBStard ( pder, Bstar );
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

  *first = 0.25*(L[0]*(L[0]*M[0]+L[1]*(M[1]+M[1]))+L[1]*L[1]*M[2])*N;
        /* the second functional */
          /* Frobenius norm of the Laplacian gradient */
  g = DotProduct3d ( &pder[9], &pder[9] ) + DotProduct3d ( &pder[10], &pder[10] );
          /* Frobenius norm of the projection of the Laplacian gradient */
          /* on the normal vector */
  *second = g - (Bstar[9]*Bstar[9] + Bstar[10]*Bstar[10])/detG;
} /*_g2bl_UFuncSQIntegrandd*/

void g2bl_UFuncSQd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv,
                     int pitch, point3d *cp,
                     double tC,
                     int isq, int jsq,
                     double *ftab )
{
  int      fcpn; /* first control point number */
  int      i, j;
  double   sum0, sum1, f, g;
  vector3d pder[11];

  fcpn = isq*pitch + jsq;
  sum0 = 0.0;
  for ( i = 0; i < nkn; i++ ) {
    sum1 = 0.0;
    for ( j = 0; j < nkn; j++ ) {
        /* compute the parameterization derivatives */
      _g2bl_UCompPDerd ( nkn, Nitab, pitch, cp, fcpn, i, j, pder );
        /* compute the integrand terms */
      _g2bl_UFuncSQIntegrandd ( pder, &f, &g );
        /* add the quadrature term */
      sum1 += (f + tC*g)*qcoeff[j];
    }
    sum0 += sum1*qcoeff[i];
  }
        /* store the integral over the square in the array */
  ftab[0] = sum0;
} /*g2bl_UFuncSQd*/

/* ///////////////////////////////////////////////////////////////////////// */
void g2bl_UFuncGradSQd ( int nkn, const double *qcoeff, double *Nitab,
                         int lastknotu, int lastknotv,
                         int pitch, point3d *cp,
                         double tC,
                         int isq, int jsq, int ip0, int ip1, int jp0, int jp1,
                         double *ftab, double *gtab )
{
  int      fcpn;
  vector3d pder[11];
  double   Gstar[9], DGstar[9*3*16], Bstar[11], DBstar[11*3*16];
  double   detG, detGu, detGv, tH, tHu, tHv;
  double   DdetG[3], DdetGu[3], DdetGv[3], DtH[3], DtHu[3], DtHv[3];
  double   L[2], M[3], N, DL[6], DM[9], DN[3];
  double   MLT[2], LMLT, DMLT[6], DLMLT[3], LDMLT[3];
  double   *DGk, *DBk;
  double   BB, aa, bb, *Ni, DA[3], DBB[3];
  double   sum0, sum1, f, g;
  double   Df[SQUAREGDS], Dg[SQUAREGDS], sumD1[SQUAREGDS];
  int      i, j, ii, jj, kk, kk3, ip, jp;

  fcpn = isq*pitch + jsq;
  memset ( Df, 0, SQUAREGDS*sizeof(double) );
  memset ( Dg, 0, SQUAREGDS*sizeof(double) );
  sum0 = 0.0;
  memset ( gtab, 0, SQUAREGDS*sizeof(double) );
  for ( i = 0; i < nkn; i++ ) {
    sum1 = 0.0;
    memset ( sumD1, 0, SQUAREGDS*sizeof(double) );
    for ( j = 0; j < nkn; j++ ) {
        /* compute the parameterization derivatives */
      _g2bl_UCompPDerd ( nkn, Nitab, pitch, cp, fcpn, i, j, pder );
        /* compute the first fundamental form and its derivatives */
      _g2bl_UCompGStard ( pder, Gstar );
      _g2bl_UCompDGStard ( nkn, Nitab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                           isq, jsq, i, j, pder, DGstar );
        /* compute the second fundamental form and its derivatives */
      _g2bl_UCompBStard ( pder, Bstar );
      _g2bl_UCompDBStard ( nkn, Nitab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                           isq, jsq, i, j, pder, DBstar );
        /* the first functional */
      M[2] = g11;
      M[1] = -g12;
      M[0] = g22; 
      detG = g11*g22 - g12*g12;
      N = detG*detG, N = 1.0/(N*N*detG*sqrt(detG));
      detGu = g11u*g22 + g11*g22u - 2.0*g12*g12u;  
      detGv = g11v*g22 + g11*g22v - 2.0*g12*g12v;  
      tH = g11*tb22 + tb11*g22 - 2.0*g12*tb12;     
      tHu = g11u*tb22 + g11*tb22u + tb11u*g22 + tb11*g22u -
            2.0*(g12u*tb12 + g12*tb12u);
      tHv = g11v*tb22 + g11*tb22v + tb11v*g22 + tb11*g22v -
            2.0*(g12v*tb12 + g12*tb12v);
      L[0] = detG*tHu - 1.5*detGu*tH;   
      L[1] = detG*tHv - 1.5*detGv*tH;   
      MLT[0] = M[0]*L[0] + M[1]*L[1];
      MLT[1] = M[1]*L[0] + M[2]*L[1];
      LMLT = L[0]*MLT[0] + L[1]*MLT[1];
      f = 0.25*LMLT*N;

      if ( tC > 0.0 ) {
        /* the second functional */
          /* Frobenius norm of the Laplacian gradient */
        g = DotProduct3d ( &pder[9], &pder[9] ) +
            DotProduct3d ( &pder[10], &pder[10] );
          /* Frobenius norm of the projection of the Laplacian gradient */
          /* on the normal vector */
        BB = Bstar[9]*Bstar[9] + Bstar[10]*Bstar[10];
        g -= BB/detG;
        /* add the quadrature term */
        sum1 += (f + tC*g)*qcoeff[j];
      }
      else
        sum1 += f*qcoeff[j];

        /* compute the functional derivatives */
      for ( ii = 0, ip = isq;  ii <  4;  ii++, ip++ )
        if ( ip >= ip0 && ip < ip1 )
          for ( jj = 0, jp = jsq;  jj < 4;  jj++, jp++ )
            if ( jp >= jp0 && jp < jp1 ) {
              kk = 4*ii + jj;
              kk3 = 3*kk;
              /*kp = ip*pitch + jp;*/
          /* compute derivatives of L, M, N with respect */
          /* to the coordinates of the point cp[kp] */
              DGk = &DGstar[3*9*kk];
              DBk = &DBstar[3*11*kk];
            /* derivatives of det G */
              pkn_MVectorLinCombd ( 3, 3, DdetG,
                      &DGk[0], g22, &DGk[6], g11, &DGk[3], -2.0*g12 );
            /* derivatives of (det G)_u */
              pkn_MVectorLinCombd ( 6, 3, DdetGu,
                      &DGk[9], g22, &DGk[6], g11u, &DGk[0], g22u,
                      &DGk[15], g11, &DGk[3], -2.0*g12u, &DGk[12], -2.0*g12 );
            /* derivatives of (det G)_v */
              pkn_MVectorLinCombd ( 6, 3, DdetGv,
                      &DGk[18], g22, &DGk[6], g11v, &DGk[0], g22v,
                      &DGk[24], g11, &DGk[3], -2.0*g12v, &DGk[21], -2.0*g12 );
            /* derivatives of tH */
              pkn_MVectorLinCombd ( 6, 3, DtH,
                      &DGk[0], tb22, &DBk[6], g11, &DBk[0], g22,
                      &DGk[6], tb11, &DGk[3], -2.0*tb12, &DBk[3], -2.0*g12 );
            /* derivatives of tH_u */
              pkn_MVectorLinCombd ( 12, 3, DtHu,
                      &DGk[9], tb22, &DBk[6], g11u, &DGk[0], tb22u,
                      &DBk[15], g11, &DBk[9], g22, &DGk[6], tb11u,
                      &DBk[0], g22u, &DGk[15], tb11, &DGk[12], -2.0*tb12,
                      &DBk[3], -2.0*g12u, &DGk[3], -2.0*tb12u, &DBk[12], -2.0*g12 );
            /* derivatives of tH_v */
              pkn_MVectorLinCombd ( 12, 3, DtHv,
                      &DGk[18], tb22, &DBk[6], g11v, &DGk[0], tb22v,
                      &DBk[24], g11, &DBk[18], g22, &DGk[6], tb11v,
                      &DBk[0], g22v, &DGk[24], tb11, &DGk[21], -2.0*tb12,
                      &DBk[3], -2.0*g12v, &DGk[3], -2.0*tb12v, &DBk[21], -2.0*g12 );
            /* derivatives of N */
              MultVector3d ( -5.5*N/detG, (vector3d*)DdetG, (vector3d*)&DN[0] );
            /* derivatives of M */
              memcpy ( &DM[0], &DGk[6], sizeof(vector3d) );
              SetVector3d ( (vector3d*)&DM[3], -DGk[3], -DGk[4], -DGk[5] );
              memcpy ( &DM[6], &DGstar[3*9*kk], sizeof(vector3d) );
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
              pkn_MVectorLinCombd ( 3, 3, &Df[kk3],
                                    DLMLT, 2.0*N, LDMLT, N, &DN[0], LMLT );
              if ( tC > 0.0 ) {
          /* compute derivatives of g */
            /* derivative of the Frobenius norm of Laplacian gradient */
                Ni = &Nitab[((kk*nkn + i)*nkn + j)*9];
                aa = Ni[5] + Ni[7];
                bb = Ni[6] + Ni[8];
                Dg[kk3]   = 2.0*(aa*pder[9].x + bb*pder[10].x);
                Dg[kk3+1] = 2.0*(aa*pder[9].y + bb*pder[10].y);
                Dg[kk3+2] = 2.0*(aa*pder[9].z + bb*pder[10].z);
            /* derivative of the projection on the surface normal */
                MultVector3d ( -1.0/(detG*detG), (vector3d*)DdetG, (vector3d*)DA );
                AddVector3Md ( (vector3d*)&Dg[kk3], (vector3d*)DA, -BB,
                               (vector3d*)&Dg[kk3] );
                MultVector3d ( Bstar[9], (vector3d*)&DBk[9*3], (vector3d*)DBB );
                AddVector3Md ( (vector3d*)DBB, (vector3d*)&DBk[10*3], Bstar[10],
                               (vector3d*)DBB );
                AddVector3Md ( (vector3d*)&Dg[kk3], (vector3d*)DBB, -2.0/detG,
                               (vector3d*)&Dg[kk3] );

              }
            }
        /* integrate the functional derivatives */
      if ( tC > 0.0 )
        pkn_MatrixLinCombd ( 1, 3*16, 0, Df, 0.25, 0, Dg, tC, 0, Df );
      pkn_AddMatrixMd ( 1, 3*16, 0, sumD1, 0, Df, qcoeff[j], 0, sumD1 );
    }
    sum0 += sum1*qcoeff[i];
    pkn_AddMatrixMd ( 1, 3*16, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
  }
  *ftab = sum0;
} /*g2bl_UFuncGradSQd*/

/* ///////////////////////////////////////////////////////////////////////// */
void g2bl_UFuncGradHessianSQd ( int nkn, const double *qcoeff, double *Nitab, 
                                double *Nijtab, double *Mijtab,
                                int lastknotu, int lastknotv,
                                int pitch, point3d *cp,
                                double tC,
                                int isq, int jsq, int ip0, int ip1, int jp0, int jp1,
                                double *ftab, double *gtab, double *htab )
{
  int      fcpn;
  vector3d pder[11];
  double   Gstar[9], DGstar[9*3*16], DDGstar[9*136];
  double   Bstar[11], DBstar[11*3*16], DDBstar[11*3*136];
  double   detG, detGu, detGv, tH, tHu, tHv;
  double   DdetG[3*16], DdetGu[3*16], DdetGv[3*16];
  double   DtH[3*16], DtHu[3*16], DtHv[3*16];
  double   DDdetG[9], DDdetGu[9], DDdetGv[9];
  double   DDtH[9], DDtHu[9], DDtHv[9];
  double   L[2], M[3], N, DL[6*16], DM[9*16], DN[3*16], DDL[18], DDM[3], DDN[9];
  double   MLT[2], LMLT, DMLT[6], DLMLT[3], LDMLT[3];
  double   *DGa, *DBa, *DGb, *DBb, *gstar, *bstar;
  double   BB, aa, bb, cc, dd;
  double   *Na, *Nb, DA[3*16], DDA[9], DBB[3*16];
  double   sum0, sum1, f, g, a;
  double   Df[SQUAREGDS], Dg[SQUAREGDS], sumD1[SQUAREGDS];
  double   DDf[9], DDg[9], sumDD1[SQUAREHDS];
  double   aux[3], DDaux[9];
  int      i, j, l, m, n, ddi;
  int      ia, ja, ka, ka3, iap, jap/*, kap*/;
  int      ib, jb, kb, kb3, ibp, jbp/*, kbp*/;

  fcpn = isq*pitch + jsq;
  memset ( Df, 0, SQUAREGDS*sizeof(double) );
  memset ( Dg, 0, SQUAREGDS*sizeof(double) );

  sum0 = 0.0;
  memset ( gtab, 0, SQUAREGDS*sizeof(double) );
  memset ( htab, 0, SQUAREHDS*sizeof(double) );
  for ( i = 0; i < nkn; i++ ) {
    sum1 = 0.0;
    memset ( sumD1, 0, SQUAREGDS*sizeof(double) );
    memset ( sumDD1, 0, SQUAREHDS*sizeof(double) );
    for ( j = 0; j < nkn; j++ ) {
        /* compute the parameterization derivatives */
      _g2bl_UCompPDerd ( nkn, Nitab, pitch, cp, fcpn, i, j, pder );
        /* compute the first fundamental form and its derivatives */
      _g2bl_UCompGStard ( pder, Gstar );
      _g2bl_UCompDGStard ( nkn, Nitab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                           isq, jsq, i, j, pder, DGstar );
      _g2bl_UCompDDGStard ( nkn, Nijtab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                            isq, jsq, i, j, DDGstar );
        /* compute the second fundamental form and its derivatives */
      _g2bl_UCompBStard ( pder, Bstar );
      _g2bl_UCompDBStard ( nkn, Nitab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                           isq, jsq, i, j, pder, DBstar );
      _g2bl_UCompDDBStard ( nkn, Mijtab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                            isq, jsq, i, j, pder, DDBstar );

        /* the first functional */
      M[2] = g11;
      M[1] = -g12;
      M[0] = g22; 
      detG = g11*g22 - g12*g12;
      N = detG*detG, N = 1.0/(N*N*detG*sqrt(detG));
      detGu = g11u*g22 + g11*g22u - 2.0*g12*g12u;  
      detGv = g11v*g22 + g11*g22v - 2.0*g12*g12v;  
      tH = g11*tb22 + tb11*g22 - 2.0*g12*tb12;     
      tHu = g11u*tb22 + g11*tb22u + tb11u*g22 + tb11*g22u -
            2.0*(g12u*tb12 + g12*tb12u);
      tHv = g11v*tb22 + g11*tb22v + tb11v*g22 + tb11*g22v -
            2.0*(g12v*tb12 + g12*tb12v);
      L[0] = detG*tHu - 1.5*detGu*tH;   
      L[1] = detG*tHv - 1.5*detGv*tH;   
      MLT[0] = M[0]*L[0] + M[1]*L[1];
      MLT[1] = M[1]*L[0] + M[2]*L[1];
      LMLT = L[0]*MLT[0] + L[1]*MLT[1];
      f = 0.25*LMLT*N;

      if ( tC > 0.0 ) {
        /* the second functional */
          /* Frobenius norm of the Laplacian gradient */
        g = DotProduct3d ( &pder[9], &pder[9] ) +
            DotProduct3d ( &pder[10], &pder[10] );
          /* Frobenius norm of the projection of the Laplacian gradient */
          /* on the normal vector */
        BB = Bstar[9]*Bstar[9] + Bstar[10]*Bstar[10];
        g -= BB/detG;
        /* add the quadrature term */
        sum1 += (f + tC*g)*qcoeff[j];
      }
      else {
        sum1 += f*qcoeff[j];
        BB = 0.0;
      }

        /* compute the functional derivatives */
      for ( ia = 0, iap = isq;  ia < 4;  ia++, iap++ )
        if ( iap >= ip0 && iap < ip1 )
          for ( ja = 0, jap = jsq;  ja < 4;  ja++, jap++ )
            if ( jap >= jp0 && jap < jp1 ) {
              ka = 4*ia + ja;
              ka3 = 3*ka;
              /*kap = iap*pitch + jap;*/
          /* compute derivatives of L, M, N with respect */
          /* to the coordinates of the point cp[kap] */
              DGa = &DGstar[9*ka3];
              DBa = &DBstar[11*ka3];
            /* derivatives of det G */
              pkn_MVectorLinCombd ( 3, 3, &DdetG[ka3],
                      &DGa[0], g22, &DGa[6], g11, &DGa[3], -2.0*g12 );
            /* derivatives of (det G)_u */
              pkn_MVectorLinCombd ( 6, 3, &DdetGu[ka3],
                      &DGa[9], g22, &DGa[6], g11u, &DGa[0], g22u,
                      &DGa[15], g11, &DGa[3], -2.0*g12u, &DGa[12], -2.0*g12 );
            /* derivatives of (det G)_v */
              pkn_MVectorLinCombd ( 6, 3, &DdetGv[ka3],
                      &DGa[18], g22, &DGa[6], g11v, &DGa[0], g22v,
                      &DGa[24], g11, &DGa[3], -2.0*g12v, &DGa[21], -2.0*g12 );
            /* derivatives of tH */
              pkn_MVectorLinCombd ( 6, 3, &DtH[ka3],
                      &DGa[0], tb22, &DBa[6], g11, &DBa[0], g22,
                      &DGa[6], tb11, &DGa[3], -2.0*tb12, &DBa[3], -2.0*g12 );
            /* derivatives of tH_u */
              pkn_MVectorLinCombd ( 12, 3, &DtHu[ka3],
                      &DGa[9], tb22, &DBa[6], g11u, &DGa[0], tb22u,
                      &DBa[15], g11, &DBa[9], g22, &DGa[6], tb11u,
                      &DBa[0], g22u, &DGa[15], tb11, &DGa[12], -2.0*tb12,
                      &DBa[3], -2.0*g12u, &DGa[3], -2.0*tb12u, &DBa[12], -2.0*g12 );
            /* derivatives of tH_v */
              pkn_MVectorLinCombd ( 12, 3, &DtHv[ka3],
                      &DGa[18], tb22, &DBa[6], g11v, &DGa[0], tb22v,
                      &DBa[24], g11, &DBa[18], g22, &DGa[6], tb11v,
                      &DBa[0], g22v, &DGa[24], tb11, &DGa[21], -2.0*tb12,
                      &DBa[3], -2.0*g12v, &DGa[3], -2.0*tb12v, &DBa[21], -2.0*g12 );
            /* derivatives of N */
              MultVector3d ( -5.5*N/detG, (vector3d*)&DdetG[ka3], (vector3d*)&DN[ka3] );
            /* derivatives of M */
              memcpy ( &DM[3*ka3], &DGa[6], sizeof(vector3d) );
              SetVector3d ( (vector3d*)&DM[3*ka3+3], -DGa[3], -DGa[4], -DGa[5] );
              memcpy ( &DM[3*ka3+6], &DGstar[3*9*ka], sizeof(vector3d) );
            /* derivatives of L */
              pkn_MVectorLinCombd ( 4, 3, &DL[2*ka3],
                      &DdetG[ka3], tHu, &DtHu[ka3], detG,
                      &DdetGu[ka3], -1.5*tH, &DtH[ka3], -1.5*detGu );
              pkn_MVectorLinCombd ( 4, 3, &DL[2*ka3+3],
                      &DdetG[ka3], tHv, &DtHv[ka3], detG,
                      &DdetGv[ka3], -1.5*tH, &DtH[ka3], -1.5*detGv );
          /* compute derivatives of f */
              pkn_MVectorLinCombd ( 2, 3, &DLMLT[0], &DL[2*ka3], MLT[0], &DL[2*ka3+3], MLT[1] );
              pkn_MVectorLinCombd ( 2, 3, &DMLT[0], &DM[3*ka3], L[0], &DM[3*ka3+3], L[1] );
              pkn_MVectorLinCombd ( 2, 3, &DMLT[3], &DM[3*ka3+3], L[0], &DM[3*ka3+6], L[1] );
              pkn_MVectorLinCombd ( 2, 3, LDMLT, &DMLT[0], L[0], &DMLT[3], L[1] );
              pkn_MVectorLinCombd ( 3, 3, &Df[ka3],
                                    DLMLT, 2.0*N, LDMLT, N, &DN[ka3], LMLT );
              if ( tC > 0.0 ) {
          /* compute derivatives of g */
            /* derivative of the Frobenius norm of Laplacian gradient */
                Na = &Nitab[((ka*nkn + i)*nkn + j)*9];
                aa = Na[5] + Na[7];
                bb = Na[6] + Na[8];
                Dg[ka3]   = 2.0*(aa*pder[9].x + bb*pder[10].x);
                Dg[ka3+1] = 2.0*(aa*pder[9].y + bb*pder[10].y);
                Dg[ka3+2] = 2.0*(aa*pder[9].z + bb*pder[10].z);
            /* derivative of the projection on the surface normal */
                MultVector3d ( -1.0/(detG*detG), (vector3d*)&DdetG[ka3],
                               (vector3d*)&DA[ka3] );
                AddVector3Md ( (vector3d*)&Dg[ka3], (vector3d*)&DA[ka3], -BB,
                               (vector3d*)&Dg[ka3] );
                MultVector3d ( Bstar[9], (vector3d*)&DBa[9*3], (vector3d*)&DBB[ka3] );
                AddVector3Md ( (vector3d*)&DBB[ka3], (vector3d*)&DBa[10*3], Bstar[10],
                               (vector3d*)&DBB[ka3] );
                AddVector3Md ( (vector3d*)&Dg[ka3], (vector3d*)&DBB[ka3], -2.0/detG,
                               (vector3d*)&Dg[ka3] );
              }
            }
        /* integrate the functional derivatives */
      if ( tC > 0.0 )
        pkn_MatrixLinCombd ( 1, 3*16, 0, Df, 0.25, 0, Dg, tC, 0, Df );
      pkn_AddMatrixMd ( 1, 3*16, 0, sumD1, 0, Df, qcoeff[j], 0, sumD1 );

        /* compute the second order functional derivatives */
      for ( ia = 0, iap = isq;  ia < 4;  ia++, iap++ )
        if ( iap >= ip0 && iap < ip1 )
          for ( ja = 0, jap = jsq;  ja < 4;  ja++, jap++ )
            if ( jap >= jp0 && jap < jp1 ) {
              ka = 4*ia + ja;
              ka3 = 3*ka;
              /*kap = iap*pitch + jap;*/
              Na = &Nitab[((ka*nkn + i)*nkn + j)*9];
              DGa = &DGstar[9*ka3];
              DBa = &DBstar[11*ka3];
              for ( ib = 0, ibp = isq;  ib < 4;  ib++, ibp++ )
                if ( ibp >= ip0 && ibp < ip1 )
                  for ( jb = 0, jbp = jsq;  jb < 4;  jb++, jbp++ ) {
                    kb = 4*ib + jb;
                    if ( jbp >= jp0 && jbp < jp1 && kb <= ka ) {
                      kb3 = 3*kb;
                      /*kbp = ibp*pitch + jbp;*/
                      Nb = &Nitab[((kb*nkn + i)*nkn + j)*9];
                      DGb = &DGstar[9*kb3];
                      DBb = &DBstar[11*kb3];
                      gstar = &DDGstar[9*pkn_SymMatIndex(ka,kb)];
                      bstar = &DDBstar[3*11*pkn_SymMatIndex(ka,kb)];

                      DDdetG[0] = DDdetG[4] = DDdetG[8] =
                        g22*gstar[0] + g11*gstar[2] - 2.0*g12*gstar[1];
                      DDdetG[0] +=
                        DGb[0]*DGa[6] + DGb[6]*DGa[0] - 2.0*DGb[3]*DGa[3];
                      DDdetG[1] =
                        DGb[0]*DGa[7] + DGb[6]*DGa[1] - 2.0*DGb[3]*DGa[4];
                      DDdetG[2] =
                        DGb[0]*DGa[8] + DGb[6]*DGa[2] - 2.0*DGb[3]*DGa[5];
                      DDdetG[3] =
                        DGb[1]*DGa[6] + DGb[7]*DGa[0] - 2.0*DGb[4]*DGa[3];
                      DDdetG[4] +=
                        DGb[1]*DGa[7] + DGb[7]*DGa[1] - 2.0*DGb[4]*DGa[4];
                      DDdetG[5] =
                        DGb[1]*DGa[8] + DGb[7]*DGa[2] - 2.0*DGb[4]*DGa[5];
                      DDdetG[6] =
                        DGb[2]*DGa[6] + DGb[8]*DGa[0] - 2.0*DGb[5]*DGa[3];
                      DDdetG[7] =
                        DGb[2]*DGa[7] + DGb[8]*DGa[1] - 2.0*DGb[5]*DGa[4];
                      DDdetG[8] +=
                        DGb[2]*DGa[8] + DGb[8]*DGa[2] - 2.0*DGb[5]*DGa[5];

                      DDdetGu[0] = DDdetGu[4] = DDdetGu[8] =
                        g22*gstar[3] + g11u*gstar[2] + g22u*gstar[0] +
                        g11*gstar[5] - 2.0*(g12*gstar[4] + g12u*gstar[1]);
                      DDdetGu[0] +=
                        DGb[9]*DGa[6] + DGb[6]*DGa[9] + DGb[0]*DGa[15] +
                        DGb[15]*DGa[0] - 2.0*(DGb[3]*DGa[12] + DGb[12]*DGa[3]);
                      DDdetGu[1] =
                        DGb[9]*DGa[7] + DGb[6]*DGa[10] + DGb[0]*DGa[16] +
                        DGb[15]*DGa[1] - 2.0*(DGb[3]*DGa[13] + DGb[12]*DGa[4]);
                      DDdetGu[2] =
                        DGb[9]*DGa[8] + DGb[6]*DGa[11] + DGb[0]*DGa[17] +
                        DGb[15]*DGa[2] - 2.0*(DGb[3]*DGa[14] + DGb[12]*DGa[5]);
                      DDdetGu[3] =
                        DGb[10]*DGa[6] + DGb[7]*DGa[9] + DGb[1]*DGa[15] +
                        DGb[16]*DGa[0] - 2.0*(DGb[4]*DGa[12] + DGb[13]*DGa[3]);
                      DDdetGu[4] +=
                        DGb[10]*DGa[7] + DGb[7]*DGa[10] + DGb[1]*DGa[16] +
                        DGb[16]*DGa[1] - 2.0*(DGb[4]*DGa[13] + DGb[13]*DGa[4]);
                      DDdetGu[5] =
                        DGb[10]*DGa[8] + DGb[7]*DGa[11] + DGb[1]*DGa[17] +
                        DGb[16]*DGa[2] - 2.0*(DGb[4]*DGa[14] + DGb[13]*DGa[5]);
                      DDdetGu[6] =
                        DGb[11]*DGa[6] + DGb[8]*DGa[9] + DGb[2]*DGa[15] +
                        DGb[17]*DGa[0] - 2.0*(DGb[5]*DGa[12] + DGb[14]*DGa[3]);
                      DDdetGu[7] =
                        DGb[11]*DGa[7] + DGb[8]*DGa[10] + DGb[2]*DGa[16] +
                        DGb[17]*DGa[1] - 2.0*(DGb[5]*DGa[13] + DGb[14]*DGa[4]);
                      DDdetGu[8] +=
                        DGb[11]*DGa[8] + DGb[8]*DGa[11] + DGb[2]*DGa[17] +
                        DGb[17]*DGa[2] - 2.0*(DGb[5]*DGa[14] + DGb[14]*DGa[5]);

                      DDdetGv[0] = DDdetGv[4] = DDdetGv[8] =
                        g22*gstar[6] + g11v*gstar[2] + g22v*gstar[0] +
                        g11*gstar[8] - 2.0*(g12*gstar[7] + g12v*gstar[1]);
                      DDdetGv[0] +=
                        DGb[18]*DGa[6] + DGb[6]*DGa[18] + DGb[0]*DGa[24] +
                        DGb[24]*DGa[0] - 2.0*(DGb[3]*DGa[21] + DGb[21]*DGa[3]);
                      DDdetGv[1] =
                        DGb[18]*DGa[7] + DGb[6]*DGa[19] + DGb[0]*DGa[25] +
                        DGb[24]*DGa[1] - 2.0*(DGb[3]*DGa[22] + DGb[21]*DGa[4]);
                      DDdetGv[2] =
                        DGb[18]*DGa[8] + DGb[6]*DGa[20] + DGb[0]*DGa[26] +
                        DGb[24]*DGa[2] - 2.0*(DGb[3]*DGa[23] + DGb[21]*DGa[5]);
                      DDdetGv[3] =
                        DGb[19]*DGa[6] + DGb[7]*DGa[18] + DGb[1]*DGa[24] +
                        DGb[25]*DGa[0] - 2.0*(DGb[4]*DGa[21] + DGb[22]*DGa[3]);
                      DDdetGv[4] +=
                        DGb[19]*DGa[7] + DGb[7]*DGa[19] + DGb[1]*DGa[25] +
                        DGb[25]*DGa[1] - 2.0*(DGb[4]*DGa[22] + DGb[22]*DGa[4]);
                      DDdetGv[5] =
                        DGb[19]*DGa[8] + DGb[7]*DGa[20] + DGb[1]*DGa[26] +
                        DGb[25]*DGa[2] - 2.0*(DGb[4]*DGa[23] + DGb[22]*DGa[5]);
                      DDdetGv[6] =
                        DGb[20]*DGa[6] + DGb[8]*DGa[18] + DGb[2]*DGa[24] +
                        DGb[26]*DGa[0] - 2.0*(DGb[5]*DGa[21] + DGb[23]*DGa[3]);
                      DDdetGv[7] =
                        DGb[20]*DGa[7] + DGb[8]*DGa[19] + DGb[2]*DGa[25] +
                        DGb[26]*DGa[1] - 2.0*(DGb[5]*DGa[22] + DGb[23]*DGa[4]);
                      DDdetGv[8] +=
                        DGb[20]*DGa[8] + DGb[8]*DGa[20] + DGb[2]*DGa[26] +
                        DGb[26]*DGa[2] - 2.0*(DGb[5]*DGa[23] + DGb[23]*DGa[5]);

                      DDtH[0] = DDtH[4] = DDtH[8] =
                        tb22*gstar[0] + tb11*gstar[2] - 2.0*tb12*gstar[1];
                      a = g11*bstar[8] + g22*bstar[2] - 2.0*g12*bstar[5]; 
                      DDtH[1] = a;   DDtH[3] = -a;
                      a = g11*bstar[7] + g22*bstar[1] - 2.0*g12*bstar[4];
                      DDtH[2] = -a;  DDtH[6] = a;
                      a = g11*bstar[6] + g22*bstar[0] - 2.0*g12*bstar[3];
                      DDtH[5] = a;   DDtH[7] = -a;
                      DDtH[0] +=
                        DGb[0]*DBa[6] + DGb[6]*DBa[0] + DBb[0]*DGa[6] + DBb[6]*DGa[0] -
                        2.0*(DGb[3]*DBa[3] + DBb[3]*DGa[3]);
                      DDtH[1] +=
                        DGb[0]*DBa[7] + DGb[6]*DBa[1] + DBb[0]*DGa[7] + DBb[6]*DGa[1] -
                        2.0*(DGb[3]*DBa[4] + DBb[3]*DGa[4]);
                      DDtH[2] +=
                        DGb[0]*DBa[8] + DGb[6]*DBa[2] + DBb[0]*DGa[8] + DBb[6]*DGa[2] -
                        2.0*(DGb[3]*DBa[5] + DBb[3]*DGa[5]);
                      DDtH[3] +=
                        DGb[1]*DBa[6] + DGb[7]*DBa[0] + DBb[1]*DGa[6] + DBb[7]*DGa[0] -
                        2.0*(DGb[4]*DBa[3] + DBb[4]*DGa[3]);
                      DDtH[4] +=
                        DGb[1]*DBa[7] + DGb[7]*DBa[1] + DBb[1]*DGa[7] + DBb[7]*DGa[1] -
                        2.0*(DGb[4]*DBa[4] + DBb[4]*DGa[4]);
                      DDtH[5] +=
                        DGb[1]*DBa[8] + DGb[7]*DBa[2] + DBb[1]*DGa[8] + DBb[7]*DGa[2] -
                        2.0*(DGb[4]*DBa[5] + DBb[4]*DGa[5]);
                      DDtH[6] +=
                        DGb[2]*DBa[6] + DGb[8]*DBa[0] + DBb[2]*DGa[6] + DBb[8]*DGa[0] -
                        2.0*(DGb[5]*DBa[3] + DBb[5]*DGa[3]);
                      DDtH[7] +=
                        DGb[2]*DBa[7] + DGb[8]*DBa[1] + DBb[2]*DGa[7] + DBb[8]*DGa[1] -
                        2.0*(DGb[5]*DBa[4] + DBb[5]*DGa[4]);
                      DDtH[8] +=
                        DGb[2]*DBa[8] + DGb[8]*DBa[2] + DBb[2]*DGa[8] + DBb[8]*DGa[2] -
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
                      DDtHu[0] +=
                        DGb[9]*DBa[6] + DGb[0]*DBa[15] + DGb[6]*DBa[9] + DGb[15]*DBa[0] +
                        DBb[9]*DGa[6] + DBb[0]*DGa[15] + DBb[6]*DGa[9] + DBb[15]*DGa[0] -
                        2.0*(DGb[12]*DBa[3] + DGb[3]*DBa[12] +
                             DBb[12]*DGa[3] + DBb[3]*DGa[12]);
                      DDtHu[1] +=
                        DGb[9]*DBa[7] + DGb[0]*DBa[16] + DGb[6]*DBa[10] + DGb[15]*DBa[1] +
                        DBb[9]*DGa[7] + DBb[0]*DGa[16] + DBb[6]*DGa[10] + DBb[15]*DGa[1] -
                        2.0*(DGb[12]*DBa[4] + DGb[3]*DBa[13] +
                             DBb[12]*DGa[4] + DBb[3]*DGa[13]);
                      DDtHu[2] +=
                        DGb[9]*DBa[8] + DGb[0]*DBa[17] + DGb[6]*DBa[11] + DGb[15]*DBa[2] +
                        DBb[9]*DGa[8] + DBb[0]*DGa[17] + DBb[6]*DGa[11] + DBb[15]*DGa[2] -
                        2.0*(DGb[12]*DBa[5] + DGb[3]*DBa[14] +
                             DBb[12]*DGa[5] + DBb[3]*DGa[14]);
                      DDtHu[3] +=
                        DGb[10]*DBa[6] + DGb[1]*DBa[15] + DGb[7]*DBa[9] + DGb[16]*DBa[0] +
                        DBb[10]*DGa[6] + DBb[1]*DGa[15] + DBb[7]*DGa[9] + DBb[16]*DGa[0] -
                        2.0*(DGb[13]*DBa[3] + DGb[4]*DBa[12] +
                             DBb[13]*DGa[3] + DBb[4]*DGa[12]);
                      DDtHu[4] +=
                        DGb[10]*DBa[7] + DGb[1]*DBa[16] + DGb[7]*DBa[10] + DGb[16]*DBa[1] +
                        DBb[10]*DGa[7] + DBb[1]*DGa[16] + DBb[7]*DGa[10] + DBb[16]*DGa[1] -
                        2.0*(DGb[13]*DBa[4] + DGb[4]*DBa[13] +
                             DBb[13]*DGa[4] + DBb[4]*DGa[13]);
                      DDtHu[5] +=
                        DGb[10]*DBa[8] + DGb[1]*DBa[17] + DGb[7]*DBa[11] + DGb[16]*DBa[2] +
                        DBb[10]*DGa[8] + DBb[1]*DGa[17] + DBb[7]*DGa[11] + DBb[16]*DGa[2] -
                        2.0*(DGb[13]*DBa[5] + DGb[4]*DBa[14] +
                             DBb[13]*DGa[5] + DBb[4]*DGa[14]);
                      DDtHu[6] +=
                        DGb[11]*DBa[6] + DGb[2]*DBa[15] + DGb[8]*DBa[9] + DGb[17]*DBa[0] +
                        DBb[11]*DGa[6] + DBb[2]*DGa[15] + DBb[8]*DGa[9] + DBb[17]*DGa[0] -
                        2.0*(DGb[14]*DBa[3] + DGb[5]*DBa[12] +
                             DBb[14]*DGa[3] + DBb[5]*DGa[12]);
                      DDtHu[7] +=
                        DGb[11]*DBa[7] + DGb[2]*DBa[16] + DGb[8]*DBa[10] + DGb[17]*DBa[1] +
                        DBb[11]*DGa[7] + DBb[2]*DGa[16] + DBb[8]*DGa[10] + DBb[17]*DGa[1] -
                        2.0*(DGb[14]*DBa[4] + DGb[5]*DBa[13] +
                             DBb[14]*DGa[4] + DBb[5]*DGa[13]);
                      DDtHu[8] +=
                        DGb[11]*DBa[8] + DGb[2]*DBa[17] + DGb[8]*DBa[11] + DGb[17]*DBa[2] +
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
                      DDtHv[0] +=
                        DGb[18]*DBa[6] + DGb[0]*DBa[24] + DGb[6]*DBa[18] + DGb[24]*DBa[0] +
                        DBb[18]*DGa[6] + DBb[0]*DGa[24] + DBb[6]*DGa[18] + DBb[24]*DGa[0] -
                        2.0*(DGb[21]*DBa[3] + DGb[3]*DBa[21] +
                             DBb[21]*DGa[3] + DBb[3]*DGa[21]);
                      DDtHv[1] +=
                        DGb[18]*DBa[7] + DGb[0]*DBa[25] + DGb[6]*DBa[19] + DGb[24]*DBa[1] +
                        DBb[18]*DGa[7] + DBb[0]*DGa[25] + DBb[6]*DGa[19] + DBb[24]*DGa[1] -
                        2.0*(DGb[21]*DBa[4] + DGb[3]*DBa[22] +
                             DBb[21]*DGa[4] + DBb[3]*DGa[22]);
                      DDtHv[2] +=
                        DGb[18]*DBa[8] + DGb[0]*DBa[26] + DGb[6]*DBa[20] + DGb[24]*DBa[2] +
                        DBb[18]*DGa[8] + DBb[0]*DGa[26] + DBb[6]*DGa[20] + DBb[24]*DGa[2] -
                        2.0*(DGb[21]*DBa[5] + DGb[3]*DBa[23] +
                             DBb[21]*DGa[5] + DBb[3]*DGa[23]);
                      DDtHv[3] +=
                        DGb[19]*DBa[6] + DGb[1]*DBa[24] + DGb[7]*DBa[18] + DGb[25]*DBa[0] +
                        DBb[19]*DGa[6] + DBb[1]*DGa[24] + DBb[7]*DGa[18] + DBb[25]*DGa[0] -
                        2.0*(DGb[22]*DBa[3] + DGb[4]*DBa[21] +
                             DBb[22]*DGa[3] + DBb[4]*DGa[21]);
                      DDtHv[4] +=
                        DGb[19]*DBa[7] + DGb[1]*DBa[25] + DGb[7]*DBa[19] + DGb[25]*DBa[1] +
                        DBb[19]*DGa[7] + DBb[1]*DGa[25] + DBb[7]*DGa[19] + DBb[25]*DGa[1] -
                        2.0*(DGb[22]*DBa[4] + DGb[4]*DBa[22] +
                             DBb[22]*DGa[4] + DBb[4]*DGa[22]);
                      DDtHv[5] +=
                        DGb[19]*DBa[8] + DGb[1]*DBa[26] + DGb[7]*DBa[20] + DGb[25]*DBa[2] +
                        DBb[19]*DGa[8] + DBb[1]*DGa[26] + DBb[7]*DGa[20] + DBb[25]*DGa[2] -
                        2.0*(DGb[22]*DBa[5] + DGb[4]*DBa[23] +
                             DBb[22]*DGa[5] + DBb[4]*DGa[23]);
                      DDtHv[6] +=
                        DGb[20]*DBa[6] + DGb[2]*DBa[24] + DGb[8]*DBa[18] + DGb[26]*DBa[0] +
                        DBb[20]*DGa[6] + DBb[2]*DGa[24] + DBb[8]*DGa[18] + DBb[26]*DGa[0] -
                        2.0*(DGb[23]*DBa[3] + DGb[5]*DBa[21] +
                             DBb[23]*DGa[3] + DBb[5]*DGa[21]);
                      DDtHv[7] +=
                        DGb[20]*DBa[7] + DGb[2]*DBa[25] + DGb[8]*DBa[19] + DGb[26]*DBa[1] +
                        DBb[20]*DGa[7] + DBb[2]*DGa[25] + DBb[8]*DGa[19] + DBb[26]*DGa[1] -
                        2.0*(DGb[23]*DBa[4] + DGb[5]*DBa[22] +
                             DBb[23]*DGa[4] + DBb[5]*DGa[22]);
                      DDtHv[8] +=
                        DGb[20]*DBa[8] + DGb[2]*DBa[26] + DGb[8]*DBa[20] + DGb[26]*DBa[2] +
                        DBb[20]*DGa[8] + DBb[2]*DGa[26] + DBb[8]*DGa[20] + DBb[26]*DGa[2] -
                        2.0*(DGb[23]*DBa[5] + DGb[5]*DBa[23] +
                             DBb[23]*DGa[5] + DBb[5]*DGa[23]);

          /* the second order derivatives of N */
                      pkn_MultMatrixd ( 3, 1, 1, &DdetG[kb3],
                                        3, 0, &DdetG[ka3], 3, DDN );
                      pkn_MatrixLinCombd ( 1, 9, 0, DDN, 13.0,
                                           0, DDdetG, -2.0*detG, 0, DDN );
                      pkn_MultMatrixNumd ( 1, 9, 0, DDN,
                                           11.0*N/(4.0*detG*detG), 0, DDN );
          /* the second order derivatives of M */
            /* the 3 x 3 matrices of the second order derivatives of the */
            /* three coefficients of M are diagonal, with all diagonal */
            /* coefficients equal, so only one number is needed to represent */
            /* each of them */
                      DDM[0] = gstar[2];
                      DDM[1] = -gstar[1];
                      DDM[2] = gstar[0];
          /* the second order derivatives of L */
                      pkn_MatrixLinCombd ( 1, 9, 0, DDdetG, tHu,
                                           0, DDtHu, detG, 0, &DDL[0] );
                      pkn_MultMatrixAddd ( 3, 1, 1, &DdetG[kb3],
                                           3, 0, &DtHu[ka3], 3, &DDL[0] );
                      pkn_MultMatrixAddd ( 3, 1, 1, &DtHu[kb3],
                                           3, 0, &DdetG[ka3], 3, &DDL[0] );
                      pkn_MatrixLinCombd ( 1, 9, 0, DDdetGu, tH,
                                           0, DDtH, detGu, 0, DDaux );
                      pkn_MultMatrixAddd ( 3, 1, 1, &DdetGu[kb3],
                                           3, 0, &DtH[ka3], 3, DDaux );
                      pkn_MultMatrixAddd ( 3, 1, 1, &DtH[kb3],
                                           3, 0, &DdetGu[ka3], 3, DDaux );
                      pkn_AddMatrixMd ( 1, 9, 0, &DDL[0], 0, DDaux, -1.5,
                                        0, &DDL[0] );

                      pkn_MatrixLinCombd ( 1, 9, 0, DDdetG, tHv,
                                           0, DDtHv, detG, 0, &DDL[9] );
                      pkn_MultMatrixAddd ( 3, 1, 1, &DdetG[kb3],
                                           3, 0, &DtHv[ka3], 3, &DDL[9] );
                      pkn_MultMatrixAddd ( 3, 1, 1, &DtHv[kb3],
                                           3, 0, &DdetG[ka3], 3, &DDL[9] );
                      pkn_MatrixLinCombd ( 1, 9, 0, DDdetGv, tH,
                                           0, DDtH, detGv, 0, DDaux );
                      pkn_MultMatrixAddd ( 3, 1, 1, &DdetGv[kb3],
                                           3, 0, &DtH[ka3], 3, DDaux );
                      pkn_MultMatrixAddd ( 3, 1, 1, &DtH[kb3],
                                           3, 0, &DdetGv[ka3], 3, DDaux );
                      pkn_AddMatrixMd ( 1, 9, 0, &DDL[9], 0, DDaux, -1.5,
                                        0, &DDL[9] );

          /* the second order derivatives of f */
                      for ( l = n = 0;  l < 3;  l++ )
                        for ( m = 0;  m < 3;  m++, n++ ) {
                          if ( l == m )
                            a = L[0]*(DDM[0]*L[0]+DDM[1]*L[1]) +
                                L[1]*(DDM[1]*L[0]+DDM[2]*L[1]);
                          else
                            a = 0.0;
                          DDf[3*m+l] =
                            (2.0*(DDL[n]*MLT[0]+DDL[9+n]*MLT[1] +
                                  DL[6*kb+l]*(DM[9*ka+m]*L[0]+DM[9*ka+3+m]*L[1]) +
                                  DL[6*kb+3+l]*(DM[9*ka+3+m]*L[0]+DM[9*ka+6+m]*L[1]) +
                                  DL[6*ka+m]*(DM[9*kb+l]*L[0]+DM[9*kb+3+l]*L[1]) +
                                  DL[6*ka+3+m]*(DM[9*kb+3+l]*L[0]+DM[9*kb+6+l]*L[1]) +
                                  DL[6*kb+l]*(M[0]*DL[6*ka+m]+M[1]*DL[6*ka+3+m]) +
                                  DL[6*kb+3+l]*(M[1]*DL[6*ka+m]+M[2]*DL[6*ka+3+m])) +
                             a)*N +
                            (2.0*(DL[6*kb+l]*MLT[0]+DL[6*kb+3+l]*MLT[1]) +
                             L[0]*(DM[9*kb+l]*L[0]+DM[9*kb+3+l]*L[1])+
                             L[1]*(DM[9*kb+3+l]*L[0]+DM[9*kb+6+l]*L[1]))*DN[ka3+m] +
                            (2.0*(DL[6*ka+m]*MLT[0]+DL[6*ka+3+m]*MLT[1])+
                             L[0]*(DM[9*ka+m]*L[0]+DM[9*ka+3+m]*L[1])+
                             L[1]*(DM[9*ka+3+m]*L[0]+DM[9*ka+6+m]*L[1]))*DN[kb3+l] +
                            LMLT*DDN[n];
                        }
          /* compute the second order derivatives of g */
            /* Frobenius norm of the Laplacian gradient */
                      if ( tC > 0.0 ) {
                        memset ( DDg, 0, 9*sizeof(double) );
                        aa = Na[5] + Na[7];
                        bb = Na[6] + Na[8];
                        cc = Nb[5] + Nb[7];
                        dd = Nb[6] + Nb[8];
                        DDg[0] = DDg[4] = DDg[8] = 2.0*(aa*cc + bb*dd);
            /* projection on the surface normal */
                        pkn_MultMatrixd ( 3, 1, 1, &DdetG[kb3],
                                          3, 3, &DdetG[ka3], 3, DDA );
                        pkn_MatrixLinCombd ( 1, 9, 0, DDA, 2.0/(detG*detG*detG),
                                             0, DDdetG, -1.0/(detG*detG),
                                             0, DDA );
                        pkn_AddMatrixMd ( 1, 9, 0, DDg, 0, DDA, -BB, 0, DDg );

                        pkn_MultMatrixd ( 3, 1, 1, (double*)&DBB[kb3],
                                          3, 3, (double*)&DA[ka3], 3, DDaux );
                        pkn_AddMatrixMd ( 1, 9, 0, DDg, 0, DDaux, -2.0, 0, DDg );
                        pkn_MultMatrixd ( 3, 1, 1, (double*)&DA[kb3],
                                          3, 3, (double*)&DBB[ka3], 3, DDaux );
                        pkn_AddMatrixMd ( 1, 9, 0, DDg, 0, DDaux, -2.0, 0, DDg );

                        pkn_MatrixLinCombd ( 1, 3, 0, &bstar[9*3], Bstar[9],
                                             0, &bstar[10*3], Bstar[10],
                                             0, aux );
                        pkn_MultMatrixd ( 3, 1, 1, &DBb[9*3],
                                          3, 3, &DBa[9*3], 3, DDaux );
                        pkn_MultMatrixAddd ( 3, 1, 1, &DBb[10*3],
                                             3, 3, &DBa[10*3], 3, DDaux );
                        DDaux[1] += aux[2];  DDaux[3] -= aux[2];
                        DDaux[2] -= aux[1];  DDaux[6] += aux[1];
                        DDaux[5] += aux[0];  DDaux[7] -= aux[0];
                        pkn_AddMatrixMd ( 1, 9, 0, DDg, 0, DDaux, -2.0/detG,
                                          0, DDg );
              /* for historic reasons this is the transposition of the */
              /* necessary matrix, to be fixed later, now just transpose it */
                        pkv_TransposeMatrixd ( 3, 3, 3, DDg, 3, DDaux );
                        memcpy ( DDg, DDaux, 9*sizeof(double) );
            /* integrate */
                        pkn_MatrixLinCombd ( 1, 9, 0, DDf, 0.25, 0, DDg, tC, 0, DDf );
                      }
                      ddi = 9*pkn_SymMatIndex ( ka, kb );
                      pkn_AddMatrixMd ( 1, 9, 0, &sumDD1[ddi],
                                        0, DDf, qcoeff[j], 0, &sumDD1[ddi] );
                    }
                  }
            }
    }
    sum0 += sum1*qcoeff[i];
    pkn_AddMatrixMd ( 1, 3*16, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
    pkn_AddMatrixMd ( 1, 9*8*17, 0, htab, 0, sumDD1, qcoeff[i], 0, htab );
  }
  ftab[0] = sum0;
} /*g2bl_UFuncGradHessianSQd*/

