
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2012                            */
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
#include "g2mblmlprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
void _g2mbl_SCompDGStard ( int nkn2, int nf, double *Nitab, int knot, int k, int l,
                           const vector3d *pder, int *cpind, int *nncpi,
                           vector3d *mvcpn, double *DGstar )
{
  int      fi, p, q;
  double   *Dgstar, *Ni, pun, pvn, puun, puvn, pvvn;
  vector3d *mvi;

  for ( fi = 0; fi < nf; fi++ )
    if ( nncpi[cpind[fi]] >= 0 ) {
      Dgstar = &DGstar[fi*9];
      if ( k == 4 || fi == 0 )
        Ni = &Nitab[(fi*nkn2+knot)*9];
      else {
        p = (fi-1) / 6;
        p = (l - p + k) % k;
        q = (fi-1) % 6;
        Ni = &Nitab[((6*p+q+1)*nkn2+knot)*9];
      }
      mvi = &mvcpn[cpind[fi]];
      pun = DotProduct3d ( &pu, mvi );
      pvn = DotProduct3d ( &pv, mvi );
      puun = DotProduct3d ( &puu, mvi );
      puvn = DotProduct3d ( &puv, mvi );
      pvvn = DotProduct3d ( &pvv, mvi );

      dg11 = 2.0*Ni10*pun;
      dg12 = Ni10*pvn + Ni01*pun;
      dg22 = 2.0*Ni01*pvn;

      dg11u = 2.0*(Ni20*pun + Ni10*puun);
      dg12u = Ni20*pvn + Ni01*puun + Ni10*puvn + Ni11*pun;
      dg22u = 2.0*(Ni11*pvn + Ni01*puvn);

      dg11v = 2.0*(Ni11*pun + Ni10*puvn);
      dg12v = Ni11*pvn + Ni01*puvn + Ni10*pvvn + Ni02*pun;
      dg22v = 2.0*(Ni02*pvn + Ni01*pvvn);
    }
} /*_g2mbl_SCompDGStard*/

void _g2mbl_SCompDBStard ( int nkn2, int nf, double *Nitab, int knot, int k, int l,
                           const vector3d *pder, int *cpind, int *nncpi,
                           vector3d *mvcpn, double *DBstar )
{
  int      fi, p, q;
  double   *Ni, *Dbstar;
  double   pupvn, pupuun, pvpuun, pupuvn, pvpuvn, pupvvn, pvpvvn,
           pupuuun, pvpuuun, pupuuvn, pvpuuvn, pupuvvn, pvpuvvn,
           pupvvvn, pvpvvvn, puupuvn, puupvvn, puvpvvn;
  double   w, y, z;
  vector3d *mvi, mvipu, mvipv, mvipuu, mvipuv;

  for ( fi = 0; fi < nf; fi++ )
    if ( nncpi[cpind[fi]] >= 0 ) {
      Dbstar = &DBstar[fi*11];
      if ( k == 4 || fi == 0 )
        Ni = &Nitab[(fi*nkn2+knot)*9];
      else {
        p = (fi-1) / 6;
        p = (l - p + k) % k;
        q = (fi-1) % 6;
        Ni = &Nitab[((6*p+q+1)*nkn2+knot)*9];
      }
      mvi = &mvcpn[cpind[fi]];
        /* to be optimized - many common subexpressions! */
      CrossProduct3d ( mvi, &pu, &mvipu );
      CrossProduct3d ( mvi, &pv, &mvipv );
      CrossProduct3d ( mvi, &puu, &mvipuu );
      CrossProduct3d ( mvi, &puv, &mvipuv );
      pupvn = DotProduct3d ( &mvipu, &pv );
      pupuun = DotProduct3d ( &mvipu, &puu );
      pvpuun = DotProduct3d ( &mvipv, &puu );
      pupuvn = DotProduct3d ( &mvipu, &puv );
      pvpuvn = DotProduct3d ( &mvipv, &puv );
      pupvvn = DotProduct3d ( &mvipu, &pvv );
      pvpvvn = DotProduct3d ( &mvipv, &pvv );
      pupuuun = DotProduct3d ( &mvipu, &puuu );
      pvpuuun = DotProduct3d ( &mvipv, &puuu );
      pupuuvn = DotProduct3d ( &mvipu, &puuv );
      pvpuuvn = DotProduct3d ( &mvipv, &puuv );
      pupuvvn = DotProduct3d ( &mvipu, &puvv );
      pvpuvvn = DotProduct3d ( &mvipv, &puvv );
      pupvvvn = DotProduct3d ( &mvipu, &pvvv );
      pvpvvvn = DotProduct3d ( &mvipv, &pvvv );
      puupuvn = DotProduct3d ( &mvipuu, &puv );
      puupvvn = DotProduct3d ( &mvipuu, &pvv );
      puvpvvn = DotProduct3d ( &mvipuv, &pvv );

      dtb11 = Ni20*pupvn;
      dtb12 = Ni11*pupvn;
      dtb22 = Ni02*pupvn;
      dtb11u = Ni30*pupvn;
      dtb12u = Ni21*pupvn;
      dtb22u = Ni12*pupvn;
      dtb22v = Ni03*pupvn;
/*
      Dbstar[9] = dtb11u + dtb22u;
      Dbstar[10] = dtb12u + dtb22v;
*/
      dtb11 -= Ni01*pupuun;
      dtb11u -= Ni11*pupuun;
      dtb11v = -Ni02*pupuun;

      dtb11 += Ni10*pvpuun;
      z = Ni11*pvpuun;
      y = -Ni02*pvpuun;

      dtb12 += Ni10*pvpuvn;
      z -= Ni20*pvpuvn;
      dtb22v -= Ni02*pvpuvn;

      dtb12 -= Ni01*pupuvn;
      dtb11u += Ni20*pupuvn;
      w = -Ni02*pupuvn;

      dtb22 += Ni10*pvpvvn;
      dtb22v += Ni11*pvpvvn;
      y += Ni20*pvpvvn;

      dtb22 -= Ni01*pupvvn;
      dtb11v += Ni20*pupvvn;
      w += Ni11*pupvvn;

      dtb11u += Ni10*pvpuuun - Ni01*pupuuun;

      dtb11u -= Ni10*puupuvn;
      z += Ni01*puupuvn;

      dtb12u += Ni10*pvpuuvn;

      dtb12u -= Ni01*pupuuvn;

      w -= Ni10*puvpvvn;
      dtb22v -= Ni01*puvpvvn;

      dtb22u += Ni10*pvpuvvn;

      dtb22u -= Ni01*pupuvvn;
      dtb12v = w + dtb22u;
      dtb22u -= w;

      dtb11v += z;
      dtb11v += dtb12u;
      dtb12u -= z;

      dtb22v += Ni10*pvpvvvn;

      dtb22v -= Ni01*pupvvvn;

      dtb22u -= Ni01*puupvvn;
      dtb11v -= Ni10*puupvvn;

      dtb22u += y;
/*
      Dbstar[9] += Ni10*(pvpuuun + pvpuvvn);
      Dbstar[9] -= Ni01*(pupuuun + pupuvvn);

      Dbstar[10] += Ni10*(pvpuuvn + pvpvvvn);
      Dbstar[10] -= Ni01*(pupuuvn + pupvvvn);
*/
    }
} /*_g2mbl_SCompDBStard*/

void _g2mbl_SFuncGradSQIntegrandd ( int nkn2, int nf, double *Nitab,
                                    int knot, int k, int l, vector3d pder[11],
                                    int *cpind, int *nncpi,
                                    double Gstar[9], double *DGstar,
                                    double Bstar[11], double *DBstar,
                                    double *F, double *DF )
{
  double detG, detGu, detGv, tH, tHu, tHv;
  double DdetG, DdetGu, DdetGv, DtH, DtHu, DtHv;
  double L[2], M[3], N, DL[2], DM[3], DN;
  double MLT[2], LMLT, DMLT[2], DLMLT, LDMLT;
  int    fn;
  double *DGf, *DBf;

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
  *F = 0.25*LMLT*N;

        /* compute the functional derivatives */
  for ( fn = 0;  fn < nf;  fn++ )
    if ( nncpi[cpind[fn]] >= 0 ) {
      DGf = &DGstar[9*fn];
      DBf = &DBstar[11*fn];
          /* derivative of det G */
      DdetG = g22*DGf[0] + g11*DGf[2] - 2.0*g12*DGf[1];
          /* derivative of (det G)_u */
      DdetGu = g22*DGf[3] + g11u*DGf[2] + g22u*DGf[0] + g11*DGf[5]
               -2.0*(g12u*DGf[1] + g12*DGf[4]);
          /* derivative of (det G)_u */
      DdetGv = g22*DGf[6] + g11v*DGf[2] + g22v*DGf[0] + g11*DGf[8]
               -2.0*(g12v*DGf[1] + g12*DGf[7]);
          /* derivative of tH */
      DtH = tb22*DGf[0] + g11*DBf[2] + g22*DBf[0] + tb11*DGf[2]
            -2.0*(tb12*DGf[1] + g12*DBf[1]);
          /* derivative of tH_u */
      DtHu = tb22*DGf[3] + g11u*DBf[2] + tb22u*DGf[0] + g11*DBf[5] + g22*DBf[3]
             + tb11u*DGf[2] + g22u*DBf[0] + tb11*DGf[5]
             -2.0*(tb12*DGf[4] + g12u*DBf[1] + tb12u*DGf[1] + g12*DBf[4] );
          /* derivative of tH_v */
      DtHv = tb22*DGf[6] + g11v*DBf[2] + tb22v*DGf[0] + g11*DBf[8] + g22*DBf[6]
             + tb11v*DGf[2] + g22v*DBf[0] + tb11*DGf[8]
             -2.0*(tb12*DGf[7] + g12v*DBf[1] + tb12v*DGf[1] + g12*DBf[7] );
          /* derivative of N */
      DN = -5.5*N*DdetG/detG;
          /* derivative of M */
      DM[0] = DGf[2];
      DM[1] = -DGf[1];
      DM[2] = DGf[0];
          /* derivative of L */
      DL[0] = tHu*DdetG + detG*DtHu - 1.5*(tH*DdetGu + detGu*DtH);
      DL[1] = tHv*DdetG + detG*DtHv - 1.5*(tH*DdetGv + detGv*DtH);
          /* derivative of f */
      DLMLT = MLT[0]*DL[0] + MLT[1]*DL[1];
      DMLT[0] = L[0]*DM[0] + L[1]*DM[1];
      DMLT[1] = L[0]*DM[1] + L[1]*DM[2];
      LDMLT = L[0]*DMLT[0] + L[1]*DMLT[1];
      DF[fn] = 0.25*((2.0*DLMLT + LDMLT)*N + LMLT*DN);
    }
} /*_g2mbl_SFuncGradSQIntegrandd*/

void _g2mbl_SCompDDGstard ( int nkn2, int nf, double *Nijtab,
                            int knot, int k, int l,
                            int *cpind, int *nncpi,
                            vector3d *mvcpn, double *DDGstar )
{
  int      fi, fj, fii, fjj, p, q;
  double   *Nij, spr;
  vector3d *mvi, *mvj;

  for ( fi = 0; fi < nf; fi++ )
    if ( nncpi[cpind[fi]] >= 0 ) {
      mvi = &mvcpn[cpind[fi]];
      if ( k != 4 && fi > 0 ) {
        p = (fi-1) / 6;
        p = (l-p+k) % k;
        q = (fi-1) % 6;
        fii = 6*p+q+1;
      }
      else fii = fi;
      for ( fj = 0; fj <= fi; fj++ )
        if ( nncpi[cpind[fj]] >= 0 ) {
          mvj = &mvcpn[cpind[fj]];
          if ( k != 4 && fj > 0 ) {
            p = (fj-1) / 6;
            p = (l-p+k) % k;
            q = (fj-1) % 6;
            fjj = 6*p+q+1;
          }
          else fjj = fj;
          spr = DotProduct3d ( mvi, mvj );
          Nij = &Nijtab[(pkn_SymMatIndex(fii,fjj)*nkn2+knot)*9];
          pkn_MultMatrixNumd ( 1, 9, 0, Nij, spr,
                               0, &DDGstar[9*pkn_SymMatIndex(fi,fj)] );
        }
    }
} /*_g2mbl_SCompDDGstard*/

void _g2mbl_SCompDDBstard ( int nkn2, int nf, double *Mijtab, int knot, int k, int l,
                            vector3d *pder, int *cpind, int *nncpi,
                            vector3d *mvcpn, double *DDBstar )
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
/*#define ddtlg1  DDB[9]*/
/*#define ddtlg2  DDB[10]*/
  int      fi, fj, fii, fjj, p, q;
  double   *Mij, *DDB;
  vector3d *mvi, *mvj, cpr, w;
  double   DDtb301001, DDtb211001, DDtb121001, DDtb031001, DDtb201011,
           DDtb112001, DDtb022001, DDtb201002, DDtb111002, DDtb021101;

  for ( fi = 0; fi < nf; fi++ )
    if ( nncpi[cpind[fi]] >= 0 ) {
      mvi = &mvcpn[cpind[fi]];
      if ( k == 4 || fi == 0 )
        fii = fi;
      else {
        p = (fi-1) / 6;
        p = (l-p+k) % k;
        q = (fi-1) % 6;
        fii = 6*p+q+1;
      }
      DDB = &DDBstar[11*pkn_SymMatIndex(fi,fi)];
      memset ( DDB, 0, 11*sizeof(double) );
      for ( fj = fi+1; fj < nf; fj++ )
        if ( nncpi[cpind[fj]] >= 0 ) {
          mvj = &mvcpn[cpind[fj]];
          if ( k == 4 )
            fjj = fj;
          else {
            p = (fj-1) / 6;
            p = (l-p+k) % k;
            q = (fj-1) % 6;
            fjj = 6*p+q+1;
          }
          DDB = &DDBstar[11*pkn_SymMatIndex(fi,fj)];
          if ( fii < fjj ) {
            Mij = &Mijtab[(pkn_SymMatIndex(fii,fjj-1)*nkn2+knot)*18];
            CrossProduct3d ( mvi, mvj, &cpr );
          }
          else {
            Mij = &Mijtab[(pkn_SymMatIndex(fii-1,fjj)*nkn2+knot)*18];
            CrossProduct3d ( mvj, mvi, &cpr );
          }

            /* DDtb301001 */
          MultVector3d ( Mij1001, &puuu, &w );
          AddVector3Md ( &w, &pu, Mij0130, &w );
          AddVector3Md ( &w, &pv, Mij3010, &w );
          DDtb301001 = DotProduct3d ( &cpr, &w );
            /* DDtb211001 */
          MultVector3d ( Mij1001, &puuv, &w );
          AddVector3Md ( &w, &pu, Mij0121, &w );
          AddVector3Md ( &w, &pv, Mij2110, &w );
          DDtb211001 = DotProduct3d ( &cpr, &w );
            /* DDtb121001 */
          MultVector3d ( Mij1001, &puvv, &w );
          AddVector3Md ( &w, &pu, Mij0112, &w );
          AddVector3Md ( &w, &pv, Mij1210, &w );
          DDtb121001 = DotProduct3d ( &cpr, &w );
            /* DDtb031001 */
          MultVector3d ( Mij1001, &pvvv, &w );
          AddVector3Md ( &w, &pu, Mij0103, &w );
          AddVector3Md ( &w, &pv, Mij0310, &w );
          DDtb031001 = DotProduct3d ( &cpr, &w );
            /* DDtb201011 */
          MultVector3d ( Mij1011, &puu, &w );
          AddVector3Md ( &w, &pu, Mij1120, &w );
          AddVector3Md ( &w, &puv, Mij2010, &w );
          DDtb201011 = DotProduct3d ( &cpr, &w );
            /* DDtb112001 */
          MultVector3d ( Mij2001, &puv, &w );
          AddVector3Md ( &w, &puu, Mij0111, &w );
          AddVector3Md ( &w, &pv, Mij1120, &w );
          DDtb112001 = DotProduct3d ( &cpr, &w );
            /* DDtb022001 */
          MultVector3d ( Mij2001, &pvv, &w );
          AddVector3Md ( &w, &puu, Mij0102, &w );
          AddVector3Md ( &w, &pv, Mij0220, &w );
          DDtb022001 = DotProduct3d ( &cpr, &w );
            /* DDtb201002 */
          MultVector3d ( Mij1002, &puu, &w );
          AddVector3Md ( &w, &pu, Mij0220, &w );
          AddVector3Md ( &w, &pvv, Mij2010, &w );
          DDtb201002 = DotProduct3d ( &cpr, &w );
            /* DDtb111002 */ 
          MultVector3d ( Mij1002, &puv, &w );
          AddVector3Md ( &w, &pu, Mij0211, &w );
          AddVector3Md ( &w, &pvv, Mij1110, &w );
          DDtb111002 = DotProduct3d ( &cpr, &w );
            /* DDtb021101 */
          MultVector3d ( Mij1101, &pvv, &w );
          AddVector3Md ( &w, &puv, Mij0102, &w );
          AddVector3Md ( &w, &pv, Mij0211, &w );
          DDtb021101 = DotProduct3d ( &cpr, &w );
            /* DDtb11 */
          MultVector3d ( Mij1001, &puu, &w );
          AddVector3Md ( &w, &pu, Mij0120, &w );
          AddVector3Md ( &w, &pv, Mij2010, &w );
          ddtb11 = DotProduct3d ( &cpr, &w );
            /* DDtb12 */
          MultVector3d ( Mij1001, &puv, &w );
          AddVector3Md ( &w, &pu, Mij0111, &w );
          AddVector3Md ( &w, &pv, Mij1110, &w );
          ddtb12 = DotProduct3d ( &cpr, &w );
            /* DDtb22 */
          MultVector3d ( Mij1001, &pvv, &w );
          AddVector3Md ( &w, &pu, Mij0102, &w );
          AddVector3Md ( &w, &pv, Mij0210, &w );
          ddtb22 = DotProduct3d ( &cpr, &w );
            /* DDtb11u */
          ddtb11u = DDtb301001 + DDtb201011;
            /* DDtb12u */
          ddtb12u = DDtb211001 + DDtb112001;
            /* DDtb22u */
          ddtb22u = DDtb121001 + DDtb022001 - DDtb111002;
            /* DDtb11v */
          ddtb11v = DDtb211001 - DDtb112001 + DDtb201002;
            /* DDtb12v */
          ddtb12v = DDtb121001 + DDtb111002;
            /* DDtb22v */
          ddtb22v = DDtb031001 + DDtb021101;
            /* DD(h100130+h100112) */
/*          ddtlg1 = DDtb301001 + DDtb121001; */
            /* DD(h100121+h100103) */
/*          ddtlg2 = DDtb211001 + DDtb031001; */
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
/*#undef ddtlg1*/
/*#undef ddtlg2*/
} /*_g2mbl_SCompDDBstard*/

boolean _g2mbl_SFuncGradHessSQIntegrandd ( int nkn2, int nf, double *Nitab,
                                int knot, int k, int l, vector3d *pder,
                                int *cpind, int *nncpi,
                                double *Gstar, double *DGstar, double *DDGstar,
                                double *Bstar, double *DBstar, double *DDBstar,
                                double *F, double *DF, double *DDF )
{
  void   *sp;
  double detG, detGu, detGv, tH, tHu, tHv;
  double *DdetG, *DdetGu, *DdetGv, *DtH, *DtHu, *DtHv;
  double L[2], M[3], N, *DL, *DM, *DN;
  double DDdetG, DDdetGu, DDdetGv, DDtH, DDtHu, DDtHv;
  double DDL[2], DDM[3], DDN, ND;
  double MLT[2], LMLT, DMLT[2], DLMLT, LDMLT;
  int    fi, fj;
  double *DGa, *DGb, *DBa, *DBb, *gstar, *bstar;

  sp = pkv_GetScratchMemTop ();
  DdetG = pkv_GetScratchMemd ( 12*nf );
  if ( !DdetG )
    goto failure;
  DdetGu = &DdetG[nf];
  DdetGv = &DdetGu[nf];
  DtH    = &DdetGv[nf];
  DtHu   = &DtH[nf];
  DtHv   = &DtHu[nf];
  DL     = &DtHv[nf];
  DM     = &DL[2*nf];
  DN     = &DM[3*nf];
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
  *F = 0.25*LMLT*N;

        /* compute the functional derivatives */
  ND = -5.5*N/detG;
  for ( fi = 0;  fi < nf;  fi++ )
    if ( nncpi[cpind[fi]] >= 0 ) {
      DGa = &DGstar[9*fi];
      DBa = &DBstar[11*fi];
          /* derivative of det G */
      DdetG[fi] = g22*DGa[0] + g11*DGa[2] - 2.0*g12*DGa[1];
          /* derivative of (det G)_u */
      DdetGu[fi] = g22*DGa[3] + g11u*DGa[2] + g22u*DGa[0] + g11*DGa[5]
                   -2.0*(g12u*DGa[1] + g12*DGa[4]);
          /* derivative of (det G)_u */
      DdetGv[fi] = g22*DGa[6] + g11v*DGa[2] + g22v*DGa[0] + g11*DGa[8]
                   -2.0*(g12v*DGa[1] + g12*DGa[7]);
          /* derivative of tH */
      DtH[fi] = tb22*DGa[0] + g11*DBa[2] + g22*DBa[0] + tb11*DGa[2]
                -2.0*(tb12*DGa[1] + g12*DBa[1]);
          /* derivative of tH_u */
      DtHu[fi] = tb22*DGa[3] + g11u*DBa[2] + tb22u*DGa[0] + g11*DBa[5] + g22*DBa[3]
                 + tb11u*DGa[2] + g22u*DBa[0] + tb11*DGa[5]
                 -2.0*(tb12*DGa[4] + g12u*DBa[1] + tb12u*DGa[1] + g12*DBa[4] );
          /* derivative of tH_v */
      DtHv[fi] = tb22*DGa[6] + g11v*DBa[2] + tb22v*DGa[0] + g11*DBa[8] + g22*DBa[6]
                 + tb11v*DGa[2] + g22v*DBa[0] + tb11*DGa[8]
                 -2.0*(tb12*DGa[7] + g12v*DBa[1] + tb12v*DGa[1] + g12*DBa[7] );
          /* derivative of N */
      DN[fi] = ND*DdetG[fi];
          /* derivative of M */
      DM[3*fi+0] = DGa[2];
      DM[3*fi+1] = -DGa[1];
      DM[3*fi+2] = DGa[0]; 
          /* derivative of L */
      DL[2*fi+0] = tHu*DdetG[fi] + detG*DtHu[fi] - 1.5*(tH*DdetGu[fi] + detGu*DtH[fi]);
      DL[2*fi+1] = tHv*DdetG[fi] + detG*DtHv[fi] - 1.5*(tH*DdetGv[fi] + detGv*DtH[fi]);
          /* derivative of f */
      DLMLT = MLT[0]*DL[2*fi+0] + MLT[1]*DL[2*fi+1];
      DMLT[0] = L[0]*DM[3*fi+0] + L[1]*DM[3*fi+1];
      DMLT[1] = L[0]*DM[3*fi+1] + L[1]*DM[3*fi+2];
      LDMLT = L[0]*DMLT[0] + L[1]*DMLT[1];
      DF[fi] = 0.25*((2.0*DLMLT + LDMLT)*N + LMLT*DN[fi]);
    }

        /* compute the second order derivatives */
  ND /= -2.0*detG;
  for ( fi = 0;  fi < nf;  fi++ )
    if ( nncpi[cpind[fi]] >= 0 ) {
      DGa = &DGstar[9*fi];
      DBa = &DBstar[11*fi];
      for ( fj = 0;  fj <= fi;  fj++ )
        if ( nncpi[cpind[fj]] >= 0 ) {
          DGb = &DGstar[9*fj];
          DBb = &DBstar[11*fj];
          gstar = &DDGstar[9*pkn_SymMatIndex(fi,fj)];
          bstar = &DDBstar[11*pkn_SymMatIndex(fi,fj)];

            /* derivative of detG */
          DDdetG = g22*gstar[0] + g11*gstar[2] - 2.0*g12*gstar[1]
                   + DGb[2]*DGa[0] + DGb[0]*DGa[2] - 2.0*DGb[1]*DGa[1];
            /* derivative of detG_u */
          DDdetGu = g22*gstar[3] + DGb[2]*DGa[3] + DGb[3]*DGa[2] + g11u*gstar[2]
                    + g22u*gstar[0] + DGb[5]*DGa[0] + DGb[0]*DGa[5] + g11*gstar[5]
                    - 2.0*(g12u*gstar[1] + DGb[4]*DGa[1] + DGb[1]*DGa[4] + g12*gstar[4]);
            /* derivative of detG_v */
          DDdetGv = g22*gstar[6] + DGb[2]*DGa[6] + DGb[6]*DGa[2] + g11v*gstar[2]
                    + g22v*gstar[0] + DGb[8]*DGa[0] + DGb[0]*DGa[8] + g11*gstar[8]
                    - 2.0*(g12v*gstar[1] + DGb[7]*DGa[1] + DGb[1]*DGa[7] + g12*gstar[7]);
            /* derivative of tH */
          DDtH = tb22*gstar[0] + DBb[2]*DGa[0] + bstar[2]*g11 + DGb[0]*DBa[2]
                 + bstar[0]*g22 + DGb[2]*DBa[0] + tb11*gstar[2] + DBb[0]*DGa[2]
                 - 2.0*(tb12*gstar[1] + DBb[1]*DGa[1] + DBa[1]*DGb[1] + bstar[1]*g12);
            /* derivative of tH_u */
          DDtHu = tb22*gstar[3] + DBb[2]*DGa[3] + bstar[2]*g11u + DGb[3]*DBa[2]
                  + tb22u*gstar[0] + DBb[5]*DGa[0] + g11*bstar[5] + DGb[0]*DBa[5]
                  + bstar[3]*g22 + DGb[2]*DBa[3] + tb11u*gstar[2] + DBb[3]*DGa[2]
                  + bstar[0]*g22u + DGb[5]*DBa[0] + tb11*gstar[5] + DBb[0]*DGa[5]
                  - 2.0*(tb12*gstar[4] + DBb[1]*DGa[4] + bstar[1]*g12u + DGb[4]*DBa[1]
                   + tb12u*gstar[1] + DBb[4]*DGa[1] + g12*bstar[4] + DGb[1]*DBa[4]);
            /* derivative of tH_v */
          DDtHv = tb22*gstar[6] + DBb[2]*DGa[6] + bstar[2]*g11v + DGb[6]*DBa[2]
                  + tb22v*gstar[0] + DBb[8]*DGa[0] + g11*bstar[8] + DGb[0]*DBa[8]
                  + bstar[6]*g22 + DGb[2]*DBa[6] + tb11v*gstar[2] + DBb[6]*DGa[2]
                  + bstar[0]*g22v + DGb[8]*DBa[0] + tb11*gstar[8] + DBb[0]*DGa[8]
                  - 2.0*(tb12*gstar[7] + DBb[1]*DGa[7] + bstar[1]*g12v + DGb[7]*DBa[1]
                   + tb12v*gstar[1] + DBb[7]*DGa[1] + g12*bstar[7] + DGb[1]*DBa[7]);
            /* the second order derivative of N */
          DDN = ND*(13.0*DdetG[fi]*DdetG[fj] - 2.0*detG*DDdetG);
            /* the second order derivative of M */
          DDM[0] = gstar[2];
          DDM[1] = -gstar[1];
          DDM[2] = gstar[0];
            /* the second order derivative of L */
          DDL[0] = DDdetG*tHu + DdetG[fi]*DtHu[fj] + DdetG[fj]*DtHu[fi] + detG*DDtHu
                   - 1.5*(DDdetGu*tH + DdetGu[fi]*DtH[fj] + DdetGu[fj]*DtH[fi]
                          + detGu*DDtH);
          DDL[1] = DDdetG*tHv + DdetG[fi]*DtHv[fj] + DdetG[fj]*DtHv[fi] + detG*DDtHv
                   - 1.5*(DDdetGv*tH + DdetGv[fi]*DtH[fj] + DdetGv[fj]*DtH[fi]
                          + detGv*DDtH);
            /* the second order derivatives of f */
          DDF[pkn_SymMatIndex(fi, fj)] =
0.25*(
  (2.0*( (DDL[0]*MLT[0]+DDL[1]*MLT[1]) +
         (DL[2*fi+0]*(DM[3*fj+0]*L[0]+DM[3*fj+1]*L[1])+
          DL[2*fi+1]*(DM[3*fj+1]*L[0]+DM[3*fj+2]*L[1])) +
         (DL[2*fj+0]*(DM[3*fi+0]*L[0]+DM[3*fi+1]*L[1])+
          DL[2*fj+1]*(DM[3*fi+1]*L[0]+DM[3*fi+2]*L[1])) +
         (DL[2*fi+0]*(M[0]*DL[2*fj+0]+M[1]*DL[2*fj+1])+
          DL[2*fi+1]*(M[1]*DL[2*fj+0]+M[2]*DL[2*fj+1]))
       ) +
    (L[0]*(DDM[0]*L[0]+DDM[1]*L[1])+L[1]*(DDM[1]*L[0]+DDM[2]*L[1]))
  )*N +
  (2.0*(DL[2*fi+0]*MLT[0]+DL[2*fi+1]*MLT[1])+
   (L[0]*(DM[3*fi+0]*L[0]+DM[3*fi+1]*L[1])+
    L[1]*(DM[3*fi+1]*L[0]+DM[3*fi+2]*L[1]))
  )*DN[fj] +
  (2.0*(DL[2*fj+0]*MLT[0]+DL[2*fj+1]*MLT[1])+
   (L[0]*(DM[3*fj+0]*L[0]+DM[3*fj+1]*L[1])+
    L[1]*(DM[3*fj+1]*L[0]+DM[3*fj+2]*L[1]))
  )*DN[fi] +
  LMLT*DDN
);
/* DEBUG - to remove */
/*
DDF[pkn_SymMatIndex(fi, fj)] = DDL[1];
*/
        }
    }
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_SFuncGradHessSQIntegrandd*/

/* ///////////////////////////////////////////////////////////////////////// */
void g2mbl_SFuncGradRSQd ( int nkn, double *qcoeff, double *Nitab,
                           int *cpind, int *nncpi, point3d *mvcp, vector3d *mvcpn,
                           double *ftab, double *gtab )
{
  int      nkn2, i, j, knot;
  vector3d pder[11];
  double   Gstar[9], DGstar[9*16], Bstar[11], DBstar[11*16];
  double   sum0, sum1, f;
  double   Df[16], sumD1[16];

  nkn2 = nkn*nkn;
  memset ( Df, 0, 16*sizeof(double) );
  sum0 = 0.0;
  memset ( gtab, 0, 16*sizeof(double) );
  for ( i = knot = 0;  i < nkn;  i++ ) {
    sum1 = 0.0;
    memset ( sumD1, 0, 16*sizeof(double) );
    for ( j = 0;  j < nkn;  j++, knot++ ) {
        /* compute the parameterization derivatives */
      _g2mbl_UCompRDerd ( nkn2, Nitab, knot, cpind, mvcp, pder );
        /* compute the first fundamental form and its derivatives */
      _g2bl_UCompGStard ( pder, Gstar );
      _g2mbl_SCompDGStard ( nkn2, 16, Nitab, knot, 4, 0, pder, cpind, nncpi,
                            mvcpn, DGstar );
        /* compute the second fundamental form and its derivatives */
      _g2bl_UCompBStard ( pder, Bstar );
      _g2mbl_SCompDBStard ( nkn2, 16, Nitab, knot, 4, 0, pder, cpind, nncpi,
                            mvcpn, DBstar );
        /* compute the integrand terms */
      _g2mbl_SFuncGradSQIntegrandd ( nkn2, 16, Nitab, knot, 4, 0, pder,
                        cpind, nncpi, Gstar, DGstar, Bstar, DBstar, &f, Df );
        /* add the quadrature term */
      sum1 += f*qcoeff[j];
      pkn_AddMatrixMd ( 1, 16, 0, sumD1, 0, Df, qcoeff[j], 0, sumD1 );
    }
    sum0 += sum1*qcoeff[i];
    pkn_AddMatrixMd ( 1, 16, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
  }
  *ftab = sum0;
} /*g2mbl_SFuncGradRSQd*/

boolean g2mbl_SFuncGradSSQd ( int nkn, double *qcoeff, int k,
                              double *Nitab, double *Jac,
                              int *cpind, int *nncpi, point3d *mvcp, vector3d *mvcpn,
                              double *ftab, double *gtab )
{
  void     *sp;
  int      nkn2, i, j, l, knot, nf;
  vector3d *pder;
  double   *Gstar, *DGstar, *Bstar, *DBstar;
  double   sum0, sum1, f;
  double   *Df, *sumD1;

  sp = pkv_GetScratchMemTop ();
  nkn2 = nkn*nkn;
  nf = 6*k+1;
  pder = (point3d*)pkv_GetScratchMemd ( 53 + 22*nf );
  if ( !pder )
    goto failure;
  Gstar  = &pder[11].x;
  DGstar = &Gstar[9];
  Bstar  = &DGstar[9*nf];
  DBstar = &Bstar[11];
  Df     = &DBstar[11*nf];
  sumD1  = &Df[nf];

  memset ( Df, 0, nf*sizeof(double) );
  sum0 = 0.0;
  memset ( gtab, 0, nf*sizeof(double) );
  for ( l = 0; l < k; l++ ) {
    for ( i = knot = 0;  i < nkn;  i++ ) {
      sum1 = 0.0;
      memset ( sumD1, 0, nf*sizeof(double) );
      for ( j = 0;  j < nkn;  j++, knot++ ) {
          /* compute the parameterization derivatives */
        _g2mbl_UCompSDerd ( nkn2, Nitab, knot, k, l, cpind, mvcp, pder );
          /* compute the first fundamental form and its derivatives */
        _g2bl_UCompGStard ( pder, Gstar );
        _g2mbl_SCompDGStard ( nkn2, nf, Nitab, knot, k, l, pder, cpind, nncpi,
                              mvcpn, DGstar );
          /* compute the second fundamental form and its derivatives */
        _g2bl_UCompBStard ( pder, Bstar );
        _g2mbl_SCompDBStard ( nkn2, nf, Nitab, knot, k, l, pder, cpind, nncpi,
                              mvcpn, DBstar );
          /* compute the integrand terms */
        _g2mbl_SFuncGradSQIntegrandd ( nkn2, nf, Nitab, knot, k, l, pder,
                          cpind, nncpi, Gstar, DGstar, Bstar, DBstar, &f, Df );
          /* add the quadrature term */
        sum1 += f*Jac[knot]*qcoeff[j];
        pkn_AddMatrixMd ( 1, nf, 0, sumD1, 0, Df, Jac[knot]*qcoeff[j], 0, sumD1 );
      }
      sum0 += sum1*qcoeff[i];
      pkn_AddMatrixMd ( 1, nf, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
    }
  }
  *ftab = sum0;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_SFuncGradSSQd*/

/* ///////////////////////////////////////////////////////////////////////// */
boolean g2mbl_SFuncGradHessRSQd ( int nkn, double *qcoeff,
                                  double *Nitab, double *Nijtab, double *Mijtab,
                                  int *cpind, int *nncpi, point3d *mvcp, vector3d *mvcpn,
                                  double *ftab, double *gtab, double *htab )
{
  int      nkn2, i, j, knot;
  vector3d pder[11];
  double   Gstar[9], DGstar[9*16], DDGstar[9*136],
           Bstar[11], DBstar[11*16], DDBstar[11*136];
  double   sum0, sum1, f;
  double   Df[16], sumD1[16], DDf[136], sumDD1[136];

  nkn2 = nkn*nkn;
  memset ( Df, 0, 16*sizeof(double) );
  memset ( DDf, 0, 136*sizeof(double) );
  sum0 = 0.0;
  memset ( gtab, 0, 16*sizeof(double) );
  memset ( htab, 0, 136*sizeof(double) );
  for ( i = knot = 0;  i < nkn;  i++ ) {
    sum1 = 0.0;
    memset ( sumD1, 0, 16*sizeof(double) );
    memset ( sumDD1, 0, 136*sizeof(double) );
    for ( j = 0;  j < nkn;  j++, knot++ ) {
        /* compute the parameterization derivatives */
      _g2mbl_UCompRDerd ( nkn2, Nitab, knot, cpind, mvcp, pder );
        /* compute the first fundamental form and its derivatives */
      _g2bl_UCompGStard ( pder, Gstar );
      _g2mbl_SCompDGStard ( nkn2, 16, Nitab, knot, 4, 0, pder, cpind, nncpi,
                            mvcpn, DGstar );
      _g2mbl_SCompDDGstard ( nkn2, 16, Nijtab, knot, 4, 0, cpind, nncpi,
                             mvcpn, DDGstar );
        /* compute the second fundamental form and its derivatives */
      _g2bl_UCompBStard ( pder, Bstar );
      _g2mbl_SCompDBStard ( nkn2, 16, Nitab, knot, 4, 0, pder, cpind, nncpi,
                            mvcpn, DBstar );
      _g2mbl_SCompDDBstard ( nkn2, 16, Mijtab, knot, 4, 0, pder, cpind, nncpi,
                             mvcpn, DDBstar );

        /* compute the integrand terms */
      if ( !_g2mbl_SFuncGradHessSQIntegrandd ( nkn2, 16, Nitab, knot, 4, 0, pder,
                                         cpind, nncpi, Gstar, DGstar, DDGstar,
                                         Bstar, DBstar, DDBstar,
                                         &f, Df, DDf ) )
        return false;
        /* add the quadrature term */
      sum1 += f*qcoeff[j];
      pkn_AddMatrixMd ( 1, 16, 0, sumD1, 0, Df, qcoeff[j], 0, sumD1 );
      pkn_AddMatrixMd ( 1, 136, 0, sumDD1, 0, DDf, qcoeff[j], 0, sumDD1 );
    }
    sum0 += sum1*qcoeff[i];
    pkn_AddMatrixMd ( 1, 16, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
    pkn_AddMatrixMd ( 1, 136, 0, htab, 0, sumDD1, qcoeff[i], 0, htab );
  }
  *ftab = sum0;
  return true;
} /*g2mbl_SFuncGradHessRSQd*/

boolean g2mbl_SFuncGradHessSSQd ( int nkn, double *qcoeff, int k,
                                  double *Nitab, double *Nijtab, double *Mijtab,
                                  double *Jac,
                                  int *cpind, int *nncpi, point3d *mvcp, vector3d *mvcpn,
                                  double *ftab, double *gtab, double *htab )
{
  void     *sp;
  int      nkn2, i, j, l, knot, nf, nh;
  vector3d *pder;
  double   *Gstar, *DGstar, *DDGstar, *Bstar, *DBstar, *DDBstar;
  double   sum0, sum1, f;
  double   *Df, *sumD1, *DDf, *sumDD1;

  sp = pkv_GetScratchMemTop ();
  nkn2 = nkn*nkn;
  nf = 6*k+1;
  nh = (nf*(nf+1))/2;
  pder = (point3d*)pkv_GetScratchMemd ( 53 + 22*nf + 22*nh );
  if ( !pder )
    goto failure;
  Gstar = &pder[11].x;
  DGstar = &Gstar[9];
  DDGstar = &DGstar[9*nf];
  Bstar = &DDGstar[9*nh];
  DBstar = &Bstar[11];
  DDBstar = &DBstar[11*nf];
  Df  = &DDBstar[11*nh];
  DDf = &Df[nf];
  sumD1 = &DDf[nh];
  sumDD1 = &sumD1[nf];

  memset ( Df, 0, nf*sizeof(double) );
  memset ( DDf, 0, nh*sizeof(double) );
  sum0 = 0.0;
  memset ( gtab, 0, nf*sizeof(double) );
  memset ( htab, 0, nh*sizeof(double) );
  for ( l = 0; l < k; l++ ) {
    for ( i = knot = 0;  i < nkn;  i++ ) {
      sum1 = 0.0;
      memset ( sumD1, 0, nf*sizeof(double) );
      memset ( sumDD1, 0, nh*sizeof(double) );
      for ( j = 0;  j < nkn;  j++, knot++ ) {
        /* compute the parameterization derivatives */
        _g2mbl_UCompSDerd ( nkn2, Nitab, knot, k, l, cpind, mvcp, pder );
        /* compute the first fundamental form and its derivatives */
        _g2bl_UCompGStard ( pder, Gstar );
        _g2mbl_SCompDGStard ( nkn2, nf, Nitab, knot, k, l, pder, cpind, nncpi,
                              mvcpn, DGstar );
        _g2mbl_SCompDDGstard ( nkn2, nf, Nijtab, knot, k, l, cpind, nncpi,
                               mvcpn, DDGstar );
        /* compute the second fundamental form and its derivatives */
        _g2bl_UCompBStard ( pder, Bstar );
        _g2mbl_SCompDBStard ( nkn2, nf, Nitab, knot, k, l, pder, cpind, nncpi,
                              mvcpn, DBstar );
        _g2mbl_SCompDDBstard ( nkn2, nf, Mijtab, knot, k, l, pder, cpind, nncpi,
                               mvcpn, DDBstar );
        /* compute the integrand terms */
        if ( !_g2mbl_SFuncGradHessSQIntegrandd ( nkn2, nf, Nitab, knot, k, l, pder,
                                           cpind, nncpi, Gstar, DGstar, DDGstar,
                                           Bstar, DBstar, DDBstar,
                                           &f, Df, DDf ) )
          goto failure;
        /* add the quadrature term */
        sum1 += f*Jac[knot]*qcoeff[j];
        pkn_AddMatrixMd ( 1, nf, 0, sumD1, 0, Df, Jac[knot]*qcoeff[j], 0, sumD1 );
        pkn_AddMatrixMd ( 1, nh, 0, sumDD1, 0, DDf, Jac[knot]*qcoeff[j], 0, sumDD1 );
      }
      sum0 += sum1*qcoeff[i];
      pkn_AddMatrixMd ( 1, nf, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
      pkn_AddMatrixMd ( 1, nh, 0, htab, 0, sumDD1, qcoeff[i], 0, htab );
    }
  }
  *ftab = sum0;
  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2mbl_SFuncGradHessSSQd*/

