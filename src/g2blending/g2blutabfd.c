
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009, 2013                            */
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
boolean _g2bl_TabBasisFuncd ( int nkn, double **knots, double **coeff,
                              double **bf, double **dbf,
                              double **ddbf, double **dddbf )
{
        /* bfcp - Bernstein basis coefficients of polynomials, which describe */
        /* a cubic B-spline function with equidistant (subsequent integer) knots */
  const double bfcp[16] =
    {0.0,0.0,0.0,0.16666666666666666,
     0.16666666666666666,0.33333333333333333,0.66666666666666666,0.66666666666666666,
     0.66666666666666666,0.66666666666666666,0.33333333333333333,0.16666666666666666,
     0.16666666666666666,0.0,0.0,0.0};
  double *_knots, *_coeff, *_bf, *_dbf, *_ddbf, *_dddbf;
  int    i, j, k;

  *knots = _knots = pkv_GetScratchMemd ( nkn );
  *coeff = _coeff = pkv_GetScratchMemd ( nkn );
  *bf    = _bf    = pkv_GetScratchMemd ( 4*nkn );
  *dbf   = _dbf   = pkv_GetScratchMemd ( 4*nkn );
  *ddbf  = _ddbf  = pkv_GetScratchMemd ( 4*nkn );
  *dddbf = _dddbf = pkv_GetScratchMemd ( 4*nkn );
  if ( !_knots || !_coeff || !_bf || !_dbf || !_ddbf || !_dddbf )
    return false;

  switch ( nkn ) {
case 3:
    pkn_QuadGaussLegendre6d ( 0.0, 1.0, nkn, _knots, _coeff );
    break;
case 4:
    pkn_QuadGaussLegendre8d ( 0.0, 1.0, nkn, _knots, _coeff );
    break;
case 5:
    pkn_QuadGaussLegendre10d ( 0.0, 1.0, nkn, _knots, _coeff );
    break;
case 6:
    pkn_QuadGaussLegendre12d ( 0.0, 1.0, nkn, _knots, _coeff );
    break;
case 7:
    pkn_QuadGaussLegendre14d ( 0.0, 1.0, nkn, _knots, _coeff );
    break;
case 8:
    pkn_QuadGaussLegendre16d ( 0.0, 1.0, nkn, _knots, _coeff );
    break;
case 9:
    pkn_QuadGaussLegendre18d ( 0.0, 1.0, nkn, _knots, _coeff );
    break;
case 10:
    pkn_QuadGaussLegendre20d ( 0.0, 1.0, nkn, _knots, _coeff );
    break;
default:
    return false;
  }
  for ( j = k = 0;  k < 4;  k++ )
    for ( i = 0;  i < nkn;  i++, j++ )
      if ( !mbs_multiBCHornerDer3d ( 3, 1, 1, 0, &bfcp[4*k], _knots[i],
                                     &_bf[j], &_dbf[j], &_ddbf[j], &_dddbf[j] ) )
        return false;

  return true;
} /*_g2bl_TabBasisFuncd*/

/* ///////////////////////////////////////////////////////////////////////// */
double *_g2bl_NijIndd ( int nkn, double *Nijtab,
                        int i0, int i1, int j0, int j1,
                        int l0, int l1 )
{
  int ir, jr, knot;

  ir = 4*i0 + i1;
  jr = 4*j0 + j1;
  knot = nkn*l0 + l1;
  return &Nijtab[(pkn_SymMatIndex ( ir, jr )*nkn*nkn + knot)*9];
} /*_g2bl_NijIndd*/

double *_g2bl_MijIndd ( int nkn, double *Mijtab,
                        int i0, int i1, int j0, int j1, int l0, int l1 )
{
  int ir, jr, knot;

  ir = 4*i0 + i1;
  jr = 4*j0 + j1;
  knot = nkn*l0 + l1;
  return &Mijtab[(pkn_SymMatIndex ( ir, jr-1 )*nkn*nkn + knot)*18];
} /*_g2bl_MijIndd*/

/* ///////////////////////////////////////////////////////////////////////// */
static void _g2bl_Nid ( int nkn, int i0, int i1, int knu, int knv,
                        const double *bf, const double *dbf,
                        const double *ddbf, const double *dddbf,
                        double *Ni )
{
  int knui, knvi;

  knui = knu + nkn*(3-i0);
  knvi = knv + nkn*(3-i1);
  Ni[0] = dbf[knui]*bf[knvi];
  Ni[1] = bf[knui]*dbf[knvi];
  Ni[2] = ddbf[knui]*bf[knvi];
  Ni[3] = dbf[knui]*dbf[knvi];
  Ni[4] = bf[knui]*ddbf[knvi];
  Ni[5] = dddbf[knui]*bf[knvi];
  Ni[6] = ddbf[knui]*dbf[knvi];
  Ni[7] = dbf[knui]*ddbf[knvi];
  Ni[8] = bf[knui]*dddbf[knvi];
} /*_g2bl_Nid*/

static double _g2bl_Nijabd ( int knui, int knvi, int knuj, int knvj,
                             const double *bfua, const double *bfva,
                             const double *bfub, const double *bfvb )
{
  return bfua[knui]*bfva[knvi]*bfub[knuj]*bfvb[knvj] +
         bfua[knuj]*bfva[knvj]*bfub[knui]*bfvb[knvi];
} /*_g2bl_Nijabd*/

static void _g2bl_Nijd ( int nkn, int i0, int i1, int j0, int j1,
                         int knu, int knv,
                         const double *bf, const double *dbf, const double *ddbf,
                         double *Nij )
{
  int    knui, knvi, knuj, knvj;
  double N1110, N1101;

  knui = knu + nkn*(3-i0);
  knvi = knv + nkn*(3-i1);
  knuj = knu + nkn*(3-j0);
  knvj = knv + nkn*(3-j1);
  Nij[0] = _g2bl_Nijabd ( knui, knvi, knuj, knvj, dbf, bf, dbf, bf );
  Nij[1] = _g2bl_Nijabd ( knui, knvi, knuj, knvj, bf, dbf, dbf, bf );
  Nij[2] = _g2bl_Nijabd ( knui, knvi, knuj, knvj, bf, dbf, bf, dbf );
  Nij[3] = 2.0*_g2bl_Nijabd ( knui, knvi, knuj, knvj, ddbf, bf, dbf, bf );
  N1110  = _g2bl_Nijabd ( knui, knvi, knuj, knvj, dbf, dbf, dbf, bf );
  N1101  = _g2bl_Nijabd ( knui, knvi, knuj, knvj, dbf, dbf, bf, dbf );
  Nij[4] = _g2bl_Nijabd ( knui, knvi, knuj, knvj, ddbf, bf, bf, dbf ) +
           N1110;
  Nij[5] = 2.0*N1101;
  Nij[6] = 2.0*N1110;
  Nij[7] = _g2bl_Nijabd ( knui, knvi, knuj, knvj, bf, ddbf, dbf, bf ) +
           N1101;
  Nij[8] = 2.0*_g2bl_Nijabd ( knui, knvi, knuj, knvj, bf, ddbf, bf, dbf );
} /*_g2bl_Nijd*/

static double _g2bl_Mijabd ( int knui, int knvi, int knuj, int knvj,
                             const double *bfua, const double *bfva,
                             const double *bfub, const double *bfvb )
{
  return bfua[knui]*bfva[knvi]*bfub[knuj]*bfvb[knvj] -
         bfua[knuj]*bfva[knvj]*bfub[knui]*bfvb[knvi];
} /*_g2bl_Mijabd*/

static void _g2bl_Mijd ( int nkn, int i0, int i1, int j0, int j1,
                         int knu, int knv, const double *bf, const double *dbf,
                         const double *ddbf, const double *dddbf,
                         double *Mij )
{
  int    knui, knvi, knuj, knvj;

  knui = knu + nkn*(3-i0);
  knvi = knv + nkn*(3-i1);
  knuj = knu + nkn*(3-j0);
  knvj = knv + nkn*(3-j1);
  Mij1001 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, bf, dbf );
  Mij1020 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, ddbf, bf );
  Mij0120 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, bf, dbf, ddbf, bf );
  Mij1011 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, dbf, dbf );
  Mij0111 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, bf, dbf, dbf, dbf );
  Mij1002 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, bf, ddbf );
  Mij0102 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, bf, dbf, bf, ddbf );
  Mij1030 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, dddbf, bf );
  Mij0130 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, bf, dbf, dddbf, bf );
  Mij1021 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, ddbf, dbf );
  Mij0121 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, bf, dbf, ddbf, dbf );
  Mij1012 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, dbf, ddbf );
  Mij0112 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, bf, dbf, dbf, ddbf );
  Mij1003 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, bf, dddbf );
  Mij0103 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, bf, dbf, bf, dddbf );
  Mij2011 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, ddbf, bf, dbf, dbf );
  Mij2002 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, ddbf, bf, bf, ddbf );
  Mij1102 = _g2bl_Mijabd ( knui, knvi, knuj, knvj, dbf, dbf, bf, ddbf );
} /*_g2bl_Mijd*/

void g2bl_TabNid ( int nkn, double *bf, double *dbf, double *ddbf, double *dddbf,
                   double *Nitab )
{
  int i0, i1, l0, l1, k;

  for ( i0 = k = 0;  i0 <= 3;  i0++ )
    for ( i1 = 0;  i1 <= 3;  i1++ )
      for ( l0 = 0;  l0 < nkn;  l0++ )
        for ( l1 = 0;  l1 < nkn;  l1++, k += 9 )
          _g2bl_Nid ( nkn, i0, i1, l0, l1, bf, dbf, ddbf, dddbf, &Nitab[k] );
} /*g2bl_TabNid*/

void g2bl_TabNijd ( int nkn, double *bf, double *dbf, double *ddbf,
                    double *Nijtab )
{
  int    i0, i1, j0, j1, l0, l1;
  double *Nij;

  for ( i0 = 0; i0 <= 3; i0++ )
    for ( i1 = 0; i1 <= 3; i1++ )
      for ( j0 = i0; j0 <= 3; j0++ )
        for ( j1 = 0; j1 <= 3; j1++ )
          if ( 4*i0+i1 <= 4*j0+j1 )
            for ( l0 = 0; l0 < nkn; l0++ )
              for ( l1 = 0; l1 < nkn; l1++ ) {
                Nij = _g2bl_NijIndd ( nkn, Nijtab, i0, i1, j0, j1, l0, l1 );
                _g2bl_Nijd ( nkn, i0, i1, j0, j1, l0, l1,
                             bf, dbf, ddbf, Nij );
              }
} /*g2bl_TabNijd*/

void g2bl_TabMijd ( int nkn, double *bf, double *dbf, double *ddbf, double *dddbf,
                    double *Mijtab )
{
  int     i0, i1, j0, j1, l0, l1, ir, jr, knot;
  double  *Mij;

  for ( i0 = 0; i0 <= 3; i0 ++ )
    for ( i1 = 0; i1 <= 3; i1 ++ ) {
      ir = 4*i0+i1;
      for ( j0 = i0; j0 <= 3; j0++ )
        for ( j1 = 0; j1 <= 3; j1++ ) {
          jr = 4*j0+j1;
          if ( ir < jr ) {
            Mij = _g2bl_MijIndd ( nkn, Mijtab, i0, i1, j0, j1, 0, 0 );
            for ( l0 = knot = 0;  l0 < nkn;  l0++ )
              for ( l1 = 0;  l1 < nkn;  l1++, knot++ )
                _g2bl_Mijd ( nkn, i0, i1, j0, j1, l0, l1,
                             bf, dbf, ddbf, dddbf, &Mij[18*knot] );
          }
        }
    }
} /*g2bl_TabMijd*/

/* ///////////////////////////////////////////////////////////////////////// */
void _g2bl_UCompPDerd ( int nkn, double *Nitab,
                        int pitch, point3d *cp,
                        int fcpn, int i, int j, vector3d *pder )
{
  int      k, l;
  vector3d *_cp;
  double   *Ni;

  memset ( pder, 0, 9*sizeof(vector3d) );
  for ( k = 0; k < 4; k++ )
    for ( l = 0; l < 4; l++ ) {
      _cp = &cp[fcpn+k*pitch+l];
      Ni = &Nitab[((((k*4)+l)*nkn + i)*nkn + j)*9];
      pkn_MultMatrixAddd ( 9, 1, 1, Ni, 3, 0, &_cp->x, 3, &pder[0].x );
    }
        /* compute the Laplacian gradient */
  AddVector3d ( &puuu, &puvv, &pder[9] );
  AddVector3d ( &puuv, &pvvv, &pder[10] );
} /*_g2bl_UCompPDerd*/

/* ///////////////////////////////////////////////////////////////////////// */
void _g2bl_UCompGStard ( const vector3d *pder, double *Gstar )
{
  g11 = DotProduct3d ( &pu, &pu );
  g12 = DotProduct3d ( &pu, &pv );
  g22 = DotProduct3d ( &pv, &pv );
  g11u = 2.0*DotProduct3d ( &pu, &puu );
  g12u = DotProduct3d ( &puu, &pv ) + DotProduct3d ( &pu, &puv );
  g22u = 2.0*DotProduct3d ( &pv, &puv );
  g11v = 2.0*DotProduct3d ( &puv, &pu );
  g12v = DotProduct3d ( &puv, &pv ) + DotProduct3d ( &pu, &pvv );
  g22v = 2.0*DotProduct3d ( &pv, &pvv );
} /*_g2bl_UCompGStard*/

void _g2bl_UCompDGStard ( int nkn, double *Nitab,
                          int lastknotu, int lastknotv,
                          int ip0, int ip1, int jp0, int jp1,
                          int isq, int jsq, int i, int j,
                          const vector3d *pder, double *DGstar )
{
  int      ii, jj, kk;
  vector3d *Dgstar;
  double   *Ni;

  for ( ii = 0;  ii <= 3;  ii++ )
    if ( ii+isq >= ip0 && ii+isq < ip1 )
      for ( jj = 0;  jj <= 3;  jj++ )
        if ( jj+jsq >= jp0 && jj+jsq < jp1 ) {
          kk = 4*ii + jj;
          Dgstar = (vector3d*)&DGstar[kk*3*9];
          Ni = &Nitab[((kk*nkn + i)*nkn + j)*9];

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
} /*_g2bl_UCompDGStard*/

void _g2bl_UCompDDGStard ( int nkn, double *Nijtab,
                           int lastknotu, int lastknotv,
                           int ip0, int ip1, int jp0, int jp1,
                           int isq, int jsq, int i, int j,
                           double *DDGstar )
{
  int      i0, i1, ii, ir, j0, j1, ji, jr;
  double   *Nij;

        /* now run through the basis functions nonzero at this knot */
  for ( i0 = isq, ir = 0;  i0 < isq+4;  i0++ )
    for ( i1 = jsq;  i1 < jsq+4;  i1++, ir++ ) {
      if ( i0 >= ip0 && i0 < ip1 && i1 >= jp0 && i1 < jp1 ) {
        ii = (i0-ip0)*(lastknotv-9) + (i1-jp0);  /* number of point unknown variable */
        for ( j0 = isq, jr = 0;  j0 < isq+4;  j0++ )
          for ( j1 = jsq;  j1 < jsq+4;  j1++, jr++ ) {
            if ( j0 >= ip0 && j0 < ip1 && j1 >= jp0 && j1 < jp1 ) {
              ji = (j0-ip0)*(lastknotv-9) + (j1-jp0);  /* number of point unknown variable */
              if ( ji >= ii ) {
                Nij = _g2bl_NijIndd ( nkn, Nijtab,
                                      i0-isq, i1-jsq, j0-isq, j1-jsq, i, j );
                memcpy ( &DDGstar[9*pkn_SymMatIndex(ir,jr)], Nij, 9*sizeof(double) );
              }
            }
          }
      }
    }
} /*_g2bl_UCompDDGStard*/

/* ///////////////////////////////////////////////////////////////////////// */
void _g2bl_UCompBStard ( const vector3d *pder, double *Bstar )
{
  double huuu, huuv, huvv, hvvv;

  tb11 = det3d ( &puu, &pu, &pv );
  tb12 = det3d ( &puv, &pu, &pv );
  tb22 = det3d ( &pvv, &pu, &pv );
  huuu = det3d ( &puuu, &pu, &pv );
  huuv = det3d ( &puuv, &pu, &pv );
  huvv = det3d ( &puvv, &pu, &pv );
  hvvv = det3d ( &pvvv, &pu, &pv );
  tb11u = huuu + det3d ( &puu, &pu, &puv );
  tb12u = huuv + det3d ( &puv, &puu, &pv );
  tb22u = huvv + det3d ( &pvv, &puu, &pv ) + det3d ( &pvv, &pu, &puv );
  tb11v = huuv + det3d ( &puu, &puv, &pv ) + det3d ( &puu, &pu, &pvv );
  tb12v = huvv + det3d ( &puv, &pu, &pvv );
  tb22v = hvvv + det3d ( &pvv, &puv, &pv );
  Bstar[9]  = huuu + huvv;
  Bstar[10] = huuv + hvvv;
} /*_g2bl_UCompBStard*/

void _g2bl_UCompDBStard ( int nkn, double *Nitab,
                          int lastknotu, int lastknotv,
                          int ip0, int ip1, int jp0, int jp1,
                          int isq, int jsq, int i, int j,
                          const vector3d *pder, double *DBstar )
{
  int      ii, jj, kk;
  vector3d *Dbstar, v, w, z, y;
  double   *Ni;

  for ( ii = 0;  ii <= 3;  ii++ )
    if ( ii+isq >= ip0 && ii+isq < ip1 )
      for ( jj = 0;  jj <= 3;  jj++ )
        if ( jj+jsq >= jp0 && jj+jsq < jp1 ) {
          kk = 4*ii + jj;
          Dbstar = (vector3d*)&DBstar[kk*3*11];
          Ni = &Nitab[((kk*nkn + i)*nkn + j)*9];

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
} /*_g2bl_UCompDBStard*/

void _g2bl_UCompDDBStard ( int nkn, double *Mijtab,
                           int lastknotu, int lastknotv,
                           int ip0, int ip1, int jp0, int jp1,
                           int isq, int jsq, int i, int j,
                           const vector3d *pder, double *DDBstar )
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

  int     i0, i1, ii, ir, j0, j1, ji, jr;
  double  *Mij;
  vector3d *DDB;
  vector3d DDtb301001, DDtb211001, DDtb121001, DDtb031001, DDtb201011,
           DDtb112001, DDtb022001, DDtb201002, DDtb111002, DDtb021101;

  for ( i0 = isq, ir = 0;  i0 < isq+4;  i0++ )
    for ( i1 = jsq;  i1 < jsq+4;  i1++, ir++ ) {
      if ( i0 >= ip0 && i0 < ip1 && i1 >= jp0 && i1 < jp1 ) {
        ii = (i0-ip0)*(lastknotv-9) + (i1-jp0);
        DDB = (vector3d*)&DDBstar[11*3*pkn_SymMatIndex(ir,ir)];
        memset ( DDB, 0, 11*3*sizeof(double) );
        for ( j0 = isq, jr = 0;  j0 < isq+4;  j0++ )
          for ( j1 = jsq;  j1 < jsq+4;  j1++, jr++ ) {
            ji = (j0-ip0)*(lastknotv-9) + (j1-jp0);
            if ( ji > ii ) {
              DDB = (vector3d*)&DDBstar[11*3*pkn_SymMatIndex(ir,jr)];
              Mij = _g2bl_MijIndd ( nkn, Mijtab,
                                    i0-isq, i1-jsq, j0-isq, j1-jsq, i, j );
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
} /*_g2bl_UCompDDBStard*/

