
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Mateusz Markowski                                */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "g1blendingd.h"

#include "g1blprivated.h"

/* 
oblicza współczynniki i węzły kwadratury Gaussa-Legendre'a
nkn - liczba węzłów w kwadraturze
*bf - wskaźnik do tablicy wartości funkcji bazowych w węzłach kwadratury
 ///////////////////////////////////////////////////////////////////////// */
boolean _g1bl_TabBasisFuncd ( int nkn, double **knots, double **coeff,
                              double **bf, double **dbf,
                              double **ddbf)
{
        /* bfcp - Bernstein basis coefficients of polynomials, which describe */
        /* a quadratic B-spline function with equidistant (subsequent integer) knots */
  const double bfcp[9] =
    {0.0,0.0,0.5,
     0.5, 1, 0.5,
     0.5,0.0,0.0};
  double *_knots, *_coeff, *_bf, *_dbf, *_ddbf;
  int    i, j, k;

  *knots = _knots = pkv_GetScratchMemd ( nkn ); /* definicja w pkvaria.h alokacja pamięci na nkn double'li*/
  *coeff = _coeff = pkv_GetScratchMemd ( nkn ); /* współczyniki kwadratury */
  *bf    = _bf    = pkv_GetScratchMemd ( 3*nkn );
  *dbf   = _dbf   = pkv_GetScratchMemd ( 3*nkn );
  *ddbf  = _ddbf  = pkv_GetScratchMemd ( 3*nkn );
  if ( !_knots || !_coeff || !_bf || !_dbf || !_ddbf )
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
      /*tablicowanie węzły i współczynniki kwad.*/
    break;
default:
    exit ( 1 );
  }
  for ( j = k = 0;  k < 3;  k++ ) /* numer wielomianu */
    for ( i = 0;  i < nkn;  i++, j++ ) /* wartości ktego wielomianu
 obliczanie funkcji i pochodnych pierwszego i drugiego rzędu */
      if ( !mbs_multiBCHornerDer2d ( 2, 1, 1, 0, &bfcp[3*k], _knots[i],
                                     &_bf[j], &_dbf[j], &_ddbf[j] ) )
        return false;
  return true;
} /*_g1bl_TabBasisFuncd*/

/* ///////////////////////////////////////////////////////////////////////// */
double *_g1bl_NijIndd ( int nkn, double *Nijtab,
                        int i0, int i1, int j0, int j1,
                        int l0, int l1 )
{
  int ir, jr, knot;
/* i0, i1 - numery funkcji 
 */
  ir = 3*i0 + i1;
  jr = 3*j0 + j1;
  knot = nkn*l0 + l1; /* węzeł kwadratury*/
  return &Nijtab[(pkn_SymMatIndex ( ir, jr )*nkn*nkn + knot)*3];
} /*_g2bl_NijIndd*/

/* //////
Wypełnia wektor: Ni = [N_(1,0), N_(0,1), N_(2,0), N_(1,1), N_(0,2) ]^T
/////////////////////////////////////////////////////////////////// */
static void _g1bl_Nid ( int nkn, int i0, int i1, int knu, int knv,
                        const double *bf, const double *dbf,
                        const double *ddbf, double *Ni )
{
  int knui, knvi;

  knui = knu + nkn*(2-i0);
  knvi = knv + nkn*(2-i1);
  Ni[0] = dbf[knui]*bf[knvi];
  Ni[1] = bf[knui]*dbf[knvi];
  Ni[2] = ddbf[knui]*bf[knvi];
  Ni[3] = dbf[knui]*dbf[knvi];
  Ni[4] = bf[knui]*ddbf[knvi];
} /*_g2bl_Nid*/

static double _g1bl_Nijabd ( int knui, int knvi, int knuj, int knvj,
                             const double *bfua, const double *bfva,
                             const double *bfub, const double *bfvb )
{
  return bfua[knui]*bfva[knvi]*bfub[knuj]*bfvb[knvj] +
         bfua[knuj]*bfva[knvj]*bfub[knui]*bfvb[knvi];
} /*_g1bl_Nijabd*/


double *_g1bl_MijIndd ( int nkn, double *Mijtab,
                        int i0, int i1, int j0, int j1, int l0, int l1 )
{
  int ir, jr, knot;

  ir = 3*i0 + i1;
  jr = 3*j0 + j1;
  knot = nkn*l0 + l1;
  return &Mijtab[(pkn_SymMatIndex ( ir, jr-1 )*nkn*nkn + knot)*MijLen];
} /*_g1bl_MijIndd*/


static double _g1bl_Mijabd ( int knui, int knvi, int knuj, int knvj,
                             const double *bfua, const double *bfva,
                             const double *bfub, const double *bfvb )
{
  return bfua[knui]*bfva[knvi]*bfub[knuj]*bfvb[knvj] -
         bfua[knuj]*bfva[knvj]*bfub[knui]*bfvb[knvi];
} /*_g1bl_Mijabd*/

static void _g1bl_Mijd ( int nkn, int i0, int i1, int j0, int j1,
                         int knu, int knv, const double *bf, const double *dbf,
                         const double *ddbf, double *Mij )
{
  int    knui, knvi, knuj, knvj;

  knui = knu + nkn*(2-i0);
  knvi = knv + nkn*(2-i1);
  knuj = knu + nkn*(2-j0);
  knvj = knv + nkn*(2-j1);
  
Mij0102 = _g1bl_Mijabd ( knui, knvi, knuj, knvj, bf, dbf, bf, ddbf );
Mij0111 = _g1bl_Mijabd ( knui, knvi, knuj, knvj, bf, dbf, dbf, dbf );
Mij0120 = _g1bl_Mijabd ( knui, knvi, knuj, knvj, bf, dbf, ddbf, bf );
Mij1001 = _g1bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, bf, dbf );
Mij1002 = _g1bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, bf, ddbf );
Mij1011 = _g1bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, dbf, dbf );
Mij1020 = _g1bl_Mijabd ( knui, knvi, knuj, knvj, dbf, bf, ddbf, bf );
} /*_g1bl_Mijd*/

/* wypełniana jest tablica Mij */
void g1bl_TabMijd ( int nkn, double *bf, double *dbf, double *ddbf, double *Mijtab )
{
  int     i0, i1, j0, j1, l0, l1, ir, jr, knot;
  double  *Mij;

  for ( i0 = 0; i0 <= 2; i0 ++ )
    for ( i1 = 0; i1 <= 2; i1 ++ ) {
      ir = 3*i0+i1;
      for ( j0 = i0; j0 <= 2; j0++ )
        for ( j1 = 0; j1 <= 2; j1++ ) {
          jr = 3*j0+j1;
          if ( ir < jr ) {
            Mij = _g1bl_MijIndd ( nkn, Mijtab, i0, i1, j0, j1, 0, 0 );
            for ( l0 = knot = 0;  l0 < nkn;  l0++ )
              for ( l1 = 0;  l1 < nkn;  l1++, knot++ )
                _g1bl_Mijd ( nkn, i0, i1, j0, j1, l0, l1,
                             bf, dbf, ddbf, &Mij[MijLen*knot] );
          }
        }
    }
} /*g1bl_TabMijd*/



/*
 * wypełniana macierz N_{*}^{ij}
dla i=(i0,i1) j=(j0,j1)
 */
static void _g1bl_Nijd ( int nkn, int i0, int i1, int j0, int j1,
                         int knu, int knv,
                         const double *bf, const double *dbf, const double *ddbf,
                         double *Nij )
{
  int    knui, knvi, knuj, knvj;
/*  double N1110, N1101;*/
  knui = knu + nkn*(2-i0); /*i0 \in {0,1,2} itd */
  knvi = knv + nkn*(2-i1);
  knuj = knu + nkn*(2-j0);
  knvj = knv + nkn*(2-j1);
  Nij[0] = _g1bl_Nijabd ( knui, knvi, knuj, knvj, dbf, bf, dbf, bf );
  Nij[1] = _g1bl_Nijabd ( knui, knvi, knuj, knvj, bf, dbf, dbf, bf );
  Nij[2] = _g1bl_Nijabd ( knui, knvi, knuj, knvj, bf, dbf, bf, dbf );
} /*_g1bl_Nijd*/

/*

*/
void g1bl_TabNid ( int nkn, double *bf, double *dbf, double *ddbf,
                   double *Nitab )
{
  int i0, i1, l0, l1, k;
/*
obliczana jest dla każdego z węzłów każdego z 9ciu kwadratów 
*/
  for ( i0 = k = 0;  i0 <= 2;  i0++ )
    for ( i1 = 0;  i1 <= 2;  i1++ )
      for ( l0 = 0;  l0 < nkn;  l0++ )
        for ( l1 = 0;  l1 < nkn;  l1++, k += 5 ) /*5 - długość N_i */
          _g1bl_Nid ( nkn, i0, i1, l0, l1, bf, dbf, ddbf, &Nitab[k] );
} /*g2bl_TabNid*/

void g1bl_TabNijd ( int nkn, double *bf, double *dbf, double *ddbf,
                    double *Nijtab )
{
  int    i0, i1, j0, j1, l0, l1;
  double *Nij;

  for ( i0 = 0; i0 <= 2; i0++ )
    for ( i1 = 0; i1 <= 2; i1++ )
      for ( j0 = i0; j0 <= 2; j0++ )
        for ( j1 = 0; j1 <= 2; j1++ )
          if ( 3*i0+i1 <= 3*j0+j1 )
            for ( l0 = 0; l0 < nkn; l0++ )
              for ( l1 = 0; l1 < nkn; l1++ ) {
                Nij = _g1bl_NijIndd ( nkn, Nijtab, i0, i1, j0, j1, l0, l1 );
                _g1bl_Nijd ( nkn, i0, i1, j0, j1, l0, l1,
                             bf, dbf, ddbf, Nij );
              }
} /*g2bl_TabNijd*/
 
/* ///////// obliczanie pochodnych do drugiego stopnia włącznie
jest ich 5 - 2 pochodne 1 rzędu i 3 pochodne 2 go rzędu
fcpn - first control point 
cp - siatka kontrolna (tablica w której znajdują sie punkty kontrolne powierzchni)
i,j numer węzła w kolumnie i wierszu w odpowiednim kwadracie 
pder - tablica na pochodne
i,j < nkn
 //////////////////////////////////////////////////////////////// */
void _g1bl_UCompPDerd ( int nkn, double *Nitab,
                        int pitch, point3d *cp,
                        int fcpn, int i, int j, vector3d *pder )
{
  int      k, l;
  vector3d *_cp;
  double   *Ni;

  memset ( pder, 0, 5*sizeof(vector3d) );
  for ( k = 0; k < 3; k++ )
    for ( l = 0; l < 3; l++ ) {
      _cp = &cp[fcpn+k*pitch+l];
      Ni = &Nitab[((((k*3)+l)*nkn + i)*nkn + j)*5];
      pkn_MultMatrixAddd ( 5, 1, 1, Ni, 3, 0, &_cp->x, 3, &pder[0].x );
    }
        /* compute the Laplacian gradient */
/*  AddVector3d ( &puuu, &puvv, &pder[9] );
  AddVector3d ( &puuv, &pvvv, &pder[10] );*/
} /*_g1bl_UCompPDerd*/

/* ///////////////////////////////////////////////////////////////////////// */
void _g1bl_UCompGStard ( const vector3d *pder, double *Gstar )
{
  g11 = DotProduct3d ( &pu, &pu );
  g12 = DotProduct3d ( &pu, &pv );
  g22 = DotProduct3d ( &pv, &pv );
} /*_g1bl_UCompGStard*/

void _g1bl_UCompDGStard ( int nkn, double *Nitab,
                          int lastknotu, int lastknotv,
                          int ip0, int ip1, int jp0, int jp1,
                          int isq, int jsq, int i, int j,
                          const vector3d *pder, double *DGstar )
{
  int      ii, jj, kk;
  vector3d *Dgstar;
  double   *Ni;

  for ( ii = 0;  ii <= 2;  ii++ )
    if ( ii+isq >= ip0 && ii+isq < ip1 )
      for ( jj = 0;  jj <= 2;  jj++ )
        if ( jj+jsq >= jp0 && jj+jsq < jp1 ) {
          kk = 3*ii + jj;
          Dgstar = (vector3d*)&DGstar[kk*3*3];
          Ni = &Nitab[((kk*nkn + i)*nkn + j)*5];

          MultVector3d ( 2.0*Ni10, &pu, &dg11 );
          MultVector3d ( Ni10, &pv, &dg12 );
          AddVector3Md ( &dg12, &pu, Ni01, &dg12 );
          MultVector3d ( 2.0*Ni01, &pv, &dg22 );
        }
} /*_g1bl_UCompDGStard*/

void _g1bl_UCompDDGStard ( int nkn, double *Nijtab,
                           int lastknotu, int lastknotv,
                           int ip0, int ip1, int jp0, int jp1,
                           int isq, int jsq, int i, int j,
                           double *DDGstar )
{
  int      i0, i1, ii, ir, j0, j1, ji, jr;
  double   *Nij;
/*
isq - numer kolumny w której jest
jsq 	 wiersza
*/

        /* now run through the basis functions nonzero at this knot */
  for ( i0 = isq, ir = 0;  i0 < isq+3;  i0++ )
    for ( i1 = jsq;  i1 < jsq+3;  i1++, ir++ ) {
      if ( i0 >= ip0 && i0 < ip1 && i1 >= jp0 && i1 < jp1 ) {
        ii = (i0-ip0)*(lastknotv-6) + (i1-jp0);  /* number of point unknown variable */
/*
 6 - liczba punktów kontrolnych brzegowych
*/
        for ( j0 = isq, jr = 0;  j0 < isq+3;  j0++ )
          for ( j1 = jsq;  j1 < jsq+3;  j1++, jr++ ) {
            if ( j0 >= ip0 && j0 < ip1 && j1 >= jp0 && j1 < jp1 ) {
              ji = (j0-ip0)*(lastknotv-6) + (j1-jp0);  /* number of point unknown variable */

/*liczenie Hesjanu - symetria */
              if ( ji >= ii ) {
                Nij = _g1bl_NijIndd ( nkn, Nijtab,
                                      i0-isq, i1-jsq, j0-isq, j1-jsq, i, j );
                memcpy ( &DDGstar[3*pkn_SymMatIndex(ir,jr)], Nij, 3*sizeof(double) );
              }
            }
          }
      }
    }
} /*_g1bl_UCompDDGStard*/

/* ///////////////////////////////////////////////////////////////////////// */
void _g1bl_UCompBStard ( const vector3d *pder, double *Bstar )
{
  /*double huuu, huuv, huvv, hvvv;*/

  tb11 = det3d ( &puu, &pu, &pv );
  tb12 = det3d ( &puv, &pu, &pv );
  tb22 = det3d ( &pvv, &pu, &pv ); 
} /*_g1bl_UCompBStard*/


/* obliczanie na podstawie wzorów ze strony 19tej*/
void _g1bl_UCompDBStard ( int nkn, double *Nitab,
                          int lastknotu, int lastknotv,
                          int ip0, int ip1, int jp0, int jp1,
                          int isq, int jsq, int i, int j,
                          const vector3d *pder, double *DBstar )
{
  int      ii, jj, kk;
  vector3d *Dbstar, v, w, z, y;
  double   *Ni;

  for ( ii = 0;  ii <= 2;  ii++ )
    if ( ii+isq >= ip0 && ii+isq < ip1 )
      for ( jj = 0;  jj <= 2;  jj++ )
        if ( jj+jsq >= jp0 && jj+jsq < jp1 ) {
          kk = 3*ii + jj;
          Dbstar = (vector3d*)&DBstar[kk*3*3];
          Ni = &Nitab[((kk*nkn + i)*nkn + j)*5];

/* należy jeszcze sprawdzić dokładnie gdzie kończy się liczenie tych rzeczy \/*/
          CrossProduct3d ( &pu, &pv, &v );
          MultVector3d ( Ni20, &v, &dtb11 );
          MultVector3d ( Ni11, &v, &dtb12 );
          MultVector3d ( Ni02, &v, &dtb22 );
          /*MultVector3d ( Ni30, &v, &dtb11u );
          MultVector3d ( Ni21, &v, &dtb12u );
          MultVector3d ( Ni12, &v, &dtb22u );
          MultVector3d ( Ni03, &v, &dtb22v );
          AddVector3d ( &dtb11u, &dtb22u, &Dbstar[9] );
          AddVector3d ( &dtb12u, &dtb22v, &Dbstar[10] );*/

          CrossProduct3d ( &puu, &pu, &v );
          AddVector3Md ( &dtb11, &v, Ni01, &dtb11 );
          /*AddVector3Md ( &dtb11u, &v, Ni11, &dtb11u );
          MultVector3d ( Ni02, &v, &dtb11v );*/

          CrossProduct3d ( &pv, &puu, &v );
          AddVector3Md ( &dtb11, &v, Ni10, &dtb11 );
          MultVector3d ( Ni11, &v, &z );
          MultVector3d ( -Ni02, &v, &y );

          CrossProduct3d ( &pv, &puv, &v );
          AddVector3Md ( &dtb12, &v, Ni10, &dtb12 );
          AddVector3Md ( &z, &v, -Ni20, &z );
          /*AddVector3Md ( &dtb22v, &v, -Ni02, &dtb22v );*/

          CrossProduct3d ( &puv, &pu, &v );
          AddVector3Md ( &dtb12, &v, Ni01, &dtb12 );
          /*AddVector3Md ( &dtb11u, &v, -Ni20, &dtb11u );*/
          MultVector3d ( Ni02, &v, &w );

          CrossProduct3d ( &pv, &pvv, &v );
          AddVector3Md ( &dtb22, &v, Ni10, &dtb22 );
          /*AddVector3Md ( &dtb22v, &v, Ni11, &dtb22v );*/
          AddVector3Md ( &y, &v, Ni20, &y );

          CrossProduct3d ( &pvv, &pu, &v );
          AddVector3Md ( &dtb22, &v, Ni01, &dtb22 );
          /*AddVector3Md ( &dtb11v, &v, -Ni20, &dtb11v );*/
          AddVector3Md ( &w, &v, -Ni11, &w );

          /*CrossProduct3d ( &pv, &puuu, &v );
          AddVector3Md ( &dtb11u, &v, Ni10, &dtb11u );
          CrossProduct3d ( &puuu, &pu, &v );
          AddVector3Md ( &dtb11u, &v, Ni01, &dtb11u );*/

          CrossProduct3d ( &puv, &puu, &v );
          /*AddVector3Md ( &dtb11u, &v, Ni10, &dtb11u );*/
          AddVector3Md ( &z, &v, -Ni01, &z );

          /*CrossProduct3d ( &pv, &puuv, &v );
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
          AddVector3Md ( &Dbstar[10], &v, Ni01, &Dbstar[10] ); */
        }
} /*_g1bl_UCompDBStard*/



void _g1bl_UCompDDBStard ( int nkn, double *Mijtab,
                           int lastknotu, int lastknotv,
                           int ip0, int ip1, int jp0, int jp1,
                           int isq, int jsq, int i, int j,
                           const vector3d *pder, double *DDBstar )
{
#define ddtb11  DDB[0]
#define ddtb12  DDB[1]
#define ddtb22  DDB[2]

  int     i0, i1, ii, ir, j0, j1, ji, jr;
  double  *Mij;
  vector3d *DDB;

  for ( i0 = isq, ir = 0;  i0 < isq+3;  i0++ )
    for ( i1 = jsq;  i1 < jsq+3;  i1++, ir++ ) {
      if ( i0 >= ip0 && i0 < ip1 && i1 >= jp0 && i1 < jp1 ) {
        ii = (i0-ip0)*(lastknotv-6) + (i1-jp0); 
        DDB = (vector3d*)&DDBstar[3*3*pkn_SymMatIndex(ir,ir)];
        memset ( DDB, 0, 3*3*sizeof(double) );
        for ( j0 = isq, jr = 0;  j0 < isq+3;  j0++ )
          for ( j1 = jsq;  j1 < jsq+3;  j1++, jr++ ) {
            ji = (j0-ip0)*(lastknotv-6) + (j1-jp0);
            if ( ji > ii ) {
              DDB = (vector3d*)&DDBstar[3*3*pkn_SymMatIndex(ir,jr)];
              Mij = _g1bl_MijIndd ( nkn, Mijtab, i0-isq, i1-jsq, j0-isq, j1-jsq, i, j );
	      
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
            }
          }
      }
    }
#undef ddtb11
#undef ddtb12
#undef ddtb22

} /*_g1bl_UCompDDBStard*/
