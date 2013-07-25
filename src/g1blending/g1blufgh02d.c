
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010                                  */
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

double g1bl_biharmFuncd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv, int pitch, point3d *cp,
                     char *dirty,
                     double tC, double *ftab )
{
  int    isq, jsq, sqk;
  double f;

  pitch /= 3;     /* express the pitch in terms of points */
        /* compute the integrals in dirty squares */
  if ( dirty ) {
    for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ )
        if ( dirty[sqk] & DIRTY_FUNC ) {
          g1bl_biharmFuncSQd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch, cp,
                          tC, isq, jsq, &ftab[sqk] );
          dirty[sqk] &= ~DIRTY_FUNC;
        }
  }
  else {
    for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ )
        g1bl_biharmFuncSQd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch, cp,
                        tC, isq, jsq, &ftab[sqk] );
  }
        /* sum the integrals */
  f = 0.0;
  for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
    for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ )
      f += ftab[sqk];
  return f;
} /*g1bl_biharmFuncd*/

double g1bl_QFuncd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv, int pitch, point3d *cp,
                     char *dirty,
                     double tC, double *ftab )
{
  int    isq, jsq, sqk;
  double f;

  pitch /= 3;     /* express the pitch in terms of points */
        /* compute the integrals in dirty squares */
  if ( dirty ) {
    for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ )
        if ( dirty[sqk] & DIRTY_FUNC ) {
          g1bl_QFuncSQd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch, cp,
                          tC, isq, jsq, &ftab[sqk] );
          dirty[sqk] &= ~DIRTY_FUNC;
        }
  }
  else {
    for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ )
        g1bl_QFuncSQd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch, cp,
                        tC, isq, jsq, &ftab[sqk] );
  }
        /* sum the integrals */
  f = 0.0;
  for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
    for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ )
      f += ftab[sqk];
  return f;
} /*g1bl_QFuncd*/


/* uaktualnianie wartosci funkcjonalu gradientu i hesjanu dla brudnych kwadratow*/
/* ///////////////////////////////////////////////////////////////////////// */
double g1bl_UFuncd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv, int pitch, point3d *cp,
                     char *dirty,
                     double tC, double *ftab )
{
  int    isq, jsq, sqk;
  double f;

  pitch /= 3;     /* express the pitch in terms of points */
        /* compute the integrals in dirty squares */
  if ( dirty ) {
    for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ )
        if ( dirty[sqk] & DIRTY_FUNC ) {
          g1bl_UFuncSQd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch, cp,
                          tC, isq, jsq, &ftab[sqk] );
          dirty[sqk] &= ~DIRTY_FUNC;
        }
  }
  else {
    for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ )
        g1bl_UFuncSQd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch, cp,
                        tC, isq, jsq, &ftab[sqk] );
  }
        /* sum the integrals */
  f = 0.0;
  for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
    for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ )
      f += ftab[sqk];
  return f;
} /*g1bl_UFuncd*/

/* ////////////////////////////////
ftab - dlugosc taka jak liczba kwadratow
gtab - 9 * 3 * nsq ktore sluza do obliczenia gradientu
///////////////////////////////////////// */
void g1bl_UFuncGradd ( int nkn, const double *qcoeff, double *Nitab,
                       int lastknotu, int lastknotv,
                       int pitch, point3d *cp, char *dirty,
                       double tC, double *ftab, double *gtab,
                       double *func, double *grad )
{
  double f;
  int    nicp, nvars;
  int    isq, jsq, sqk;
  int    i, j, ip, jp, kp;

  pitch /= 3;             /* express the pitch in terms of points */
  /* poprawilem juz to z 3*3 na 3*2*/
  nicp = (lastknotu-3*2)*(lastknotv-3*2);
  nvars = 3*nicp;
        /* compute the integrals in dirty squares */
  if ( dirty ) {
    for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ ) {
        if ( dirty[sqk] & (DIRTY_FUNC | DIRTY_GRAD) ) {
          g1bl_UFuncGradSQd ( nkn, qcoeff, Nitab,
                              lastknotu, lastknotv, pitch, cp,
                              tC, isq, jsq, 2, lastknotu-4, 2, lastknotv-4,
                              &ftab[sqk], &gtab[sqk*SQUAREGDS] );
          dirty[sqk] &= ~(DIRTY_FUNC | DIRTY_GRAD);
        }
      }
  }
  else {
    for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ )
        g1bl_UFuncGradSQd ( nkn, qcoeff, Nitab,
                            lastknotu, lastknotv, pitch, cp,
                            tC, isq, jsq, 2, lastknotu-4, 2, lastknotv-4,
                            &ftab[sqk], &gtab[sqk*SQUAREGDS] );
  }
        /* sum the integrals */
  f = 0.0;
  memset ( grad, 0, nvars*sizeof(double) );
  for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
    for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ ) {
      f += ftab[sqk];
      for ( i = 0, ip = isq;  i < 3;  i++, ip++ )
        if ( ip >= 2 && ip < lastknotu-4 )
          for ( j = 0, jp = jsq;  j < 3;  j++, jp++ )
            if ( jp >= 2 && jp < lastknotv-4 ) {
	      /*kp - wewnetrzny punkt kontrolny dla ktorego liczony jest gradient */
              kp = 3*((ip-2)*(lastknotv-6) + jp-2);
	      
              pkn_AddMatrixd ( 1, 3, 0, &grad[kp],
                               0, &gtab[sqk*SQUAREGDS + 3*(3*i+j)],
                               0, &grad[kp] );
            }
    }
  *func = f;
} /*g1bl_UFuncGradd*/

/* ///////////////////////////////////////////////////
htab - do obliczenia hesjanu //////////////////////
func - wartosc funkcjonalu
grad - gradient 3*liczba punktow wewnetrznych */
void g1bl_UFuncGradHessiand ( int nkn, const double *qcoeff, double *Nitab, 
                              double *Nijtab, double *Mijtab,
                              int lastknotu, int lastknotv,
                              int pitch, point3d *cp, char *dirty,
                              double tC, double *ftab, double *gtab, double *htab,
                              double *func, double *grad,
                              int hsize, const int *prof, double **hrows )
			      /*hsize - dlugosc tablicy
			      prof  - profil
			      hrows - kolejne wiersze macierzy o nieregularnej budowie symetrycznej
			      */
{
  double f;
  int    nicp, nvars;
  int    isq, jsq, sqk;
  int    ia, ja, ka, iap, jap, kap, ib, jb, kb, ibp, jbp, kbp;
  int    i, j, bl;

  pitch /= 3;             /* express the pitch in terms of points */
  nicp = (lastknotu-3*2)*(lastknotv-3*2);
  nvars = 3*nicp;
        /* compute the integrals in dirty squares */
  if ( dirty ) {
    /* iterowanie po kwadratach jednostkowych dziedziny Omega=[0,M-m]x[0,N-n]*/
    for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ ) {
        if ( dirty[sqk] & (DIRTY_FUNC | DIRTY_GRAD | DIRTY_HESS) ) {
          g1bl_UFuncGradHessianSQd ( nkn, qcoeff, Nitab, Nijtab, Mijtab,
                  lastknotu, lastknotv, pitch, cp, tC, isq, jsq,
                  2, lastknotu-4, 2, lastknotv - 4,
                  &ftab[sqk], &gtab[sqk*SQUAREGDS], &htab[sqk*SQUAREHDS]);
          dirty[sqk] &= ~(DIRTY_FUNC | DIRTY_GRAD | DIRTY_HESS);
        }
      }
  }
  else {
    for ( isq = sqk = 0;  isq < lastknotu-DEG2;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-DEG2;  jsq++, sqk++ )
        g1bl_UFuncGradHessianSQd ( nkn, qcoeff, Nitab, Nijtab, Mijtab,
                lastknotu, lastknotv, pitch, cp, tC, isq, jsq,
                2, lastknotu-4, 2, lastknotv-4,
                &ftab[sqk], &gtab[sqk*SQUAREGDS], &htab[sqk*SQUAREHDS]);
  }
        /* sum the integrals */
  f = 0.0;
  memset ( grad, 0, nvars*sizeof(double) );
  memset ( hrows[0], 0, hsize*sizeof(double) );
  /* iteracja po kwadratach dziedziny */
  for ( isq = sqk = 0;  isq < lastknotu-4;  isq++ )
    for ( jsq = 0;  jsq < lastknotv-4;  jsq++, sqk++ ) {
      /* value of functional */
      f += ftab[sqk];
      for ( ia = 0, iap = isq;  ia < 3;  ia++, iap++ )
        if ( iap >= 2 && iap < lastknotu-4 )
          for ( ja = 0, jap = jsq;  ja < 3;  ja++, jap++ )
            if ( jap >= 2 && jap < lastknotv-4 ) {
              ka = 3*ia+ja;
	      
	      /* control point for which the derivative is calculated*/
              kap = 3*((iap-2)*(lastknotv-6) + jap-2);
              pkn_AddMatrixd ( 1, 3, 0, &grad[kap],
                               0, &gtab[sqk*SQUAREGDS + 3*ka],
                               0, &grad[kap] );
			       
              for ( ib = 0, ibp = isq;  ib <= iap;  ib++, ibp++ )
                if ( ibp >= 2 && ibp < lastknotu-4 )
                  for ( jb = 0, jbp = jsq;  jb < 3;  jb++, jbp++ )
                    if ( jbp >= 2 && jbp < lastknotv-4 ) {
                      kb = 3*ib+jb;
                      if ( kb <= ka ) {
                        kbp = 3*((ibp-2)*(lastknotv-6) + jbp-2);
                        bl = sqk*SQUAREHDS + 9*pkn_SymMatIndex ( ka, kb );
                        if ( kb < ka ) {
                          for ( i = 0; i < 3; i++ )
                            for ( j = 0; j < 3; j++ )
                              hrows[kap+i][kbp+j] += htab[bl+3*i+j];
                        }
                        else if ( kbp == kap ) {
                          for ( i = 0; i < 3; i++ )
                            for ( j = 0; j <= i; j++ )
                              hrows[kap+i][kbp+j] += htab[bl+3*i+j];
                        }
                      }
                    }
            }
    }
  *func = f;
} /*g1bl_UFuncGradHessiand*/

