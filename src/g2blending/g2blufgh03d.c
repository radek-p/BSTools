
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
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
void g2bl_ClosedUFuncGradd ( int nkn, const double *qcoeff, double *Nitab,
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
  nicp = (lastknotu-6)*(lastknotv-9);
  nvars = 3*nicp;
        /* compute the integrals in dirty squares */
  if ( dirty ) {
    for ( isq = sqk = 0;  isq < lastknotu-6;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-6;  jsq++, sqk++ ) {
        if ( dirty[sqk] & (DIRTY_FUNC | DIRTY_GRAD) ) {
          g2bl_UFuncGradSQd ( nkn, qcoeff, Nitab,
                              lastknotu, lastknotv, pitch, cp,
                              tC, isq, jsq, 0, lastknotu-3, 3, lastknotv-6,
                              &ftab[sqk], &gtab[sqk*SQUAREGDS] );
          dirty[sqk] &= ~(DIRTY_FUNC | DIRTY_GRAD);
        }
      }
  }
  else {
    for ( isq = sqk = 0;  isq < lastknotu-6;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-6;  jsq++, sqk++ )
        g2bl_UFuncGradSQd ( nkn, qcoeff, Nitab,
                            lastknotu, lastknotv, pitch, cp,
                            tC, isq, jsq, 0, lastknotu-3, 3, lastknotv-6,
                            &ftab[sqk], &gtab[sqk*SQUAREGDS] );
  }
        /* sum the integrals */
  f = 0.0;
  memset ( grad, 0, nvars*sizeof(double) );
  for ( isq = sqk = 0;  isq < lastknotu-6;  isq++ )
    for ( jsq = 0;  jsq < lastknotv-6;  jsq++, sqk++ ) {
      f += ftab[sqk];
      for ( i = 0, ip = isq;  i < 4;  i++, ip++ )
        if ( ip < lastknotu-6 ) {
          for ( j = 0, jp = jsq;  j < 4;  j++, jp++ )
            if ( jp >= 3 && jp < lastknotv-6 ) {
              kp = 3*(ip*(lastknotv-9) + jp-3);
              pkn_AddMatrixd ( 1, 3, 0, &grad[kp],
                               0, &gtab[sqk*SQUAREGDS + 3*(4*i+j)],
                               0, &grad[kp] );
            }
        }
        else {
          for ( j = 0, jp = jsq;  j < 4;  j++, jp++ )
            if ( jp >= 3 && jp < lastknotv-6 ) {
              kp = 3*((ip-lastknotu+6)*(lastknotv-9) + jp-3);
              pkn_AddMatrixd ( 1, 3, 0, &grad[kp],
                               0, &gtab[sqk*SQUAREGDS + 3*(4*i+j)],
                               0, &grad[kp] );
            }
        }
    }
  *func = f;
} /*g2bl_ClosedUFuncGradd*/

/* ///////////////////////////////////////////////////////////////////////// */
void g2bl_ClosedUFuncGradHessiand ( int nkn, const double *qcoeff, double *Nitab, 
                              double *Nijtab, double *Mijtab,
                              int lastknotu, int lastknotv,
                              int pitch, point3d *cp, char *dirty,
                              double tC, double *ftab, double *gtab, double *htab,
                              double *func, double *grad,
                              int hsize, const int *prof, double **hrows )
{
  double f;
  int    nicp, nvars;
  int    isq, jsq, sqk;
  int    ia, ja, ka, iap, jap, kap, ib, jb, kb, ibp, jbp, kbp;
  int    i, j, bl;

  pitch /= 3;             /* express the pitch in terms of points */
  nicp = (lastknotu-6)*(lastknotv-9);
  nvars = 3*nicp;
        /* compute the integrals in dirty squares */
  if ( dirty ) {
    for ( isq = sqk = 0;  isq < lastknotu-6;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-6;  jsq++, sqk++ ) {
        if ( dirty[sqk] & (DIRTY_FUNC | DIRTY_GRAD | DIRTY_HESS) ) {
#ifdef __DEBUG
printf ( "%5d\b\b\b\b\b", sqk );
#endif
          g2bl_UFuncGradHessianSQd ( nkn, qcoeff, Nitab, Nijtab, Mijtab,
                  lastknotu, lastknotv, pitch, cp, tC, isq, jsq,
                  0, lastknotu-3, 3, lastknotv-6,
                  &ftab[sqk], &gtab[sqk*SQUAREGDS], &htab[sqk*SQUAREHDS] );
          dirty[sqk] &= ~(DIRTY_FUNC | DIRTY_GRAD | DIRTY_HESS);
        }
      }
  }
  else {
    for ( isq = sqk = 0;  isq < lastknotu-6;  isq++ )
      for ( jsq = 0;  jsq < lastknotv-6;  jsq++, sqk++ ) {
#ifdef __DEBUG
printf ( "%5d\b\b\b\b\b", sqk );
#endif
        g2bl_UFuncGradHessianSQd ( nkn, qcoeff, Nitab, Nijtab, Mijtab,
                lastknotu, lastknotv, pitch, cp, tC, isq, jsq,
                0, lastknotu-3, 3, lastknotv-6,
                &ftab[sqk], &gtab[sqk*SQUAREGDS], &htab[sqk*SQUAREHDS] );
      }
  }
        /* sum the integrals */
  f = 0.0;
  memset ( grad, 0, nvars*sizeof(double) );
  memset ( hrows[0], 0, hsize*sizeof(double) );
  for ( isq = sqk = 0;  isq < lastknotu-6;  isq++ )
    for ( jsq = 0;  jsq < lastknotv-6;  jsq++, sqk++ ) {
      f += ftab[sqk];
      for ( ia = 0, iap = isq;  ia < 4;  ia++, iap++ )
        if ( iap < lastknotu-6 ) {
          for ( ja = 0, jap = jsq;  ja < 4;  ja++, jap++ )
            if ( jap >= 3 && jap < lastknotv-6 ) {
              ka = 4*ia+ja;
              kap = 3*(iap*(lastknotv-9) + jap-3);
              pkn_AddMatrixd ( 1, 3, 0, &grad[kap],
                               0, &gtab[sqk*SQUAREGDS + 3*(4*ia+ja)],
                               0, &grad[kap] );
              for ( ib = 0, ibp = isq;  ib <= iap;  ib++, ibp++ )
                if ( ibp < lastknotu-6 )
                  for ( jb = 0, jbp = jsq;  jb < 4;  jb++, jbp++ )
                    if ( jbp >= 3 && jbp < lastknotv-6 ) {
                      kb = 4*ib+jb;
                      if ( kb <= ka ) {
                        kbp = 3*(ibp*(lastknotv-9) + jbp-3);
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
        else {
          for ( ja = 0, jap = jsq;  ja < 4;  ja++, jap++ )
            if ( jap >= 3 && jap < lastknotv-6 ) {
              ka = 4*ia+ja;
              kap = 3*((iap-lastknotu+6)*(lastknotv-9) + jap-3);
              pkn_AddMatrixd ( 1, 3, 0, &grad[kap],
                               0, &gtab[sqk*SQUAREGDS + 3*(4*ia+ja)],
                               0, &grad[kap] );
              for ( ib = 0, ibp = isq;  ib <= iap;  ib++, ibp++ )
                if ( ibp < lastknotu-6 ) {
                  for ( jb = 0, jbp = jsq;  jb < 4;  jb++, jbp++ )
                    if ( jbp >= 3 && jbp < lastknotv-6 ) {
                      kb = 4*ib+jb;
                      kbp = 3*(ibp*(lastknotv-9) + jbp-3);
                      bl = sqk*SQUAREHDS + 9*pkn_SymMatIndex ( ka, kb );
                      for ( i = 0; i < 3; i++ )
                        for ( j = 0; j < 3; j++ )
                          hrows[kbp+i][kap+j] += htab[bl+3*j+i];
                    }
                }
                else {
                  for ( jb = 0, jbp = jsq;  jb < 4;  jb++, jbp++ )
                    if ( jbp >= 3 && jbp < lastknotv-6 ) {
                      kb = 4*ib+jb;
                      if ( kb <= ka ) {
                        kbp = 3*((ibp-lastknotu+6)*(lastknotv-9) + jbp-3);
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
        }
    }
  *func = f;
} /*g2bl_ClosedUFuncGradHessiand*/

