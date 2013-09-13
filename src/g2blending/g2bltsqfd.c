
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
/* The procedures in this file evaluate the functionals measuring the        */
/* surface badness; the first one is auxiliary.                              */
/* ///////////////////////////////////////////////////////////////////////// */
void g2bl_FuncTSQFsqd ( int nkn, const double *qcoeff, double *Nitab,
                        int lastknotu, int lastknotv,
                        int pitch, point3d *cp,
                        double tC,
                        int isq, int jsq,
                        double *sT, double *sS, double *sQ, double *sF )
{
  int      fcpn; /* first control point number */
  vector3d pder[11];                       /* parameterization derivatives */
  double   Gstar[9], detG, detGu, detGv;  /* first form and related data */
  double   Bstar[11], tH, tHu, tHv;       /* second form and related data */
  double   L[2], M[3], N, f, g;
  int      i, j;
  double   sum0T, sum1T, sum0S, sum1S, sum0Q, sum1Q;

  fcpn = isq*pitch + jsq;
  sum0T = sum0S = sum0Q = 0.0;
  for ( i = 0; i < nkn; i++ ) {
    sum1T = sum1S = sum1Q = 0.0;
    for ( j = 0; j < nkn; j++ ) {
        /* compute the parameterization derivatives */
      _g2bl_UCompPDerd ( nkn, Nitab, pitch, cp, fcpn, i, j, pder );
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
      tHu = g11u*tb22 + g11*tb22u + tb11u*g22 + tb11*g22u -
            2.0*(g12u*tb12 + g12*tb12u);
      tHv = g11v*tb22 + g11*tb22v + tb11v*g22 + tb11*g22v -
            2.0*(g12v*tb12 + g12*tb12v);
      L[0] = detG*tHu - 1.5*detGu*tH;
      L[1] = detG*tHv - 1.5*detGv*tH;

      f = 0.25*(L[0]*(L[0]*M[0]+L[1]*(M[1]+M[1]))+L[1]*L[1]*M[2])*N;
      sum1S += f*qcoeff[j];

        /* the second functional */
          /* Frobenius norm of the Laplacian gradient */
      g = DotProduct3d ( &pder[9], &pder[9] ) +
          DotProduct3d ( &pder[10], &pder[10] );
      sum1T += g*qcoeff[j];
          /* Frobenius norm of the projection of the Laplacian gradient */
          /* on the normal vector */
      g -= (Bstar[9]*Bstar[9] + Bstar[10]*Bstar[10])/detG;
      sum1Q += g*qcoeff[j];
    }
    sum0T += sum1T*qcoeff[i];
    sum0S += sum1S*qcoeff[i];
    sum0Q += sum1Q*qcoeff[i];
  }
        /* store the integral over the square in the array */
  *sT = sum0T;  *sS = sum0S;  *sQ = sum0Q;
  *sF = sum0S + tC*sum0Q;
} /*g2bl_FuncTSQFsqd*/

boolean g2bl_FuncTSQFd ( int nkn,
                         int lastknotu, int lastknotv, int pitch, point3d *cp,
                         double tC,
                         double *fT, double *fS, double *fQ, double *fF )
{
  void   *sp;
  double *qknots, *qcoeff, *Nitab;
  double *bf, *dbf, *ddbf, *dddbf;
  int    isq, jsq, sqk, size;
  double sT, sS, sQ, sF;

  sp = pkv_GetScratchMemTop ();
  pitch /= 3;     /* express the pitch in terms of points */

  size = g2bl_NiSize ( nkn );
  Nitab = pkv_GetScratchMemd ( size );
  if ( !Nitab )
    goto failure;
  if ( !_g2bl_TabBasisFuncd ( nkn, &qknots, &qcoeff, &bf, &dbf, &ddbf, &dddbf ) )
    goto failure;
  g2bl_TabNid ( nkn, bf, dbf, ddbf, dddbf, Nitab );
  *fT = *fS = *fQ = *fF = 0;
  for ( isq = sqk = 0;  isq < lastknotu-6;  isq++ )
    for ( jsq = 0;  jsq < lastknotv-6;  jsq++, sqk++ ) {
      g2bl_FuncTSQFsqd ( nkn, qcoeff, Nitab, lastknotu, lastknotv, pitch, cp,
                         tC, isq, jsq, &sT, &sS, &sQ, &sF );
      *fT += sT;  *fS += sS;  *fQ += sQ;  *fF += sF;
    }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*g2bl_UFuncd*/

