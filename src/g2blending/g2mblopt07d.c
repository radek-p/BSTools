
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
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

boolean _g2mbl_MultAxd ( int n, void *usrdata, const double *x, double *Ax )
{
  mesh_lmt_optdata *md;
  nzHbl_rowdesc    *iHbl;
  int              *cHbl;
  double           *Hbl, *bHbl;
  int              i, i3, b, b9, b0, b1, j3, nvcp;

  md = (mesh_lmt_optdata*)usrdata;
  nvcp = md->nvcp;
  iHbl = md->iHbl;
  cHbl = md->cHbl;
  Hbl = md->Hbl;
  memset ( Ax, 0, n*sizeof(double) );
  for ( i = i3 = 0;  i < nvcp;  i++, i3 += 3 ) {
    b0 = iHbl[i].firsthbl;
    b1 = b0 + iHbl[i].nhbl - 1;
    for ( b = b0, b9 = 9*b;  b < b1;  b++, b9 += 9 ) {
      j3 = 3*cHbl[b];
      bHbl = &Hbl[b9];
        /* for speed no procedure calls, explicit instructions instead */
      Ax[i3]   += bHbl[0]*x[j3] + bHbl[1]*x[j3+1] + bHbl[2]*x[j3+2];
      Ax[i3+1] += bHbl[3]*x[j3] + bHbl[4]*x[j3+1] + bHbl[5]*x[j3+2];
      Ax[i3+2] += bHbl[6]*x[j3] + bHbl[7]*x[j3+1] + bHbl[8]*x[j3+2];

      Ax[j3]   += bHbl[0]*x[i3] + bHbl[3]*x[i3+1] + bHbl[6]*x[i3+2];
      Ax[j3+1] += bHbl[1]*x[i3] + bHbl[4]*x[i3+1] + bHbl[7]*x[i3+2];
      Ax[j3+2] += bHbl[2]*x[i3] + bHbl[5]*x[i3+1] + bHbl[8]*x[i3+2];
    }
    bHbl = &Hbl[b9];
    Ax[i3]   += bHbl[0]*x[i3] + bHbl[1]*x[i3+1] + bHbl[2]*x[i3+2];
    Ax[i3+1] += bHbl[3]*x[i3] + bHbl[4]*x[i3+1] + bHbl[5]*x[i3+2];
    Ax[i3+2] += bHbl[6]*x[i3] + bHbl[7]*x[i3+1] + bHbl[8]*x[i3+2];
  }
  return true;
} /*_g2mbl_MultAxd*/

