
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2012                            */
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
#include "g2mblmlprivated.h"

/* ///////////////////////////////////////////////////////////////////////// */
boolean _g2mbl_MLSetupElemConstd ( mesh_ml_optdata *d,
                                   double dM, double dO, double C )
{
  void         *sp;
  int          nv, ndomelems;
  point3d      *mvcp;
  meshdom_elem *domelem;
  int          *nncpi, *vncpi, nvcp;
  int          i, j;
  double       r;

  sp = pkv_GetScratchMemTop ();
  nv        = d->nv;
  mvcp      = d->mvcp;
  ndomelems = d->ndomelems;
  domelem   = d->domelem;
        /* set up the constant c for each element */
  if ( dM <= 0.0 ) {
          /* compute the "surface diameter" */
    nncpi = pkv_GetScratchMemi ( nv );
    if ( !nncpi )
      goto failure;
    nvcp = d->nvcp;
    vncpi = d->bd[0].vncpi;
    for ( i = j = 0;  i < nv;  i++ )
      if ( j < nvcp && vncpi[j] == i )
        nncpi[i] = j++;
      else
        nncpi[i] = -1;
    for ( i = 1; i < nv; i++ )
      if ( nncpi[i] == -1 )
        for ( j = 0; j < i; j++ )
          if ( nncpi[j] == -1 ) {
            r = Point3Distanced ( &mvcp[i], &mvcp[j] );
            dM = max ( dM, r );
          }
    if ( dM <= 0 )
      goto failure;
    pkv_SetScratchMemTop ( sp );  /* deallocate nncpi */
  }
  dM *= dM;

  if ( dO <= 0.0 ) {
          /* compute the "domain diameter" */
    j = 0;
    for ( i = 0; i < ndomelems; i++ )
      if ( domelem[i].type == 4 )
        j ++;
      else
        j += domelem[i].type;
    dO = (double)j;
  }
  else
    dO *= dO;
          /* assign the  constant */
  r = C*dO*dO/(dM*dM);
  for ( i = 0; i < ndomelems; i++ )
    domelem[i].C = r;

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSetupElemConstd*/


