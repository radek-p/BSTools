
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <sys/times.h>
#include <unistd.h>

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
#include "msgpool.h"

/* ///////////////////////////////////////////////////////////////////////// */
static void _g2mbl_MLSFindCPNormal ( int d, int *vertnum, int *mtab,
                                     void *usrptr )
{
  mesh_ml_optdata *data;
  point3d         *mvcp, cp[9];
  vector3d        *mvcpn, v1, v2, v3;
  int             *nncpi, i;

  data  = (mesh_ml_optdata*)usrptr;
  mvcp  = data->mvcp;
  mvcpn = data->mvcpn;
  nncpi = data->nncpi;
  if ( nncpi[vertnum[4]] >= 0 ) {
    for ( i = 0; i < 9; i++ )
      cp[i] = mvcp[vertnum[i]];
    pkn_MatrixLinCombd ( 1, 9, 0, &cp[0].x, 1.0/3.0, 0, &cp[3].x, 2.0/3.0,
                         0, &cp[0].x );
    pkn_MatrixLinCombd ( 1, 9, 0, &cp[3].x, 2.0/3.0, 0, &cp[6].x, 1.0/3.0,
                         0, &cp[3].x );
    pkn_MatrixLinCombd ( 2, 3, 9, &cp[0].x, 1.0/3.0, 9, &cp[1].x, 2.0/3.0,
                         9, &cp[0].x );
    pkn_MatrixLinCombd ( 2, 3, 9, &cp[1].x, 2.0/3.0, 9, &cp[2].x, 1.0/3.0,
                         9, &cp[1].x );
    MidPoint3d ( &cp[0], &cp[1], &cp[2] );
    MidPoint3d ( &cp[3], &cp[4], &cp[5] );
    SubtractPoints3d ( &cp[2], &cp[5], &v1 );
    MidPoint3d ( &cp[0], &cp[3], &cp[2] );
    MidPoint3d ( &cp[1], &cp[4], &cp[5] );
    SubtractPoints3d ( &cp[2], &cp[5], &v2 );
    CrossProduct3d ( &v1, &v2, &v3 );
    NormalizeVector3d ( &v3 );
    mvcpn[vertnum[4]] = v3;
  }
} /*_g2mbl_MLSFindCPNormal*/

static void _g2mbl_MLSFindSpecialCPNormal ( int d, int k, int *vertnum, int *mtab,
                                            void *usrptr )
{
  mesh_ml_optdata *data;
  int             *nncpi, i, j;
  point3d         *mvcp;
  vector3d        *mvcpn, v0, v1, v2, v3;

  data  = (mesh_ml_optdata*)usrptr;
  mvcp  = data->mvcp;
  mvcpn = data->mvcpn;
  nncpi = data->nncpi;
  if ( nncpi[vertnum[0]] >= 0 ) {
    memset ( &v0, 0, sizeof(vector3d) );
    for ( i = 0, j = k-1;  i < k;  j = i++ ) {
      SubtractPoints3d ( &mvcp[vertnum[2*i+1]], &mvcp[vertnum[0]], &v1 );
      SubtractPoints3d ( &mvcp[vertnum[2*j+1]], &mvcp[vertnum[0]], &v2 );
      CrossProduct3d ( &v1, &v2, &v3 );
      AddVector3d ( &v0, &v3, &v0 );
    }
    NormalizeVector3d ( &v0 );
    mvcpn[vertnum[0]] = v0;
  }
} /*_g2mbl_MLSFindSpecialCPNormal*/

boolean _g2mbl_MLSFindCPNormalsd ( mesh_ml_optdata *d )
{
  void        *sp;
  int         nv, nhe, nfac, *mvhei, *mfhei;
  BSMvertex   *mv;
  BSMhalfedge *mhe;
  BSMfacet    *mfac;

  sp = pkv_GetScratchMemTop ();
  nv = d->nv;
  nhe = d->nhe;
  nfac = d->nfac;
  mvhei = d->mvhei;
  mfhei = d->mfhei;
  mv = d->mv;
  mhe = d->mhe;
  mfac = d->mfac;

  PKV_MALLOC ( d->mvcpn, nv*sizeof(vector3d) );
  if ( !d->mvcpn ) {
    PKV_SIGNALERROR ( LIB_G2BLENDING, ERRCODE_1, ERRMSG_1 );
    goto failure;
  }
  memset ( d->mvcpn, 0, nv*sizeof(vector3d) );
  bsm_FindRegularSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                           3, d, _g2mbl_MLSFindCPNormal );
  bsm_FindSpecialVSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                            1, d, _g2mbl_MLSFindSpecialCPNormal );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLSFindCPNormalsd*/

