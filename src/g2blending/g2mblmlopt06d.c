
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2011                                  */
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
/* the two procedures below are used to count the mesh facets corresponding to */
/* ordinary elements and Sabin nets corresponding to special elements */
static void _g2mbl_MLCountRegularFacets ( int d, int *vertnum, int *mtab,
                                          void *usrptr )
{
  mesh_ml_optdata *data;
  int             i;
  int             *nncpi;

  data = (mesh_ml_optdata*)usrptr;
  data->tdomelems ++;
  nncpi = data->nncpi;
  for ( i = 0; i < 16; i++ )
    if ( nncpi[vertnum[i]] >= 0 )
      goto cont;
  return;

cont:
  data->ndomelems ++;
  data->ndomelcpind += 16;
  data->eltypes[1] = true;
} /*_g2mbl_MLCountRegularFacets*/

static void _g2mbl_MLCountSpecialVertices ( int d, int k, int *vertnum, int *mtab,
                                            void *usrptr )
{
  mesh_ml_optdata *data;
  int             i, ncpi;
  int             *nncpi;

  data = (mesh_ml_optdata*)usrptr;
  data->tdomelems ++;
  nncpi = data->nncpi;
  ncpi = 6*k+1;
  for ( i = 0; i < ncpi; i++ )
    if ( nncpi[vertnum[i]] >= 0 )
      goto cont;
  return;

cont:
  data->ndomelems ++;
  data->ndomelcpind += ncpi;
  data->eltypes[k-3] = true;
} /*_g2mbl_MLCountSpecialVertices*/

/* the two procedures below store the numbers of vertices (i.e. of basis */
/* functions) nonzero in ordinary and special elements */
static void _g2mbl_MLGetRegularElemVertNum ( int d, int *vertnum, int *mtab,
                                             void *usrptr )
{
  mesh_ml_optdata *data;
  int             nel, fcpi, i, bls;
  int             *nncpi;
  meshdom_elem    *domelem;

  data = (mesh_ml_optdata*)usrptr;
  nncpi = data->nncpi;
  for ( i = 0; i < 16; i++ )
    if ( nncpi[vertnum[i]] >= 0 )
      goto cont;
  return;

cont:
  nel = data->ndomelems ++;
  fcpi = data->ndomelcpind;
  data->ndomelcpind += 16;
  domelem = &data->domelem[nel];
  domelem->type = 4;  /* ordinary */
  domelem->ncp = 16;
  domelem->firstcpi = fcpi;
  domelem->hti = data->hti;
  domelem->cfn = mtab[24];
  memcpy ( &data->domelcpind[fcpi], vertnum, 16*sizeof(int) );
  bls = data->bls;
  data->hti += bls*bls*136;
} /*_g2mbl_MLGetRegularElemVertNum*/

static void _g2mbl_MLGetSpecialElemVertNum ( int d, int k, int *vertnum, int *mtab,
                                             void *usrptr )
{
  mesh_ml_optdata *data;
  int             nel, fcpi, ncpi, i, bls;
  int             *nncpi;
  meshdom_elem    *domelem;

  data = (mesh_ml_optdata*)usrptr;
  nncpi = data->nncpi;
  ncpi = 6*k+1;
  for ( i = 0; i < ncpi; i++ )
    if ( nncpi[vertnum[i]] >= 0 )
      goto cont;
  return;

cont:
  nel = data->ndomelems ++;
  fcpi = data->ndomelcpind;
  data->ndomelcpind += ncpi;
  domelem = &data->domelem[nel];
  domelem->type = k;    /* 3 or 5,...,16 */
  domelem->ncp = ncpi;
  domelem->firstcpi = fcpi;
  domelem->hti = data->hti;
  domelem->cfn = mtab[6];
  memcpy ( &data->domelcpind[fcpi], vertnum, ncpi*sizeof(int) );
  bls = data->bls;
  data->hti += bls*bls*(ncpi*(ncpi+1))/2;
} /*_g2mbl_MLGetSpecialElemVertNum*/

boolean _g2mbl_MLFindElementsd ( mesh_ml_optdata *d, int bls, boolean dnkn )
{
  void         *sp;
  int          nv, nhe, nfac, *mvhei, *mfhei;
  BSMvertex    *mv;
  BSMhalfedge  *mhe;
  BSMfacet     *mfac;
  int          ndomelems, ndomelcpind;
  int          *nncpi;
  int          nvcp, *vncpi;
  int          i, size;

  sp = pkv_GetScratchMemTop ();
        /* get the mesh */
  nv    = d->nv;
  mv    = d->mv;
  mvhei = d->mvhei;
  nhe   = d->nhe;
  mhe   = d->mhe;
  nfac  = d->nfac;
  mfac  = d->mfac;
  mfhei = d->mfhei;
        /* determine the variable vertices */
  PKV_MALLOC ( d->nncpi, nv*sizeof(int) );
  if ( !d->nncpi )
    goto failure;
  nncpi = d->nncpi;
  for ( i = 0; i < nv; i++ )
    nncpi[i] = -1;
  nvcp = d->bd[0].nvcp;
  vncpi = d->bd[0].vncpi;
  for ( i = 0; i < nvcp; i++ )
    nncpi[vncpi[i]] = i;
        /* count all elements */
  d->bls = bls;
  d->tdomelems = d->ndomelems = d->ndomelcpind = 0;
  bsm_FindRegularSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                           4, d, _g2mbl_MLCountRegularFacets );
  bsm_FindSpecialVSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                            2, d, _g2mbl_MLCountSpecialVertices );
  ndomelems   = d->ndomelems;
  ndomelcpind = d->ndomelcpind;
        /* allocate arrays for the elements */
  PKV_MALLOC ( d->domelem, ndomelems*sizeof(meshdom_elem) +
                           ndomelcpind*sizeof(int) );
  if ( !d->domelem )
    goto failure;
  d->domelcpind = (int*)&d->domelem[ndomelems];
  d->ndomelems = d->ndomelcpind = d->hti = 0;
  bsm_FindRegularSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                           4, d, _g2mbl_MLGetRegularElemVertNum );
  bsm_FindSpecialVSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                            2, d, _g2mbl_MLGetSpecialElemVertNum );
  size = ndomelems + bls*ndomelcpind + d->hti;
  if ( dnkn )
    size += ndomelems + bls*ndomelcpind;

printf ( "ftabsize = %d\n", size );

  PKV_MALLOC ( d->ftab1, size*sizeof(double) );
  if ( !d->ftab1 )
    goto failure;
  if ( dnkn ) {
    d->ftab2 = &d->ftab1[ndomelems];
    d->gtab1 = &d->ftab2[ndomelems];
    d->gtab2 = &d->gtab1[bls*ndomelcpind];
    d->htab1 = &d->gtab2[bls*ndomelcpind];
  }
  else {
    d->ftab2 = d->ftab1;
    d->gtab2 = d->gtab1 = &d->ftab1[ndomelems];
    d->htab1 = &d->gtab1[bls*ndomelcpind];
  }
  memset ( d->ftab1, 0, size*sizeof(double) );

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  d->domelem = NULL;
  pkv_SetScratchMemTop ( sp );
  return false;
} /*_g2mbl_MLFindElementsd*/

