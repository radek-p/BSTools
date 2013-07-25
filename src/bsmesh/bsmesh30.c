
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"
#include "pknum.h"
#include "bsmesh.h"
#include "bsmprivate.h"

static void _bsm_CountSpecialFSubnets ( int d, int k, int *vertnum, int *mtab,
                                        void *usrptr )
{
  index2 *cnt;

  cnt = (index2*)usrptr;
  cnt->i ++;
  cnt->j += k*(d+1)*(d+1);
} /*_bsm_CountSpecialFSubnets*/

boolean bsm_CountSpecialFSubnets ( int nv, BSMvertex *mv, int *mvhei,
                                   int nhe, BSMhalfedge *mhe,
                                   int nfac, BSMfacet *mfac, int *mfhei,
                                   byte snet_rad,
                                   int *nspecials, int *nspvert )
{
  index2 cnt;

  cnt.i = cnt.j = 0;
  if ( !bsm_FindSpecialFSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                  snet_rad, &cnt, _bsm_CountSpecialFSubnets ) )
    return false;
  *nspecials = cnt.i;
  *nspvert = cnt.j;
  return true;
} /*bsm_CountSpecialFSubnets*/

static void _bsm_StoreSpecialFSubnets ( int d, int k, int *vertnum, int *mtab,
                                        void *usrptr )
{
  bsm_special_elem_list *list;
  bsm_special_el        *spel;
  int                   *spvert;
  int                   nsp, nspv, n;

  list = (bsm_special_elem_list*)usrptr;
  spel = list->spel;
  spvert = list->spvert;
  nsp = list->nspecials;
  nspv = list->nspvert;
  spel[nsp].el_type = 1;  /* it's a special facet */
  spel[nsp].snet_rad = d;
  spel[nsp].degree = k;
  spel[nsp].snet_nvert = n = k*(d+1)*(d+1);
  spel[nsp].first_snet_vertex = nspv;
  memcpy ( &spvert[nspv], vertnum, n*sizeof(int) );
  list->nspecials ++;
  list->nspvert += n;
} /*_bsm_StoreSpecialFSubnets*/

boolean bsm_FindSpecialFSubnetLists ( int nv, BSMvertex *mv, int *mvhei,
                                      int nhe, BSMhalfedge *mhe,
                                      int nfac, BSMfacet *mfac, int *mfhei,
                                      byte snet_rad,
                                      boolean append,
                                      bsm_special_elem_list *list )
{
  if ( list->spel && list->spvert ) {
    if ( !append )
      list->nspecials = list->nspvert = 0;
    if ( bsm_FindSpecialFSubnets ( nv, mv, mvhei, nhe, mhe, nfac, mfac, mfhei,
                                    snet_rad, list, _bsm_StoreSpecialFSubnets ) )
      return true;
  }
  return false;
} /*bsm_FindSpecialFSubnetLists*/

