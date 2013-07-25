
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkgeom.h"
#include "linkage.h"

kl_link *NewLink ( int id, void (*drawlink) ( struct kl_link *lnk, trans3d *tr ) )
{
  kl_link *lk;

  lk = malloc ( sizeof(kl_link) );
  if ( lk ) {
    lk->id = id;
    lk->drawlink = drawlink;
    lk->firstjoint = NULL;
    IdentTrans3d ( &lk->M );
  }
  return lk;
} /*NewLink*/

kl_joint *NewJoint ( int id, kl_link *first, kl_link *second,
                     void (*articulate) ( struct kl_joint *jnt ) )
{
  kl_joint *jt;

  jt = malloc ( sizeof(kl_joint) );
  if ( jt ) {
    jt->id = id;
    jt->up = first;
    if ( first ) {
      jt->next = first->firstjoint;
      first->firstjoint = jt;
    }
    else
      jt->next = NULL;
    jt->lnk = second;
    jt->articulate = articulate;
    IdentTrans3d ( &jt->L );
    IdentTrans3d ( &jt->R );
  }
  return jt;
} /*NewJoint*/

void DisplayLinkage ( kl_joint *jnt, trans3d *m )
{
  trans3d a, b;
  kl_link   *lnk;

  lnk = jnt->lnk;
  if ( lnk ) {
    if ( jnt->articulate )
      jnt->articulate ( jnt );
    CompTrans3d ( &b, m, &jnt->L );
    CompTrans3d ( &a, &b, &jnt->R  );
    CompTrans3d ( &b, &a, &lnk->M  );
    if ( lnk->drawlink )
      lnk->drawlink ( lnk, &b );
    for ( jnt = lnk->firstjoint;  jnt;  jnt = jnt->next ) {
      DisplayLinkage ( jnt, &a );
    }
  }
} /*DisplayLinkage*/

