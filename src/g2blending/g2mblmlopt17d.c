
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2012                                  */
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

#define _DEBUG
#define NEW_METHOD

/* ////////////////////////////////////////////////////////////////////////// */
/* The old algorithm is as follows: make an iteration for a block and get to  */
/* the next one if advance is true, which means that the Newton method step   */
/* for this block was successful.                                             */
/* The new algorithm is the following: make one iteration for each block at   */
/* the current level (changing them in a cyclic way), until the Newton method */
/* for all these blocks is successful. Then advance to the next level.        */
/* The new algorithm seems to work better; actually the old method needed     */
/* very good starting points for the optimization of the entire mesh.         */
/* ////////////////////////////////////////////////////////////////////////// */

#ifdef NEW_METHOD
static short _get_level ( short bl )
{
  short lev;

  if ( bl < 0 )
    return -1;
  bl = (bl+1) >> 1;
  lev = 0;
  while ( bl ) {
    bl >>= 1;
    lev ++;
  }
  return lev;
} /*_get_level*/
#endif

int _g2mbl_MLNextBlockNumd ( mesh_ml_optdata *d, boolean advance )
{
#ifdef NEW_METHOD
  short currentblock, nextblock, currentlevel, nextlevel;

  currentblock = d->currentblock;
  if ( currentblock == 0 && !advance )
    return 0;
  nextblock = currentblock - 1;
  currentlevel = _get_level ( currentblock );
  nextlevel = _get_level ( nextblock );
  d->nextlevel = min ( nextlevel, d->nextlevel );
  if ( nextblock >= 0 && currentblock < d->dirtyblock ) {
    if ( nextlevel < currentlevel )
      nextblock = (0x01 << (currentlevel+1)) - 2;
  }
  else if ( advance ) {
    if ( !currentblock )
      d->dirtyblock = -1;
    else if ( currentblock == d->dirtyblock ) {
      d->dirtyblock = 0;
      if ( d->nextlevel < currentlevel )
        nextblock = (0x01 << (d->nextlevel+1)) - 2;
    }
  }
  if ( !advance ) {
#ifdef _DEBUG
printf ( " <-" );
#endif
    d->dirtyblock = currentblock;
    if ( nextlevel < currentlevel )
      nextblock = 2*currentblock;
  }
  return (d->currentblock = nextblock);
#else
  if ( advance )
    d->currentblock --;
  return d->currentblock;
#endif
} /*_g2mbl_MLNextBlockNumd*/

