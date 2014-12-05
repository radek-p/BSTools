
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2014                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <math.h>

#include "pkvaria.h"

pkv_queue *pkv_InitQueue ( int nmax, int itemsize )
{
  pkv_queue *q;
  int       size;

  size = (nmax+1)*itemsize+sizeof(pkv_queue) - 1;
  PKV_MALLOC ( q, size );
  if ( q ) {
    q->nmax = nmax;
    q->itemsize = itemsize;
    q->fr = q->en = 0;
  }
  return q;
} /*pkv_InitQueue*/

pkv_queue *pkv_InitScratchQueue ( int nmax, int itemsize )
{
  pkv_queue *q;
  int       size;

  size = (nmax+1)*itemsize+sizeof(pkv_queue) - 1;
  q = pkv_GetScratchMem ( size );
  if ( q ) {
    q->nmax = nmax;
    q->itemsize = itemsize;
    q->fr = q->en = 0;
  }
  return q;
} /*pkv_InitScratchQueue*/

void pkv_ResetQueue ( pkv_queue *q )
{
  q->fr = q->en = 0;
} /*pkv_ResetQueue*/

boolean pkv_QueueEmpty ( pkv_queue *q )
{
  return q->fr == q->en;
} /*pkv_QueueEmpty*/

boolean pkv_QueueFull ( pkv_queue *q )
{
  if ( q->fr == 0 )
    return q->en == q->nmax;
  else
    return q->fr == q->en+1;
} /*pkv_QueueFull*/

boolean pkv_QueueInsert ( pkv_queue *q, void *item )
{
  if ( !pkv_QueueFull ( q ) ) {
    memcpy ( &q->qtab[q->en*q->itemsize], item, q->itemsize );
    q->en ++;
    if ( q->en > q->nmax )
      q->en = 0;
    return true;
  }
  else
    return false;
} /*pkv_QueueInsert*/

boolean pkv_QueueGetFirst ( pkv_queue *q, void *item )
{
  if ( !pkv_QueueEmpty ( q ) ) {
    memcpy ( item, &q->qtab[q->fr*q->itemsize], q->itemsize );
    return true;
  }
  else
    return false;
} /*pkv_QueueGetFirst*/

boolean pkv_QueueRemoveFirst ( pkv_queue *q, void *item )
{
  if ( !pkv_QueueEmpty ( q ) ) {
    if ( item )
      memcpy ( item, &q->qtab[q->fr*q->itemsize], q->itemsize );
    q->fr ++;
    if ( q->fr > q->nmax )
      q->fr = 0;
    return true;
  }
  else
    return false;
} /*pkv_QueueRemoveFirst*/

