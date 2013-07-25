
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pkvaria.h"

void *DMalloc ( size_t size )
{
  char *ptr;

  ptr = malloc ( size+16 );
#if __WORDSIZE == 64
  printf ( "malloc ( %ld )", size );
#else
  printf ( "malloc ( %d )", size );
#endif
  memset ( ptr, 0, size+16 );
  *((size_t*)ptr) = size;
#if __WORDSIZE == 64
  printf ( " -> %lx\n", (size_t)ptr );
#else
  printf ( " -> %x\n", (size_t)ptr );
#endif
  return (void*)(ptr+8);
} /*DMalloc*/

void DFree ( void *ptr )
{
  int    i;
  size_t size;
  char   bad0, bad1, *p;

  p = (char*)ptr;
#if __WORDSIZE == 64
  p -= 16;
  size = *((size_t*)p);
  printf ( "free ( %lx ), size %ld", (size_t)p, size );
  bad0 = bad1 = 0;
  for ( i = 8; i < 15; i++ )
    if ( p[i] ) { bad0 = 1;  break; }
  for ( i = size+16; i < size+32; i++ )
    if ( p[i] ) { bad1 = 1;  break; }
  if ( bad0 )
    printf ( " bad0" );
  if ( bad1 )
    printf ( " bad1" );
#else
  p -= 8;
  size = *((size_t*)p);
  printf ( "free ( %x ), size %d", (size_t)p, size );
  bad0 = bad1 = 0;
  for ( i = 4; i < 7; i++ )
    if ( p[i] ) { bad0 = 1;  break; }
  for ( i = size+8; i < size+16; i++ )
    if ( p[i] ) { bad1 = 1;  break; }
  if ( bad0 )
    printf ( " bad0" );
  if ( bad1 )
    printf ( " bad1" );
#endif
  printf ( "\n" );
  free ( (void*)p );
} /*DFree*/
