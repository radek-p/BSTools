
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>

#include "pkvaria.h"

#include "proc_regmem.h"

/*#define DEBUG*/

static void *blocklist[MAX_MEM_BLOCKS];
static int nblocks = 0;

void MemBlockListInit ( void )
{
  nblocks = 0;
} /*MemBlockListInit*/

void MemBlockRegister ( void *ptr, boolean alloc )
{
  int i;

  if ( alloc ) {
#ifdef DEBUG
printf ( "alloc, ptr = %x, n = %d\n", (int)ptr, nblocks  );
#endif
    if ( nblocks < MAX_MEM_BLOCKS )
      blocklist[nblocks++] = ptr;
    else {
      printf ( "memory block register overflow\n" );
      exit ( 1 );
    }
  }
  else {
#ifdef DEBUG
printf ( "free, ptr = %x, n = %d\n", (int)ptr, nblocks  );
#endif
    for ( i = 0; i < nblocks; i++ )
      if ( ptr == blocklist[i] ) {
        blocklist[i] = blocklist[--nblocks];
        return;
      }
    printf ( "memory block register underflow or absent block\n" );
    exit ( 1 );
  }
} /*MemBlockRegister*/

void MemBlockFreeAll ( void )
{
  int i;
  boolean cr;

  cr = pkv_critical;
  pkv_critical = true;
  for ( i = 0; i < nblocks; i++ )
    free ( blocklist[i] );
  nblocks = 0;
  pkv_critical = cr;
} /*MemBlockFreeAll*/

