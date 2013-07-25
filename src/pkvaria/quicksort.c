
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2011                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <math.h>
#include <string.h>

#include "pkvaria.h"

/* /////////////////////////////////////////// */
/* QuickSort - this implementation does not process the sequence */
/* directly (by referring to an array), but it rather uses the   */
/* application-supplied procedures, less to compare two elements */
/* and swap to sw*/

static void rQuickSort ( int i, int j, void *usrptr,
                         boolean (*less)( int i, int j, void *usrptr ),
                         void (*swap)( int i, int j, void *usrptr ) )
{
  int p, q;

  while ( j-i > 8 ) {
        /* finding the pivot */
    p = (i+j)/2;
        /* partition */
    swap ( i, p, usrptr );
    p = i;
    for ( q = i+1; q <= j; q++ )
      if ( less ( q, i, usrptr ) )
        swap ( ++p, q, usrptr );
    swap ( i, p, usrptr );
        /* one recursive call, followed by repetition in the while loop. */
        /* the if instruction ensures the smallest recursion depth. */
    if ( j-p < p-i ) {
      rQuickSort ( p+1, j, usrptr, less, swap );
      j = p-1;
    }
    else {
      rQuickSort ( i, p-1, usrptr, less, swap );
      i = p+1;
    }
  }
} /*rQuickSort*/

void pkv_QuickSort ( int n, void *usrptr,
                     boolean (*less) ( int i, int j, void *usrptr ),
                     void (*swap) ( int i, int j, void *usrptr ) )
{
  int i, j, k;

  if ( n > 1 ) {
        /* use the true QuickSort first */
    rQuickSort ( 0, n-1, usrptr, less, swap );
        /* then finish the work by using InsertionSort */
    for ( i = 1; i < n; i++ ) {
      k = i;  j = k-1;
      while ( j >= 0 )
        if ( less ( k, j, usrptr ) ) {
          swap ( j, k, usrptr );
          k = j--;
        }
        else break;
    }
  }
} /*pkv_QuickSort*/

