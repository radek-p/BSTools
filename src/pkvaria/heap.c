
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2009                                  */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>

#include "pkvaria.h"

/* ///////////////////////////////////////////////////////////////////////// */
/* heap priority queue implementation; the array is used to store pointers,  */
/* and an application-provided procedure compares the element priorities.    */
/* This procedure returns true if the object pointed by the first parameter  */
/* has a greater priority. */

/* pkn_UpHeap - the array contains l+1 pointers, all except for the last one */
/* point to heap elements ordered correctly */
int pkv_UpHeap ( void *a[], int l, boolean (*cmp)(void*,void*) )
{
  int  i;
  void *b;

  while ( l > 0 ) {
    i = (l-1) / 2;
    if ( cmp ( a[l], a[i] ) ) {
      b = a[l];  a[l] = a[i];  a[i] = b;
      l = i;
    }
    else
      break;
  }
  return l;  /* return the final position of the initially l-th element */
} /*pkv_UpHeap*/

/* pkn_DownHeap - the array contains l+1 pointers, the f-th points to an */
/* object, which might be at the wrong place */
int pkv_DownHeap ( void *a[], int l, int f, boolean (*cmp)(void*,void*) )
{
  int  i, j;
  void *b;

  i = f+f+1;
  while ( i <= l ) {
    j = i+1;
    if ( j <= l )
      if ( cmp ( a[j], a[i] ) ) i = j;
    if ( cmp ( a[i], a[f] ) ) {
      b = a[f];  a[f] = a[i];  a[i] = b;
      f = i;
      i = f+f+1;
    }
    else break;
  }
  return f;
} /*pkv_DownHeap*/

int pkv_HeapInsert ( void *a[], int *l, void *newelem,
                     boolean (*cmp)(void*,void*) )
{
  int ll;

  *l = ll = *l + 1;
  a[ll] = newelem;
  return pkv_UpHeap ( a, ll, cmp );
} /*pkv_HeapInsert*/

void pkv_HeapRemove ( void *a[], int *l, int el,
                      boolean (*cmp)(void*,void*) )
{
  int ll;

  ll = *l;
  if ( el >= 0 && el < ll ) {
    a[el] = a[ll--];
    *l = ll;
    pkv_DownHeap ( a, ll, el, cmp );
  }
  else
    *l = ll-1;
} /*pkv_HeapRemove*/

void pkv_HeapOrder ( void *a[], int n, boolean (*cmp)(void*,void*) )
{
  int i;

  for ( i = n/2-1; i >= 0; i-- )
    pkv_DownHeap ( a, n-1, i, cmp );
} /*pkv_HeapOrder*/

void pkv_HeapSort ( void *a[], int n, boolean (*cmp)(void*,void*) )
{
  int  i;
  void *b;

  for ( i = n/2-1; i >= 0; i-- )
    pkv_DownHeap ( a, n-1, i, cmp );
  for ( i = n-1; i > 0; i-- ) {
    b = a[0];  a[0] = a[i];  a[i] = b;
    pkv_DownHeap ( a, i-1, 0, cmp );
  }
} /*pkv_HeapSort*/

