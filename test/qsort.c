
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "pkvaria.h"

#define N 40
int a[N] = {  1,  9, 10,  2,  0,  3,  6, 17, 19, 18,
             20, 29, 21, 28, 13, 34, 36,  5, 32,  7,
             30, 39, 31, 38, 23, 14, 16, 35,  8, 37,
             12, 11, 15,  4, 22, 27, 33, 26, 24, 25 };

int cnt = 0;

static boolean Less ( int i, int j, void *usrptr )
{
  cnt ++;
  return (boolean)(a[i] < a[j]);
} /*Less*/

static void Swap ( int i, int j, void *usrptr )
{
  int b;

  b = a[i];  a[i] = a[j];  a[j] = b;
} /*Swap*/

int main ( void )
{
  int i;

  pkv_QuickSort ( N, NULL, Less, Swap );
  for ( i = 0; i < N; i++ )
    printf ( " %3d", a[i] );
  printf ( "\n" );
  printf ( "cnt=%d\n", cnt );
  exit ( 0 );
} /*main*/

