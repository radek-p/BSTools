
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"

#define N1  20
#define N2 /*256*/ 50

short       stab[N2+3];
int         itab[N2+3];
float       ftab[N2+3];
double      dtab[N2+3];
long double ldtab[N2+3];

int main ( void )
{
  int i;
  pkv_InitScratchMem ( 655360 );

  printf ( "short:\n" );
  for ( i = 0; i < N1; i += 4 ) {
    stab[i] = 2*i+100;  stab[i+1] = -2*i-1+100;
    stab[i+2] = 2*i+1+100;  stab[i+3] = -2*i+100;
  }
  if ( pkv_SortFast ( sizeof(short), ID_SIGNED_INT, sizeof(short), 0,
                      N1, stab ) != SORT_OK ) {
    printf ( "qq\n" );
    exit ( 1 );
  }
  for ( i = 0; i < N1; i++ )
    printf ( "%d ", stab[i] );
  printf ( "\n" );
  for ( i = 0; i < N2; i += 4 ) {
    stab[i] = 2*i+100;  stab[i+1] = -2*i-1+100;
    stab[i+2] = 2*i+1+100;  stab[i+3] = -2*i+100;
  }
  if ( pkv_SortFast ( sizeof(short), ID_SIGNED_INT, sizeof(short), 0,
                      N2, stab ) != SORT_OK ) {
    printf ( "qq\n" );
    exit ( 1 );
  }
  for ( i = 0; i < N2; i++ )
    printf ( "%d ", stab[i] );
  printf ( "\n" );

  printf ( "int:\n" );
  for ( i = 0; i < N1; i += 4 ) {
    itab[i] = 2*i;  itab[i+1] = -2*i-1;
    itab[i+2] = 2*i+1;  itab[i+3] = -2*i;
  }
  if ( pkv_SortFast ( sizeof(int), ID_SIGNED_INT, sizeof(int), 0,
                      N1, itab ) != SORT_OK ) {
    printf ( "qq\n" );
    exit ( 1 );
  }
  for ( i = 0; i < N1; i++ )
    printf ( "%d ", itab[i] );
  printf ( "\n" );
  for ( i = 0; i < N2; i += 4 ) {
    itab[i] = 2*i;  itab[i+1] = -2*i-1;
    itab[i+2] = 2*i+1;  itab[i+3] = -2*i;
  }
  if ( pkv_SortFast ( sizeof(int), ID_SIGNED_INT, sizeof(int), 0,
                      N2, itab ) != SORT_OK ) {
    printf ( "qq\n" );
    exit ( 1 );
  }
  for ( i = 0; i < N2; i++ )
    printf ( "%d ", itab[i] );
  printf ( "\n" );

  printf ( "float:\n" );
  for ( i = 0; i < N1; i += 4 ) {
    ftab[i] = 2*i;  ftab[i+1] = -2*i-1;
    ftab[i+2] = 2*i+1;  ftab[i+3] = -2*i;
  }
  if ( pkv_SortFast ( sizeof(float), ID_IEEE754_FLOAT, sizeof(float), 0,
                      N1, ftab ) != SORT_OK ) {
    printf ( "qq\n" );
    exit ( 1 );
  }
  for ( i = 0; i < N1; i++ )
    printf ( "%6.1f ", ftab[i] );
  printf ( "\n" );
  for ( i = 0; i < N2; i += 4 ) {
    ftab[i] = 2*i;  ftab[i+1] = -2*i-1;
    ftab[i+2] = 2*i+1;  ftab[i+3] = -2*i;
  }
  if ( pkv_SortFast ( sizeof(float), ID_IEEE754_FLOAT, sizeof(float), 0,
                      N2, ftab ) != SORT_OK ) {
    printf ( "qq\n" );
    exit ( 1 );
  }
  for ( i = 0; i < N2; i++ )
    printf ( "%6.1f ", ftab[i] );
  printf ( "\n" );

  printf ( "double:\n" );
  for ( i = 0; i < N1; i += 4 ) {
    dtab[i] = 2*i;  dtab[i+1] = -2*i-1;
    dtab[i+2] = 2*i+1;  dtab[i+3] = -2*i;
  }
  if ( pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(double), 0,
                      N1, dtab ) != SORT_OK ) {
    printf ( "qq\n" );
    exit ( 1 );
  }
  for ( i = 0; i < N1; i++ )
    printf ( "%6.1f ", dtab[i] );
  printf ( "\n" );
  for ( i = 0; i < N2; i += 4 ) {
    dtab[i] = 2*i;  dtab[i+1] = -2*i-1;
    dtab[i+2] = 2*i+1;  dtab[i+3] = -2*i;
  }
  if ( pkv_SortFast ( sizeof(double), ID_IEEE754_DOUBLE, sizeof(double), 0,
                      N2, dtab ) != SORT_OK ) {
    printf ( "qq\n" );
    exit ( 1 );
  }
  for ( i = 0; i < N2; i++ )
    printf ( "%6.1f ", dtab[i] );
  printf ( "\n" );

  printf ( "long double:\n" );
  for ( i = 0; i < N1; i += 4 ) {
    ldtab[i] = 2*i;  ldtab[i+1] = -2*i-1;
    ldtab[i+2] = 2*i+1;  ldtab[i+3] = -2*i;
  }
  if ( pkv_SortFast ( sizeof(long double), ID_IEEE754_LONG_DOUBLE,
                      sizeof(long double), 0, N1, ldtab ) != SORT_OK ) {
    printf ( "qq\n" );
    exit ( 1 );
  }
  for ( i = 0; i < N1; i++ )
    printf ( "%6.1Lf ", ldtab[i] );
  printf ( "\n" );
  for ( i = 0; i < N2; i += 4 ) {
    ldtab[i] = 2*i;  ldtab[i+1] = -2*i-1;
    ldtab[i+2] = 2*i+1;  ldtab[i+3] = -2*i;
  }
  if ( pkv_SortFast ( sizeof(long double), ID_IEEE754_LONG_DOUBLE,
                      sizeof(long double), 0, N2, ldtab ) != SORT_OK ) {
    printf ( "qq\n" );
    exit ( 1 );
  }
  for ( i = 0; i < N2; i++ )
    printf ( "%6.1Lf ", ldtab[i] );
  printf ( "\n" );

  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

