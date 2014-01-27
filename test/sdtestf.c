
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"

#define N 3

#define X x[0]
#define Y x[1]
#define Z x[2]
boolean myfunc ( int n, void *usrdata, float *x, float *f )
{
  float xx, yy, zz, a, b;

  xx = X*X;
  yy = Y*Y;
  zz = Z*Z;
  a = (xx+yy+zz-1.0);
  b = (xx+yy+zz-4.0);
  *f = a*b + 0.1*(X+yy+zz) + 10.0;
  return true;
} /*myfunc*/

boolean myfuncg ( int n, void *usrdata, float *x, float *f, float *g )
{
  float xx, yy, zz, a, b;

  xx = X*X;
  yy = Y*Y;
  zz = Z*Z;
  a = (xx+yy+zz-1.0);
  b = (xx+yy+zz-4.0);
  *f = a*b + 0.1*(X+yy+zz) + 10.0;

  g[0] = 2.0*X*(a+b)+0.1;
  g[1] = 2.0*Y*(a+b)+0.2*Y;
  g[2] = 2.0*Z*(a+b)+0.2*Z;
  return true;
} /*myfuncg*/
#undef X
#undef Y
#undef Z

void TestGrad ( int n, void *usrdata,
                pkn_NLMTevalfuncf funcf, pkn_NLMTevalfuncgf funcfg,
                float *x, float *f, float *grf, float *grt )
{
#define EPS 1.0e-4
  float s, _f1, _f2;
  int   i;

  funcf ( n, usrdata, x, f );
  funcfg ( n, usrdata, x, &_f1, grf );
  for ( i = 0; i < n; i++ ) {
    s = x[i];
    x[i] += EPS;
    funcf ( n, usrdata, x, &_f1 );
    x[i] = s-EPS;
    funcf ( n, usrdata, x, &_f2 );
    grt[i] = (_f1-_f2)/(EPS+EPS);
    x[i] = s;
  }
#undef EPS
} /*TestGrad*/

int main ( void )
{
  float  x[N] = {0.0,1.0,1.5};
  float  nu, f;
  float  grf[N], grt[N];
  int     i;

  pkv_InitScratchMem ( 655360 );

  printf ( "test gradientu:\n" );
  TestGrad ( N, NULL, myfunc, myfuncg, x, &f, grf, grt );
  for ( i = 0; i < N; i++ )
    printf ( "%14.8f, %14.8f\n", grf[i], grt[i] );
  printf ( "\n" );

  nu = -1.0;
  for ( i = 0; i < 300; i++ ) {
    switch ( pkn_SDIterf ( N, NULL, x, myfunc, myfuncg, NULL, NULL,
                           0.0, 1.0e-6, 1.0e-6, &nu ) ) {
  case PKN_SD_ERROR:
      printf ( "error\n" );
      goto way_out;
  case PKN_SD_CONTINUE:
      if ( !(i % 10) ) {
        myfunc ( N, NULL, x, &f );
        printf ( "  %3d: x=(%f,%f,%f), f = %f, nu = %f\n",
                 i, x[0], x[1], x[2], f, nu );
      }
      break;
  case PKN_SD_FOUND_ZEROGRAD:
      printf ( "found gradient zero\n" );
      goto way_out;
  case PKN_SD_FOUND_BARRIER:
      printf ( "found barrier\n" );
      goto way_out;
  case PKN_SD_CROSSED_LIMIT:
      printf ( "crossed limit\n" );
      goto way_out;
  case PKN_SD_NO_PROGRESS:
      printf ( "no progress\n" );
      goto way_out;
  default:
      printf ( "something wrong\n" );
      goto way_out;
    }
  }
way_out:
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

