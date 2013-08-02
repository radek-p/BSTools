
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
boolean myfunc ( int n, void *usrdata, double *x, double *f )
{
  double xx, yy, zz, a, b;

  xx = X*X;
  yy = Y*Y;
  zz = Z*Z;
  a = (xx+yy+zz-1.0);
  b = (xx+yy+zz-4.0);
  *f = a*b + 0.1*(X+yy+zz) + 10.0;
  return true;
} /*myfunc*/

boolean myfuncg ( int n, void *usrdata, double *x, double *f, double *g )
{
  double xx, yy, zz, a, b;

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

boolean myfuncgh ( int n, void *usrdata,
                   double *x, double *f, double *g, double *h )
{
  double xx, yy, zz, a, b;

  xx = X*X;
  yy = Y*Y;
  zz = Z*Z;
  a = (xx+yy+zz-1.0);
  b = (xx+yy+zz-4.0);
  *f = a*b + 0.1*(X+yy+zz) + 10.0;

  g[0] = 2.0*X*(a+b)+0.1;
  g[1] = 2.0*Y*(a+b)+0.2*Y;
  g[2] = 2.0*Z*(a+b)+0.2*Z;

  h[0] = 2.0*(a+b)+8.0*xx;
  h[1] = 8.0*X*Y;
  h[2] = 2.0*(a+b)+8.0*yy+0.2;
  h[3] = 8.0*X*Z;
  h[4] = 8.0*Y*Z;
  h[5] = 2.0*(a+b)+8.0*zz+0.2;
  return true;
} /*myfuncgh*/
#undef X
#undef Y
#undef Z

void TestGrad ( int n, void *usrdata,
                pkn_NLMTevalfuncd funcf, pkn_NLMTevalfuncgd funcfg,
                double *x, double *f, double *grf, double *grt )
{
#define EPS 1.0e-8
  double s, _f;
  int    i;

  funcf ( n, usrdata, x, f );
  funcfg ( n, usrdata, x, &_f, grf );
  for ( i = 0; i < n; i++ ) {
    s = x[i];
    x[i] += EPS;
    funcf ( n, usrdata, x, &_f );
    grt[i] = (_f-*f)/EPS;
    x[i] = s;
  }
#undef EPS
} /*TestGrad*/

void TestHess ( int n, void *usrdata,
                pkn_NLMTevalfuncgd funcfg, pkn_NLMTevalfuncghd funcfgh,
                double *x, double *f, double *grf, double *hf, double *ht )
{
#define EPS 1.0e-8
  void   *sp;
  double s, _f, *grt;
  int    i, j;

  sp = pkv_GetScratchMemTop ();
  grt = pkv_GetScratchMemd ( n );
  if ( !grt )
    return;
  funcfgh ( n, usrdata, x, f, grf, hf );
  for ( i = 0; i < n; i++ ) {
    s = x[i];
    x[i] += EPS;
    funcfg ( n, usrdata, x, &_f, grt );
    for( j = 0; j <= i; j++ )
      ht[pkn_LowerTrMatIndex(i,j)] = (grt[j]-grf[j])/EPS;
    x[i] = s;
  }
  pkv_SetScratchMemTop ( sp );
#undef EPS
}/*TestHess*/

int main ( void )
{
  double  x[N] = {0.0,1.0,1.5};
  double  nu, f;
  double  grf[N], grt[N], hf[(N*(N+1))/2], ht[(N*(N+1))/2];
  int     i;

  pkv_InitScratchMem ( 655360 );

  printf ( "test gradientu:\n" );
  TestGrad ( N, NULL, myfunc, myfuncg, x, &f, grf, grt );
  for ( i = 0; i < N; i++ )
    printf ( "%14.8f, %14.8f\n", grf[i], grt[i] );
  printf ( "test hesjanu:\n" );
  TestHess ( N, NULL, myfuncg, myfuncgh, x, &f, grf, hf, ht );
  for ( i = 0; i < (N*(N+1))/2; i++ )
    printf ( "%14.8f, %14.8f\n", hf[i], ht[i] );
  printf ( "\n" );

  nu = -1.0;
  for ( i = 0; i < 20; i++ ) {
    switch ( pkn_NLMIterd ( N, NULL, x, myfunc, myfuncg, myfuncgh, NULL, NULL,
                         0.0, 1.0e-6, 1.0e-6, &nu ) ) {
  case PKN_LMT_ERROR:
      printf ( "error\n" );
      goto way_out;
  case PKN_LMT_CONTINUE_LM:
      myfunc ( N, NULL, x, &f );
      printf ( "- %3d: x=(%f,%f,%f), f = %f, nu = %f\n", i, x[0], x[1], x[2], f, nu );
      break;
  case PKN_LMT_CONTINUE_N:
      myfunc ( N, NULL, x, &f );
      printf ( "+ %3d: x=(%f,%f,%f), f = %f, nu = %f\n", i, x[0], x[1], x[2], f, nu );
      break;
  case PKN_LMT_FOUND_MINIMUM:
      printf ( "found minimum\n" );
      goto way_out;
  case PKN_LMT_FOUND_ZEROGRAD:
      printf ( "found gradient zero\n" );
      goto way_out;
  case PKN_LMT_FOUND_BARRIER:
      printf ( "found barrier\n" );
      goto way_out;
    }
  }
way_out:
  pkv_DestroyScratchMem ();
  exit ( 0 );
} /*main*/

