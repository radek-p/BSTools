
/* this program may help to prepare data, i.e. initial closed curves */
/* being knots in R^3, for minimization of integral Menger curvature */
/* -just uncomment one of sets of parameters below, recompile it and */
/* execute */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "bsmesh.h"
#include "bsfile.h"

#define DEG  3

/*
#define LKN 14
#define N1   1
#define N2   1
#define R1   1.24264068711929
#define R2   0.0
char fn[] = "c11_14.bs";
*/
/*
#define LKN 39
#define N1   2
#define N2   3
#define R1   1.0
#define R2   0.5
char fn[] = "c23_39.bs";
*/
/*
#define LKN 54
#define N1   3
#define N2   5
#define R1   1.0
#define R2   0.5
char fn[] = "c35_54.bs";
*/
/*
#define LKN 62
#define N1   4
#define N2   5
#define R1   1.0
#define R2   0.5
char fn[] = "c45_62.bs";
*/
/*
#define LKN 20
#define N1   1
#define N2   0
#define R1   1.0
#define R2   0.0
char fn[] = "c10_20.bs";
*/
/*
#define LKN 54
#define N1   2
#define N2   5
#define R1   1.0
#define R2   0.5
char fn[] = "c25_54.bs";
*/
/*
#define LKN 54
#define N1   2
#define N2   7
#define R1   1.0
#define R2   0.5
char fn[] = "c27_54.bs";
*/
/*
#define LKN 78
#define N1   2
#define N2   9
#define R1   1.0
#define R2   0.2
char fn[] = "c29_78.bs";
*/

#define LKN 94
#define N1   2
#define N2   11
#define R1   1.0
#define R2   0.2
char fn[] = "c211_94.bs";


double  knots[LKN+1];
point3d cpoints[LKN-DEG];

int main ( void )
{
  int     i, ncp;
  double  alpha, beta;
  trans3d tr1, tr2;
  point3d p, q;

  for ( i = 0; i <= LKN; i++ )
    knots[i] = (double)i;
  ncp = LKN-2*DEG;
  alpha = 2.0*PI*(double)N1/(double)ncp;
  beta = 2.0*PI*(double)N2/(double)ncp;
  IdentTrans3d ( &tr1 );
  IdentTrans3d ( &tr2 );
  RotYTrans3d ( &tr2, beta );
  SetPoint3d ( &p, R2, 0.0, 0.0 );
  for ( i = 0; i < ncp; i++ ) {
    TransPoint3d ( &tr2, &p, &p );
    SetPoint3d ( &q, p.x+R1, p.y, p.z );
    TransPoint3d ( &tr1, &q, &cpoints[i] );
    RotZTrans3d ( &tr1, alpha );
  }
  memcpy ( &cpoints[ncp], cpoints, DEG*sizeof(point3d) );
  if ( bsf_OpenOutputFile ( fn, false ) ) {
    bsf_WriteBSplineCurved ( 3, 3, false, DEG, LKN, knots, true,
                             (double*)cpoints, NULL, -1, NULL, NULL );
    bsf_CloseOutputFile ();
    printf ( "%s\n", fn );
  }
  else {
    printf ( "klops\n" );
    exit ( 1 );
  }

  exit ( 0 );
} /*main*/

