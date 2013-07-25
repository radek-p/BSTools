
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

static void Solve2x2d ( const vector2d *a1, const vector2d *a2,
                        const vector2d *b, double *x1, double *x2 )
{
  double aa[4], bb[2];

  aa[0] = a1->x;  aa[1] = a2->x;  bb[0] = b->x;
  aa[2] = a1->y;  aa[3] = a2->y;  bb[1] = b->y;
  pkn_multiGaussSolveLinEqd ( 2, aa, 1, 1, bb );
  *x1 = bb[0];
  *x2 = bb[1];
} /*Solve2x2d*/

/* compute function values from the first order derivatives */
/* compatibility equation */
#define SolveCompatibilityEq1d(ru,rs,qv,b1,c1) \
  Solve2x2d ( ru, rs, qv, b1, c1 )

static void SolveCompatibilityEq2ad (
            const vector2d *ru, const vector2d *rs,
            const vector2d *ruu, const vector2d *rus,
            const vector2d *qv, const vector2d *qt,
            const vector2d *qvv, const vector2d *qvt,
            double b1, double c1, double f1, double g1,
            double *b1d, double c1d, double *f1d, double g1d )
{
/* compute derivatives of b1, f1 from mixed partial derivatives */
/* compatibility equation */
  vector2d b;

  SetVector2d ( &b,
                c1d*rs->x+b1*ruu->x+c1*rus->x-g1d*qt->x-f1*qvv->x-g1*qvt->x,
                c1d*rs->y+b1*ruu->y+c1*rus->y-g1d*qt->y-f1*qvv->y-g1*qvt->y );
  Solve2x2d ( ru, qv, &b, b1d, f1d );
  *b1d = -(*b1d);
} /*SolveCompatibilityEq2ad*/

static void SolveCompatibilityEq2cd (
            const vector2d *ru, const vector2d *rs, const vector2d *rus,
            const vector2d *qv, const vector2d *qt, const vector2d *qvt,
            double c1, double g1,
            double *b1d, double c1d, double *f1d, double g1d )
{
/* compute derivatives of b1, f1 from mixed partial derivatives */
/* compatibility equation; in this version b1 and f1 are 0 */
  vector2d b;

  SetVector2d ( &b, c1d*rs->x+c1*rus->x-g1d*qt->x-g1*qvt->x,
                    c1d*rs->y+c1*rus->y-g1d*qt->y-g1*qvt->y );
  Solve2x2d ( ru, qv, &b, b1d, f1d );
  *b1d = -(*b1d);
} /*SolveCompatibilityEq2cd*/

static void SolveCompatibilityEq2dd (
            const vector2d *ru, const vector2d *rs, const vector2d *rus,
            const vector2d *qv, const vector2d *qt, const vector2d *qvt,
            double c1, double g1, double b1d, double c1d, double *f1d, double *g1d )
{
/* compute derivatives of f1, g1 from the mixed partial derivatives */
/* compatibility equation; in this version b1 and f1 are 0 */
  vector2d b;

  SetVector2d ( &b, ru->x*b1d+rs->x*c1d+rus->x*c1-qvt->x*g1,
                    ru->y*b1d+rs->y*c1d+rus->y*c1-qvt->y*g1 );
  Solve2x2d ( qv, qt, &b, f1d, g1d );
} /*SolveCompatibilityEq2dd*/

