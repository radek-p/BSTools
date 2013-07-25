
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

static void Solve2x2f ( const vector2f *a1, const vector2f *a2,
                        const vector2f *b, float *x1, float *x2 )
{
  float aa[4], bb[2];

  aa[0] = a1->x;  aa[1] = a2->x;  bb[0] = b->x;
  aa[2] = a1->y;  aa[3] = a2->y;  bb[1] = b->y;
  pkn_multiGaussSolveLinEqf ( 2, aa, 1, 1, bb );
  *x1 = bb[0];
  *x2 = bb[1];
} /*Solve2x2f*/

/* compute function values from the first order derivatives */
/* compatibility equation */
#define SolveCompatibilityEq1f(ru,rs,qv,b1,c1) \
  Solve2x2f ( ru, rs, qv, b1, c1 )

static void SolveCompatibilityEq2af (
            const vector2f *ru, const vector2f *rs,
            const vector2f *ruu, const vector2f *rus,
            const vector2f *qv, const vector2f *qt,
            const vector2f *qvv, const vector2f *qvt,
            float b1, float c1, float f1, float g1,
            float *b1d, float c1d, float *f1d, float g1d )
{
/* compute derivatives of b1, f1 from mixed partial derivatives */
/* compatibility equation */
  vector2f b;

  SetVector2f ( &b,
                c1d*rs->x+b1*ruu->x+c1*rus->x-g1d*qt->x-f1*qvv->x-g1*qvt->x,
                c1d*rs->y+b1*ruu->y+c1*rus->y-g1d*qt->y-f1*qvv->y-g1*qvt->y );
  Solve2x2f ( ru, qv, &b, b1d, f1d );
  *b1d = -(*b1d);
} /*SolveCompatibilityEq2af*/

static void SolveCompatibilityEq2cf (
            const vector2f *ru, const vector2f *rs, const vector2f *rus,
            const vector2f *qv, const vector2f *qt, const vector2f *qvt,
            float c1, float g1,
            float *b1d, float c1d, float *f1d, float g1d )
{
/* compute derivatives of b1, f1 from mixed partial derivatives */
/* compatibility equation; in this version b1 and f1 are 0 */
  vector2f b;

  SetVector2f ( &b, c1d*rs->x+c1*rus->x-g1d*qt->x-g1*qvt->x,
                    c1d*rs->y+c1*rus->y-g1d*qt->y-g1*qvt->y );
  Solve2x2f ( ru, qv, &b, b1d, f1d );
  *b1d = -(*b1d);
} /*SolveCompatibilityEq2cf*/

static void SolveCompatibilityEq2df (
            const vector2f *ru, const vector2f *rs, const vector2f *rus,
            const vector2f *qv, const vector2f *qt, const vector2f *qvt,
            float c1, float g1, float b1d, float c1d, float *f1d, float *g1d )
{
/* compute derivatives of f1, g1 from the mixed partial derivatives */
/* compatibility equation; in this version b1 and f1 are 0 */
  vector2f b;

  SetVector2f ( &b, ru->x*b1d+rs->x*c1d+rus->x*c1-qvt->x*g1,
                    ru->y*b1d+rs->y*c1d+rus->y*c1-qvt->y*g1 );
  Solve2x2f ( qv, qt, &b, f1d, g1d );
} /*SolveCompatibilityEq2df*/

