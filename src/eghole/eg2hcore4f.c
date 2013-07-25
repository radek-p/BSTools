
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005,2007                             */
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

static void SolveCompatibilityEq2bf (
            const vector2f *ru, const vector2f *rs,
            const vector2f *ruu, const vector2f *rus,
            const vector2f *qv, const vector2f *qt,
            const vector2f *qvv, const vector2f *qvt,
            float b1, float c1, float f1, float g1,
            float b1d, float c1d, float *f1d, float *g1d )
{
/* compute derivatives of f1, g1 from the mixed partial derivatives */
/* compatibility equation */
  vector2f b;

  SetVector2f ( &b,
                ru->x*b1d+rs->x*c1d+ruu->x*b1+rus->x*c1-qvv->x*f1-qvt->x*g1,
                ru->y*b1d+rs->y*c1d+ruu->y*b1+rus->y*c1-qvv->y*f1-qvt->y*g1 );
  Solve2x2f ( qv, qt, &b, f1d, g1d );
} /*SolveCompatibilityEq2bf*/

static void SolveCompatibilityEq3f (
            const vector2f *ru, const vector2f *rs, const vector2f *ruu,
            const vector2f *rus, const vector2f *rss, const vector2f *qvv,
            float b1, float c1, float *b2, float *c2 )
{
/* compute values of b2, c2 from the second order partial derivatives */
/* compatibility equation */
  vector2f b;

  SetVector2f ( &b,
                (float)(qvv->x-b1*b1*ruu->x-2.0*b1*c1*rus->x-c1*c1*rss->x),
                (float)(qvv->y-b1*b1*ruu->y-2.0*b1*c1*rus->y-c1*c1*rss->y) );
  Solve2x2f ( ru, rs, &b, b2, c2 );
} /*SolveCompatibilityEq3f*/

static void SolveCompatibilityEq4af (
            const vector2f *qv, const vector2f *qt, const vector2f *qvv,
            const vector2f *qvt, const vector2f *qvvv, const vector2f *qvvt,
            const vector2f *ru, const vector2f *rs, const vector2f *ruu,
            const vector2f *rus, const vector2f *rss, const vector2f *ruuu,
            const vector2f *ruus, const vector2f *russ,
            float b1, float b1d, float c1, float c1d,
            float f1, float f1d, float f1dd, float g1, float g1d, float g1dd,
            float b2, float *b2d, float c2, float *c2d )
{
/* compute derivatives of b2, c2, from the third order mixed partial */
/* derivatives compatibility equation */
  vector2f b;

  SetVector2f ( &b,
    (float)(f1dd*qv->x+g1dd*qt->x+2.0*(f1d*qvv->x+g1d*qvt->x)+f1*qvvv->x+g1*qvvt->x-
    ((2.0*b1*b1d+b2)*ruu->x+(2.0*(b1d*c1+b1*c1d)+c2)*rus->x+
    2.0*c1*c1d*rss->x+b1*b1*ruuu->x+2.0*b1*c1*ruus->x+c1*c1*russ->x)),
    (float)(f1dd*qv->y+g1dd*qt->y+2.0*(f1d*qvv->y+g1d*qvt->y)+f1*qvvv->y+g1*qvvt->y-
    ((2.0*b1*b1d+b2)*ruu->y+(2.0*(b1d*c1+b1*c1d)+c2)*rus->y+
    2.0*c1*c1d*rss->y+b1*b1*ruuu->y+2.0*b1*c1*ruus->y+c1*c1*russ->y)) );
  Solve2x2f ( ru, rs, &b, b2d, c2d );
} /*SolveCompatibilityEq4af*/

static void SolveCompatibilityEq4bf (
            const vector2f *qv, const vector2f *qt, const vector2f *qvv,
            const vector2f *qvt, const vector2f *qvvt,
            const vector2f *ru, const vector2f *rs,
            const vector2f *rus, const vector2f *rss, const vector2f *russ,
            float b1d, float c1, float c1d, float f1d, float f1dd,
            float g1, float g1d, float g1dd,
            float *b2d, float c2, float *c2d )
{
/* compute derivatives of b2, c2, from the third order mixed partial */
/* derivatives compatibility equation; this version assumes that */
/* b1 = f1 = b2 = 0 */
  vector2f b;

  SetVector2f ( &b,
    (float)(f1dd*qv->x+g1dd*qt->x+2.0*(f1d*qvv->x+g1d*qvt->x)+g1*qvvt->x-
    ((2.0*b1d*c1+c2)*rus->x+2.0*c1*c1d*rss->x+c1*c1*russ->x)),
    (float)(f1dd*qv->y+g1dd*qt->y+2.0*(f1d*qvv->y+g1d*qvt->y)+g1*qvvt->y-
    ((2.0*b1d*c1+c2)*rus->y+2.0*c1*c1d*rss->y+c1*c1*russ->y)) );
  Solve2x2f ( ru, rs, &b, b2d, c2d );
} /*SolveCompatibilityEq4bf*/

static void SolveCompatibilityEq4cf (
            const vector2f *qv, const vector2f *qt, const vector2f *qvv,
            const vector2f *qvt, const vector2f *qvvt,
            const vector2f *ru, const vector2f *rs,
            const vector2f *rus, const vector2f *rss, const vector2f *russ,
            float b1d, float c1, float c1d, float f1d, float *f1dd,
            float g1, float g1d, float *g1dd,
            float b2d, float c2, float c2d )
{
/* compute second order derivatives of f1, g1, from the third order */
/* mixed partial derivatives compatibility equation; this version assumes */
/* that b1 = f1 = b2 = 0 */
  vector2f b;

  SetVector2f ( &b,
    (float)(-2.0*f1d*qvv->x -2.0*g1d*qvt->x -g1*qvvt->x +b2d*ru->x +c2d*rs->x
    +(2.0*b1d*c1+c2)*rus->x +2.0*c1*c1d*rss->x +c1*c1*russ->x),
    (float)(-2.0*f1d*qvv->y -2.0*g1d*qvt->y -g1*qvvt->y +b2d*ru->y +c2d*rs->y
    +(2.0*b1d*c1+c2)*rus->y +2.0*c1*c1d*rss->y +c1*c1*russ->y) );
  Solve2x2f ( qv, qt, &b, f1dd, g1dd );
} /*SolveCompatibilityEq4cf*/

static void SolveCompatibilityEq5af (
            const vector2f *qv, const vector2f *qt, const vector2f *qvv,
            const vector2f *qvt, const vector2f *qtt, const vector2f *qvvv,
            const vector2f *qvvt, const vector2f *qvtt, const vector2f *qvvvv,
            const vector2f *qvvvt, const vector2f *qvvtt,
            const vector2f *ru, const vector2f *rs, const vector2f *ruu,
            const vector2f *rus, const vector2f *rss, const vector2f *ruuu,
            const vector2f *ruus, const vector2f *russ, const vector2f *ruuuu,
            const vector2f *ruuus, const vector2f *ruuss,
            float b1, float b1d, float b1dd, float c1, float c1d, float c1dd,
            float f1, float f1d, float f1dd, float g1, float g1d, float g1dd,
            float b2, float b2d, float *b2dd, float c2, float c2d, float c2dd,
            float f2, float f2d, float *f2dd, float g2, float g2d, float g2dd )
{
/* compute second order derivatives of b2, f2, from the fourth order mixed */
/* partial derivatives compatibility equation */
  vector2f b;

  SetVector2f ( &b,
     (float)(g2dd*qt->x+2.0*(f1*f1dd+f1d*f1d+f2d)*qvv->x+
     2.0*(f1dd*g1+2.0*f1d*g1d+f1*g1dd+g2d)*qvt->x+
     2.0*(g1*g1dd+g1d*g1d)*qtt->x+(4.0*f1*f1d+f2)*qvvv->x+
     (4.0*(f1d*g1+f1*g1d)+g2)*qvvt->x+4.0*g1*g1d*qvtt->x+
     f1*f1*qvvvv->x+2.0*f1*g1*qvvvt->x+g1*g1*qvvtt->x -
     (c2dd*rs->x+2.0*(b1*b1dd+b1d*b1d+b2d)*ruu->x+
     2.0*(b1dd*c1+2.0*b1d*c1d+b1*c1dd+c2d)*rus->x+
     2.0*(c1*c1dd+c1d*c1d)*rss->x+(4.0*b1*b1d+b2)*ruuu->x+
     (4.0*(b1d*c1+b1*c1d)+c2)*ruus->x+4.0*c1*c1d*russ->x+
     b1*b1*ruuuu->x+2.0*b1*c1*ruuus->x+c1*c1*ruuss->x)),
     (float)(g2dd*qt->y+2.0*(f1*f1dd+f1d*f1d+f2d)*qvv->y+
     2.0*(f1dd*g1+2.0*f1d*g1d+f1*g1dd+g2d)*qvt->y+
     2.0*(g1*g1dd+g1d*g1d)*qtt->y+(4.0*f1*f1d+f2)*qvvv->y+
     (4.0*(f1d*g1+f1*g1d)+g2)*qvvt->y+4.0*g1*g1d*qvtt->y+
     f1*f1*qvvvv->y+2.0*f1*g1*qvvvt->y+g1*g1*qvvtt->y -
     (c2dd*rs->y+2.0*(b1*b1dd+b1d*b1d+b2d)*ruu->y+
     2.0*(b1dd*c1+2.0*b1d*c1d+b1*c1dd+c2d)*rus->y+
     2.0*(c1*c1dd+c1d*c1d)*rss->y+(4.0*b1*b1d+b2)*ruuu->y+
     (4.0*(b1d*c1+b1*c1d)+c2)*ruus->y+4.0*c1*c1d*russ->y+
     b1*b1*ruuuu->y+2.0*b1*c1*ruuus->y+c1*c1*ruuss->y)) );
  Solve2x2f ( ru, qv, &b, b2dd, f2dd );
  *f2dd = -(*f2dd);
} /*SolveCompatibilityEq5af*/

static void SolveCompatibilityEq5bf (
            const vector2f *qv, const vector2f *qt, const vector2f *qvv,
            const vector2f *qvt, const vector2f *qtt,
            const vector2f *qvvt, const vector2f *qvtt, const vector2f *qvvtt,
            const vector2f *ru, const vector2f *rs, const vector2f *ruu,
            const vector2f *rus, const vector2f *rss,
            const vector2f *ruus, const vector2f *russ, const vector2f *ruuss,
            float b1d, float b1dd, float c1, float c1d, float c1dd,
            float f1d, float f1dd, float g1, float g1d, float g1dd,
            float b2d, float *b2dd, float c2, float c2d, float *c2dd,
            float f2d, float f2dd, float g2, float g2d, float g2dd )
{
/* compute second order derivatives of b2, c2, from the fourth order mixed */
/* partial derivatives compatibility equation; this version assumes that */
/* b1 = f1 = b2 = f2 = 0 */
  vector2f b;

  SetVector2f ( &b,
     (float)(f2dd*qv->x+g2dd*qt->x+2.0*(f1d*f1d+f2d)*qvv->x+
     2.0*(f1dd*g1+2.0*f1d*g1d+g2d)*qvt->x+2.0*(g1*g1dd+g1d*g1d)*qtt->x+
     (4.0*f1d*g1+g2)*qvvt->x+4.0*g1*g1d*qvtt->x+g1*g1*qvvtt->x -
     (2.0*(b1d*b1d+b2d)*ruu->x+2.0*(b1dd*c1+2.0*b1d*c1d+c2d)*rus->x+
     2.0*(c1*c1dd+c1d*c1d)*rss->x+(4.0*b1d*c1+c2)*ruus->x+
     4.0*c1*c1d*russ->x+c1*c1*ruuss->x)),
     (float)(f2dd*qv->y+g2dd*qt->y+2.0*(f1d*f1d+f2d)*qvv->y+
     2.0*(f1dd*g1+2.0*f1d*g1d+g2d)*qvt->y+2.0*(g1*g1dd+g1d*g1d)*qtt->y+
     (4.0*f1d*g1+g2)*qvvt->y+4.0*g1*g1d*qvtt->y+g1*g1*qvvtt->y -
     (2.0*(b1d*b1d+b2d)*ruu->y+2.0*(b1dd*c1+2.0*b1d*c1d+c2d)*rus->y+
     2.0*(c1*c1dd+c1d*c1d)*rss->y+(4.0*b1d*c1+c2)*ruus->y+
     4.0*c1*c1d*russ->y+c1*c1*ruuss->y)) );
  Solve2x2f ( ru, rs, &b, b2dd, c2dd );
} /*SolveCompatibilityEq5bf*/

static void SolveCompatibilityEq5cf (
            const vector2f *qv, const vector2f *qt, const vector2f *qvv,
            const vector2f *qvt, const vector2f *qtt, const vector2f *qvvt,
            const vector2f *qvtt, const vector2f *qvvtt,
            const vector2f *ru, const vector2f *rs, const vector2f *ruu,
            const vector2f *rus, const vector2f *rss, const vector2f *ruus,
            const vector2f *russ, const vector2f *ruuss,
            float b1d, float b1dd, float c1, float c1d, float c1dd,
            float f1d, float f1dd, float g1, float g1d, float g1dd,
            float b2d, float *b2dd, float c2, float c2d, float c2dd,
            float f2d, float *f2dd, float g2, float g2d, float g2dd )
{
/* compute second order derivatives of b2, f2, from the fourth order mixed */
/* partial derivatives compatibility equation; this version assumes that */
/* b1 = f1 = b2 = f2 = 0 */
  vector2f b;

  SetVector2f ( &b,
     (float)(g2dd*qt->x+2.0*(f1d*f1d+f2d)*qvv->x+2.0*(f1dd*g1+2.0*f1d*g1d+g2d)*qvt->x+
     2.0*(g1*g1dd+g1d*g1d)*qtt->x+(4.0*f1d*g1+g2)*qvvt->x+4.0*g1*g1d*qvtt->x+
     g1*g1*qvvtt->x -
     (c2dd*rs->x+2.0*(b1d*b1d+b2d)*ruu->x+2.0*(b1dd*c1+2.0*b1d*c1d+c2d)*rus->x+
     2.0*(c1*c1dd+c1d*c1d)*rss->x+(4.0*b1d*c1+c2)*ruus->x+4.0*c1*c1d*russ->x+
     c1*c1*ruuss->x)),
     (float)(g2dd*qt->y+2.0*(f1d*f1d+f2d)*qvv->y+2.0*(f1dd*g1+2.0*f1d*g1d+g2d)*qvt->y+
     2.0*(g1*g1dd+g1d*g1d)*qtt->y+(4.0*f1d*g1+g2)*qvvt->y+4.0*g1*g1d*qvtt->y+
     g1*g1*qvvtt->y -
     (c2dd*rs->y+2.0*(b1d*b1d+b2d)*ruu->y+2.0*(b1dd*c1+2.0*b1d*c1d+c2d)*rus->y+
     2.0*(c1*c1dd+c1d*c1d)*rss->y+(4.0*b1d*c1+c2)*ruus->y+4.0*c1*c1d*russ->y+
     c1*c1*ruuss->y)) );
  Solve2x2f ( ru, qv, &b, b2dd, f2dd );
  *f2dd = -(*f2dd);
} /*SolveCompatibilityEq5cf*/

