
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005,2007                             */
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

static void SolveCompatibilityEq2bd (
            const vector2d *ru, const vector2d *rs,
            const vector2d *ruu, const vector2d *rus,
            const vector2d *qv, const vector2d *qt,
            const vector2d *qvv, const vector2d *qvt,
            double b1, double c1, double f1, double g1,
            double b1d, double c1d, double *f1d, double *g1d )
{
/* compute derivatives of f1, g1 from the mixed partial derivatives */
/* compatibility equation */
  vector2d b;

  SetVector2d ( &b,
                ru->x*b1d+rs->x*c1d+ruu->x*b1+rus->x*c1-qvv->x*f1-qvt->x*g1,
                ru->y*b1d+rs->y*c1d+ruu->y*b1+rus->y*c1-qvv->y*f1-qvt->y*g1 );
  Solve2x2d ( qv, qt, &b, f1d, g1d );
} /*SolveCompatibilityEq2bd*/

static void SolveCompatibilityEq3d (
            const vector2d *ru, const vector2d *rs, const vector2d *ruu,
            const vector2d *rus, const vector2d *rss, const vector2d *qvv,
            double b1, double c1, double *b2, double *c2 )
{
/* compute values of b2, c2 from the second order partial derivatives */
/* compatibility equation */
  vector2d b;

  SetVector2d ( &b,
                (double)(qvv->x-b1*b1*ruu->x-2.0*b1*c1*rus->x-c1*c1*rss->x),
                (double)(qvv->y-b1*b1*ruu->y-2.0*b1*c1*rus->y-c1*c1*rss->y) );
  Solve2x2d ( ru, rs, &b, b2, c2 );
} /*SolveCompatibilityEq3d*/

static void SolveCompatibilityEq4ad (
            const vector2d *qv, const vector2d *qt, const vector2d *qvv,
            const vector2d *qvt, const vector2d *qvvv, const vector2d *qvvt,
            const vector2d *ru, const vector2d *rs, const vector2d *ruu,
            const vector2d *rus, const vector2d *rss, const vector2d *ruuu,
            const vector2d *ruus, const vector2d *russ,
            double b1, double b1d, double c1, double c1d,
            double f1, double f1d, double f1dd, double g1, double g1d, double g1dd,
            double b2, double *b2d, double c2, double *c2d )
{
/* compute derivatives of b2, c2, from the third order mixed partial */
/* derivatives compatibility equation */
  vector2d b;

  SetVector2d ( &b,
    (double)(f1dd*qv->x+g1dd*qt->x+2.0*(f1d*qvv->x+g1d*qvt->x)+f1*qvvv->x+g1*qvvt->x-
    ((2.0*b1*b1d+b2)*ruu->x+(2.0*(b1d*c1+b1*c1d)+c2)*rus->x+
    2.0*c1*c1d*rss->x+b1*b1*ruuu->x+2.0*b1*c1*ruus->x+c1*c1*russ->x)),
    (double)(f1dd*qv->y+g1dd*qt->y+2.0*(f1d*qvv->y+g1d*qvt->y)+f1*qvvv->y+g1*qvvt->y-
    ((2.0*b1*b1d+b2)*ruu->y+(2.0*(b1d*c1+b1*c1d)+c2)*rus->y+
    2.0*c1*c1d*rss->y+b1*b1*ruuu->y+2.0*b1*c1*ruus->y+c1*c1*russ->y)) );
  Solve2x2d ( ru, rs, &b, b2d, c2d );
} /*SolveCompatibilityEq4ad*/

static void SolveCompatibilityEq4bd (
            const vector2d *qv, const vector2d *qt, const vector2d *qvv,
            const vector2d *qvt, const vector2d *qvvt,
            const vector2d *ru, const vector2d *rs,
            const vector2d *rus, const vector2d *rss, const vector2d *russ,
            double b1d, double c1, double c1d, double f1d, double f1dd,
            double g1, double g1d, double g1dd,
            double *b2d, double c2, double *c2d )
{
/* compute derivatives of b2, c2, from the third order mixed partial */
/* derivatives compatibility equation; this version assumes that */
/* b1 = f1 = b2 = 0 */
  vector2d b;

  SetVector2d ( &b,
    (double)(f1dd*qv->x+g1dd*qt->x+2.0*(f1d*qvv->x+g1d*qvt->x)+g1*qvvt->x-
    ((2.0*b1d*c1+c2)*rus->x+2.0*c1*c1d*rss->x+c1*c1*russ->x)),
    (double)(f1dd*qv->y+g1dd*qt->y+2.0*(f1d*qvv->y+g1d*qvt->y)+g1*qvvt->y-
    ((2.0*b1d*c1+c2)*rus->y+2.0*c1*c1d*rss->y+c1*c1*russ->y)) );
  Solve2x2d ( ru, rs, &b, b2d, c2d );
} /*SolveCompatibilityEq4bd*/

static void SolveCompatibilityEq4cd (
            const vector2d *qv, const vector2d *qt, const vector2d *qvv,
            const vector2d *qvt, const vector2d *qvvt,
            const vector2d *ru, const vector2d *rs,
            const vector2d *rus, const vector2d *rss, const vector2d *russ,
            double b1d, double c1, double c1d, double f1d, double *f1dd,
            double g1, double g1d, double *g1dd,
            double b2d, double c2, double c2d )
{
/* compute second order derivatives of f1, g1, from the third order */
/* mixed partial derivatives compatibility equation; this version assumes */
/* that b1 = f1 = b2 = 0 */
  vector2d b;

  SetVector2d ( &b,
    (double)(-2.0*f1d*qvv->x -2.0*g1d*qvt->x -g1*qvvt->x +b2d*ru->x +c2d*rs->x
    +(2.0*b1d*c1+c2)*rus->x +2.0*c1*c1d*rss->x +c1*c1*russ->x),
    (double)(-2.0*f1d*qvv->y -2.0*g1d*qvt->y -g1*qvvt->y +b2d*ru->y +c2d*rs->y
    +(2.0*b1d*c1+c2)*rus->y +2.0*c1*c1d*rss->y +c1*c1*russ->y) );
  Solve2x2d ( qv, qt, &b, f1dd, g1dd );
} /*SolveCompatibilityEq4cd*/

static void SolveCompatibilityEq5ad (
            const vector2d *qv, const vector2d *qt, const vector2d *qvv,
            const vector2d *qvt, const vector2d *qtt, const vector2d *qvvv,
            const vector2d *qvvt, const vector2d *qvtt, const vector2d *qvvvv,
            const vector2d *qvvvt, const vector2d *qvvtt,
            const vector2d *ru, const vector2d *rs, const vector2d *ruu,
            const vector2d *rus, const vector2d *rss, const vector2d *ruuu,
            const vector2d *ruus, const vector2d *russ, const vector2d *ruuuu,
            const vector2d *ruuus, const vector2d *ruuss,
            double b1, double b1d, double b1dd, double c1, double c1d, double c1dd,
            double f1, double f1d, double f1dd, double g1, double g1d, double g1dd,
            double b2, double b2d, double *b2dd, double c2, double c2d, double c2dd,
            double f2, double f2d, double *f2dd, double g2, double g2d, double g2dd )
{
/* compute second order derivatives of b2, f2, from the fourth order mixed */
/* partial derivatives compatibility equation */
  vector2d b;

  SetVector2d ( &b,
     (double)(g2dd*qt->x+2.0*(f1*f1dd+f1d*f1d+f2d)*qvv->x+
     2.0*(f1dd*g1+2.0*f1d*g1d+f1*g1dd+g2d)*qvt->x+
     2.0*(g1*g1dd+g1d*g1d)*qtt->x+(4.0*f1*f1d+f2)*qvvv->x+
     (4.0*(f1d*g1+f1*g1d)+g2)*qvvt->x+4.0*g1*g1d*qvtt->x+
     f1*f1*qvvvv->x+2.0*f1*g1*qvvvt->x+g1*g1*qvvtt->x -
     (c2dd*rs->x+2.0*(b1*b1dd+b1d*b1d+b2d)*ruu->x+
     2.0*(b1dd*c1+2.0*b1d*c1d+b1*c1dd+c2d)*rus->x+
     2.0*(c1*c1dd+c1d*c1d)*rss->x+(4.0*b1*b1d+b2)*ruuu->x+
     (4.0*(b1d*c1+b1*c1d)+c2)*ruus->x+4.0*c1*c1d*russ->x+
     b1*b1*ruuuu->x+2.0*b1*c1*ruuus->x+c1*c1*ruuss->x)),
     (double)(g2dd*qt->y+2.0*(f1*f1dd+f1d*f1d+f2d)*qvv->y+
     2.0*(f1dd*g1+2.0*f1d*g1d+f1*g1dd+g2d)*qvt->y+
     2.0*(g1*g1dd+g1d*g1d)*qtt->y+(4.0*f1*f1d+f2)*qvvv->y+
     (4.0*(f1d*g1+f1*g1d)+g2)*qvvt->y+4.0*g1*g1d*qvtt->y+
     f1*f1*qvvvv->y+2.0*f1*g1*qvvvt->y+g1*g1*qvvtt->y -
     (c2dd*rs->y+2.0*(b1*b1dd+b1d*b1d+b2d)*ruu->y+
     2.0*(b1dd*c1+2.0*b1d*c1d+b1*c1dd+c2d)*rus->y+
     2.0*(c1*c1dd+c1d*c1d)*rss->y+(4.0*b1*b1d+b2)*ruuu->y+
     (4.0*(b1d*c1+b1*c1d)+c2)*ruus->y+4.0*c1*c1d*russ->y+
     b1*b1*ruuuu->y+2.0*b1*c1*ruuus->y+c1*c1*ruuss->y)) );
  Solve2x2d ( ru, qv, &b, b2dd, f2dd );
  *f2dd = -(*f2dd);
} /*SolveCompatibilityEq5ad*/

static void SolveCompatibilityEq5bd (
            const vector2d *qv, const vector2d *qt, const vector2d *qvv,
            const vector2d *qvt, const vector2d *qtt,
            const vector2d *qvvt, const vector2d *qvtt, const vector2d *qvvtt,
            const vector2d *ru, const vector2d *rs, const vector2d *ruu,
            const vector2d *rus, const vector2d *rss,
            const vector2d *ruus, const vector2d *russ, const vector2d *ruuss,
            double b1d, double b1dd, double c1, double c1d, double c1dd,
            double f1d, double f1dd, double g1, double g1d, double g1dd,
            double b2d, double *b2dd, double c2, double c2d, double *c2dd,
            double f2d, double f2dd, double g2, double g2d, double g2dd )
{
/* compute second order derivatives of b2, c2, from the fourth order mixed */
/* partial derivatives compatibility equation; this version assumes that */
/* b1 = f1 = b2 = f2 = 0 */
  vector2d b;

  SetVector2d ( &b,
     (double)(f2dd*qv->x+g2dd*qt->x+2.0*(f1d*f1d+f2d)*qvv->x+
     2.0*(f1dd*g1+2.0*f1d*g1d+g2d)*qvt->x+2.0*(g1*g1dd+g1d*g1d)*qtt->x+
     (4.0*f1d*g1+g2)*qvvt->x+4.0*g1*g1d*qvtt->x+g1*g1*qvvtt->x -
     (2.0*(b1d*b1d+b2d)*ruu->x+2.0*(b1dd*c1+2.0*b1d*c1d+c2d)*rus->x+
     2.0*(c1*c1dd+c1d*c1d)*rss->x+(4.0*b1d*c1+c2)*ruus->x+
     4.0*c1*c1d*russ->x+c1*c1*ruuss->x)),
     (double)(f2dd*qv->y+g2dd*qt->y+2.0*(f1d*f1d+f2d)*qvv->y+
     2.0*(f1dd*g1+2.0*f1d*g1d+g2d)*qvt->y+2.0*(g1*g1dd+g1d*g1d)*qtt->y+
     (4.0*f1d*g1+g2)*qvvt->y+4.0*g1*g1d*qvtt->y+g1*g1*qvvtt->y -
     (2.0*(b1d*b1d+b2d)*ruu->y+2.0*(b1dd*c1+2.0*b1d*c1d+c2d)*rus->y+
     2.0*(c1*c1dd+c1d*c1d)*rss->y+(4.0*b1d*c1+c2)*ruus->y+
     4.0*c1*c1d*russ->y+c1*c1*ruuss->y)) );
  Solve2x2d ( ru, rs, &b, b2dd, c2dd );
} /*SolveCompatibilityEq5bd*/

static void SolveCompatibilityEq5cd (
            const vector2d *qv, const vector2d *qt, const vector2d *qvv,
            const vector2d *qvt, const vector2d *qtt, const vector2d *qvvt,
            const vector2d *qvtt, const vector2d *qvvtt,
            const vector2d *ru, const vector2d *rs, const vector2d *ruu,
            const vector2d *rus, const vector2d *rss, const vector2d *ruus,
            const vector2d *russ, const vector2d *ruuss,
            double b1d, double b1dd, double c1, double c1d, double c1dd,
            double f1d, double f1dd, double g1, double g1d, double g1dd,
            double b2d, double *b2dd, double c2, double c2d, double c2dd,
            double f2d, double *f2dd, double g2, double g2d, double g2dd )
{
/* compute second order derivatives of b2, f2, from the fourth order mixed */
/* partial derivatives compatibility equation; this version assumes that */
/* b1 = f1 = b2 = f2 = 0 */
  vector2d b;

  SetVector2d ( &b,
     (double)(g2dd*qt->x+2.0*(f1d*f1d+f2d)*qvv->x+2.0*(f1dd*g1+2.0*f1d*g1d+g2d)*qvt->x+
     2.0*(g1*g1dd+g1d*g1d)*qtt->x+(4.0*f1d*g1+g2)*qvvt->x+4.0*g1*g1d*qvtt->x+
     g1*g1*qvvtt->x -
     (c2dd*rs->x+2.0*(b1d*b1d+b2d)*ruu->x+2.0*(b1dd*c1+2.0*b1d*c1d+c2d)*rus->x+
     2.0*(c1*c1dd+c1d*c1d)*rss->x+(4.0*b1d*c1+c2)*ruus->x+4.0*c1*c1d*russ->x+
     c1*c1*ruuss->x)),
     (double)(g2dd*qt->y+2.0*(f1d*f1d+f2d)*qvv->y+2.0*(f1dd*g1+2.0*f1d*g1d+g2d)*qvt->y+
     2.0*(g1*g1dd+g1d*g1d)*qtt->y+(4.0*f1d*g1+g2)*qvvt->y+4.0*g1*g1d*qvtt->y+
     g1*g1*qvvtt->y -
     (c2dd*rs->y+2.0*(b1d*b1d+b2d)*ruu->y+2.0*(b1dd*c1+2.0*b1d*c1d+c2d)*rus->y+
     2.0*(c1*c1dd+c1d*c1d)*rss->y+(4.0*b1d*c1+c2)*ruus->y+4.0*c1*c1d*russ->y+
     c1*c1*ruuss->y)) );
  Solve2x2d ( ru, qv, &b, b2dd, f2dd );
  *f2dd = -(*f2dd);
} /*SolveCompatibilityEq5cd*/

