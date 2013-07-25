
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2008                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

static boolean FindCentralDomPointd ( GHoleDomaind *domain, point2d *centre )
{
  void     *sp;
  GHolePrivateRecd *privateG;
  G1HolePrivateRecd *privateG1;
  int      i, hole_k;
  point2d  c;
  vector2d v;
  point2d  *surrpcbc;
  int      option;
  int      ndata, *idata;
  double   *fdata, *a, *b;

  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;

  option = privateG1->GetOption ( domain, G1HQUERY_CENTRAL_POINT, 0,
                                  &ndata, &idata, &fdata );
  privateG1->opt_cp = (char)option;
  switch ( option ) {
case G1H_DEFAULT:
    surrpcbc = privateG->surrpcbc;
    c = surrpcbc[0];
    for ( i = 1; i < hole_k; i++ )
      AddVector2d ( &c, &surrpcbc[36*i], &c );
    MultVector2d ( 1.0/(double)hole_k, &c, centre );
    return true;

case G1H_CENTRAL_POINT_ALT:
    sp = pkv_GetScratchMemTop ();
    a = pkv_GetScratchMemd ( 3*hole_k );
    if ( !a ) {
      domain->error_code = G1H_ERROR_NO_SCRATCH_MEMORY;
      return false;
    }
    b = &a[2*hole_k];
    surrpcbc = privateG->surrpcbc;   
    for ( i = 0; i < hole_k; i++ ) {
      c = surrpcbc[36*i];
      v = surrpcbc[(2*i+1)*18+1];   
      SetVector2d ( &v, -v.y, v.x );
      NormalizeVector2d ( &v );
      memcpy ( &a[2*i], &v, 2*sizeof(double) );
      b[i] = (double)DotProduct2d ( &v, &c );
    }
    pkn_multiSolveRLSQd ( hole_k, 2, a, 1, 1, b, 1, (double*)centre );
    pkv_SetScratchMemTop ( sp );
    return true;

case G1H_CENTRAL_POINT_GIVEN:
    if ( ndata == 2 ) {
      memcpy ( centre, fdata, sizeof(point2d) );
      return true;
    }
    else
      goto failure;
      
default:
failure:
    domain->error_code = G1H_ERROR_INVALID_OPTION;
    return false;
  }
} /*FindCentralDomPointd*/

static void CorrectThePartitiond ( int hole_k, int pitch, vector2d *part )
{
/* the purpose of this procedure is to ensure that the angle between the */
/* partition halflines is either PI or different from PI by at least 3deg. */
/* this procedure is not quite correct, but it should work in practical cases */
#define ATOL (3.0*PI/180.0)
  int      i, j;
  vector2d *a, *b, c;
  double   la, lb, sa, ca, ang, x, y;

  for ( i = 0, a = part;  i < hole_k-1;  i++, a += pitch ) {
    la = (double)sqrt ( DotProduct2d ( a, a ) );
    for ( j = i+1, b = a+pitch;  j < hole_k;  j++, b += pitch ) {
      lb = (double)sqrt ( DotProduct2d ( b, b ) );
      ca = (double)DotProduct2d ( a, b ) / (la*lb);
      sa = (double)det2d ( a, b ) / (la*lb);
      ang = (double)atan2 ( sa, ca );
      if ( fabs(ang-PI) < ATOL || fabs(ang+PI) < ATOL ) {
        SubtractPoints2d ( a, b, &c );
/*
printf ( "ang = %d\n", ang );
printf ( "a = (%f,%f), b = (%f,%f), c = (%f, %f)\n", a->x, a->y, b->x, b->y, c.x, c.y );
*/
        x = (double)DotProduct2d ( &c, &c );
        y = (double)DotProduct2d ( a, &c );
        MultVector2d ( y/x, &c, a );
        y = (double)DotProduct2d ( b, &c );
        MultVector2d ( y/x, &c, b );
/*
printf ( "a = (%f,%f), b = (%f,%f)\n", a->x, a->y, b->x, b->y );
*/
      }
    }
  }
#undef ATOL
} /*CorrectThePartitiond*/

static boolean FindAuxDPatchesd ( GHoleDomaind *domain )
{
  void     *sp;
  GHolePrivateRecd *privateG;
  G1HolePrivateRecd *privateG1;
  int      hole_k;
  vector2d *surrpcbc;
  vector2d *omcbc, *omcbcd, *omc, *omcd;
  int      i, j;
  double   *hole_knots;
  double   g11, h11;
  int      omcsize;
  int      option, ndata, *idata;
  double   *fdata;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  hole_knots = domain->hole_knots;
  privateG = domain->privateG;
  privateG1 = domain->privateG1;
  omcsize = (2*G1H_OMCDEG+1)*hole_k*sizeof(vector2d);
  privateG1->omcbc = omcbc = (vector2d*)malloc ( omcsize );
  privateG1->omc = omc = (vector2d*)malloc ( omcsize );
  if ( !omcbc || !omc )
    goto failure;
  omcbcd = &omcbc[(G1H_OMCDEG+1)*hole_k];
  omcd = &omc[(G1H_OMCDEG+1)*hole_k];

        /* initialize the boundary conditions to zero */
  memset ( omcbc, 0, omcsize );

        /* setup the nonzero boundary conditions */
  surrpcbc = privateG->surrpcbc;
  if ( !FindCentralDomPointd ( domain, &omcbc[0] ) )
    goto failure;

  for ( i = 0; i < hole_k; i++ )
    omcbc[i*(G1H_OMCDEG+1)] = omcbc[0];

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    g11 = (hole_knots[11*j+5]-hole_knots[11*j+4])/
          (hole_knots[11*j+4]-hole_knots[11*j+3]);
    h11 = (double)((hole_knots[11*i+6]-hole_knots[11*i+4])/
          (2.0*(hole_knots[11*i+6]-hole_knots[11*i+5])));
    omcbc[i*(G1H_OMCDEG+1)+3] = surrpcbc[36*i];
    MultVector2d ( g11, &surrpcbc[(2*i+1)*18+1], &omcbc[i*(G1H_OMCDEG+1)+4] );
/*    omcbcd[i*G1H_OMCDEG+2] = surrpcbc[(2*i+1)*18+3]; */
    MultVector2d ( h11, &surrpcbc[(2*i+1)*18+3], &omcbcd[i*G1H_OMCDEG+2] );
    MultVector2d ( g11*h11, &surrpcbc[(2*i+1)*18+4], &omcbcd[i*G1H_OMCDEG+3] );
  }

  option = privateG1->GetOption ( domain, G1HQUERY_CENTRAL_DERIVATIVES1, 0,
                                &ndata, &idata, &fdata );
  privateG1->opt_der1 = (char)option;
  switch ( option ) {
case G1H_DEFAULT:
    for ( i = 0; i < hole_k; i++ ) {
      SubtractPoints2d ( &omcbc[i*(G1H_OMCDEG+1)+3], &omcbc[i*(G1H_OMCDEG+1)], &omcbc[i*(G1H_OMCDEG+1)+1] );
      omcbcd[i*G1H_OMCDEG] = omcbcd[i*G1H_OMCDEG+2];
    }
    break;

case G1H_CENRTAL_DERIVATIVES1_ALT:
    for ( i = 0; i < hole_k; i++ ) {
      SubtractPoints2d ( &omcbc[i*(G1H_OMCDEG+1)+3], &omcbc[i*(G1H_OMCDEG+1)], &omcbc[i*(G1H_OMCDEG+1)+1] );
      SetVector2d ( &omcbcd[i*G1H_OMCDEG], -omcbc[i*(G1H_OMCDEG+1)+1].y, omcbc[i*(G1H_OMCDEG+1)+1].x );
    }
    break;

case G1H_CENTRAL_DERIVATIVES1_GIVEN:
        /* ******** no data correctness verified at the moment ******** */
    for ( i = 0; i < hole_k; i++ ) {
      memcpy ( &omcbc[i*(G1H_OMCDEG+1)+1], &fdata[2*i], sizeof(vector2d) );
      omcbcd[i*G1H_OMCDEG] = omcbcd[i*G1H_OMCDEG+2];
    }
    break;

default:
    domain->error_code = G1H_ERROR_INVALID_OPTION;
    goto failure;
  }

  CorrectThePartitiond ( hole_k, (G1H_OMCDEG+1), &omcbc[1] );

  option = privateG1->GetOption ( domain, G1HQUERY_DOMAIN_CURVES, 0,
                                &ndata, &idata, &fdata );
  switch ( option ) {
case G1H_DEFAULT:
    mbs_multiInterp2knHermiteBezd ( hole_k, 2, G1_CROSS00DEG,
         3, 10, (double*)omcbc, 2, 10, (double*)&omcbc[3], 10, (double*)omc );
    mbs_multiInterp2knHermiteBezd ( hole_k, 2, G1_CROSS00DEG-1,
         2, 8, (double*)omcbcd, 2, 8, (double*)&omcbcd[2], 8, (double*)omcd );
    break;

default:
    domain->error_code = G1H_ERROR_INVALID_OPTION;
    goto failure;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindAuxDPatchesd*/

