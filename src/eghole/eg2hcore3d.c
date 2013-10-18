
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2005, 2013                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */

static boolean FindCentralDomPointd ( GHoleDomaind *domain, point2d *centre )
{
  void     *sp;
  GHolePrivateRecd *privateG;
  G2HolePrivateRecd *privateG2;
  int      i, hole_k;
  point2d  c;
  vector2d v;
  point2d *surrpcbc;
  int      option;
  int      ndata, *idata;
  double   *fdata, *a, *b;

  hole_k = domain->hole_k;
  privateG = domain->privateG;
  privateG2 = domain->privateG2;

  option = privateG2->GetOption ( domain, G2HQUERY_CENTRAL_POINT, 0,
                                  &ndata, &idata, &fdata );
  privateG2->opt_cp = (char)option;
  switch ( option ) {
case G2H_DEFAULT:
    surrpcbc = privateG->surrpcbc;
    c = surrpcbc[0];
    for ( i = 1; i < hole_k; i++ )
      AddVector2d ( &c, &surrpcbc[36*i], &c );
    MultVector2d ( 1.0/(double)hole_k, &c, centre );
    return true;

case G2H_CENTRAL_POINT_ALT:
    sp = pkv_GetScratchMemTop ();
    a = pkv_GetScratchMemd ( 3*hole_k );
    if ( !a ) {
      domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
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

case G2H_CENTRAL_POINT_GIVEN:
    if ( ndata == 2 ) {
      memcpy ( centre, fdata, sizeof(point2d) );
      return true;
    }
    else
      goto failure;
      
default:
failure:
    domain->error_code = G2H_ERROR_INVALID_OPTION;
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
  G2HolePrivateRecd *privateG2;
  int      hole_k;
  vector2d *surrpcbc;
  vector2d *omcbc, *omcbcd, *omcbcdd, *omc, *omcd, *omcdd, *auxomc;
  int      i, j, l;
  double   *hole_knots;
  double   g11, g11p2, h11, h11p2;
  int      omcsize;
  int      option, ndata, *idata;
  double   *fdata;

  sp = pkv_GetScratchMemTop ();

  hole_k = domain->hole_k;
  hole_knots = domain->hole_knots;
  privateG = domain->privateG;
  privateG2 = domain->privateG2;
  omcsize = 3*G2H_OMCDEG*hole_k*sizeof(vector2d);
  privateG2->omcbc = omcbc = (vector2d*)malloc ( omcsize );
  privateG2->omc = omc = (vector2d*)malloc ( omcsize );
  if ( !omcbc || !omc )
    goto failure;
  omcbcd = &omcbc[(G2H_OMCDEG+1)*hole_k];
  omcbcdd = &omcbcd[G2H_OMCDEG*hole_k];
  omcd = &omc[(G2H_OMCDEG+1)*hole_k];
  omcdd = &omcd[G2H_OMCDEG*hole_k];

        /* initialize the boundary conditions to zero */
  memset ( omcbc, 0, omcsize );

        /* setup the nonzero boundary conditions */
  surrpcbc = privateG->surrpcbc;
  if ( !FindCentralDomPointd ( domain, &omcbc[0] ) )
    goto failure;

  for ( i = 0; i < hole_k; i++ )
    omcbc[i*(G2H_OMCDEG+1)] = omcbc[0];

  for ( i = 0; i < hole_k; i++ ) {
    j = (i+1) % hole_k;
    g11 = (hole_knots[11*j+5]-hole_knots[11*j+4])/
          (hole_knots[11*j+4]-hole_knots[11*j+3]);
    g11p2 = g11*g11;
    h11 = (double)((hole_knots[11*i+6]-hole_knots[11*i+4])/
                  (2.0*(hole_knots[11*i+6]-hole_knots[11*i+5])));
    h11p2 = h11*h11;
    omcbc[i*(G2H_OMCDEG+1)+5] = surrpcbc[36*i];
    MultVector2d ( g11, &surrpcbc[(2*i+1)*18+1], &omcbc[i*(G2H_OMCDEG+1)+6] );
    MultVector2d ( g11p2, &surrpcbc[(2*i+1)*18+2], &omcbc[i*(G2H_OMCDEG+1)+7] );
/*    omcbcd[i*G2H_OMCDEG+4] = surrpcbc[(2*i+1)*18+3]; */
    MultVector2d ( h11, &surrpcbc[(2*i+1)*18+3], &omcbcd[i*G2H_OMCDEG+4] );
    MultVector2d ( g11*h11, &surrpcbc[(2*i+1)*18+4], &omcbcd[i*G2H_OMCDEG+5] );
    MultVector2d ( g11p2*h11, &surrpcbc[(2*i+1)*18+5], &omcbcd[i*G2H_OMCDEG+6] );
/*    omcbcdd[i*(G2H_OMCDEG-1)+3] = surrpcbc[(2*i+1)*18+6]; */
    MultVector2d ( h11p2, &surrpcbc[(2*i+1)*18+6], &omcbcdd[i*(G2H_OMCDEG-1)+3] );
    MultVector2d ( g11*h11p2, &surrpcbc[(2*i+1)*18+7], &omcbcdd[i*(G2H_OMCDEG-1)+4] );
    MultVector2d ( g11p2*h11p2, &surrpcbc[(2*i+1)*18+8], &omcbcdd[i*(G2H_OMCDEG-1)+5] );
  }

  option = privateG2->GetOption ( domain, G2HQUERY_CENTRAL_DERIVATIVES1, 0,
                                &ndata, &idata, &fdata );
  privateG2->opt_der1 = (char)option;
  switch ( option ) {
case G2H_DEFAULT:
    for ( i = 0; i < hole_k; i++ ) {
      SubtractPoints2d ( &omcbc[i*(G2H_OMCDEG+1)+5], &omcbc[i*(G2H_OMCDEG+1)],
                         &omcbc[i*(G2H_OMCDEG+1)+1] );
      omcbcd[i*G2H_OMCDEG] = omcbcd[i*G2H_OMCDEG+4];
    }
    break;

case G2H_CENRTAL_DERIVATIVES1_ALT:
    for ( i = 0; i < hole_k; i++ ) {
      SubtractPoints2d ( &omcbc[i*(G2H_OMCDEG+1)+5], &omcbc[i*(G2H_OMCDEG+1)],
                         &omcbc[i*(G2H_OMCDEG+1)+1] );
      SetVector2d ( &omcbcd[i*G2H_OMCDEG],
                    -omcbc[i*(G2H_OMCDEG+1)+1].y, omcbc[i*(G2H_OMCDEG+1)+1].x );
    }
    break;

case G2H_CENTRAL_DERIVATIVES1_GIVEN:
        /* ******** no data correctness verified at the moment ******** */
    for ( i = 0; i < hole_k; i++ ) {
      memcpy ( &omcbc[i*(G2H_OMCDEG+1)+1], &fdata[2*i], sizeof(vector2d) );
      omcbcd[i*G2H_OMCDEG] = omcbcd[i*G2H_OMCDEG+4];
    }
    break;

default:
    domain->error_code = G2H_ERROR_INVALID_OPTION;
    goto failure;
  }

  CorrectThePartitiond ( hole_k, (G2H_OMCDEG+1), &omcbc[1] );

  option = privateG2->GetOption ( domain, G2HQUERY_DOMAIN_CURVES, 0,
                                &ndata, &idata, &fdata );
  switch ( option ) {
case G2H_DEFAULT:
    if ( !mbs_multiInterp2knHermiteBezd ( hole_k, 2, G2_CROSS00DEG,
             5, 16, (double*)omcbc, 3, 16, (double*)&omcbc[5], 16, (double*)omc ) )
      goto failure;
    if ( !mbs_multiInterp2knHermiteBezd ( hole_k, 2, G2_CROSS00DEG-1,
             4, 14, (double*)omcbcd, 3, 14, (double*)&omcbcd[4], 14, (double*)omcd ) )
      goto failure;
    if ( !mbs_multiInterp2knHermiteBezd ( hole_k, 2, G2_CROSS00DEG-2,
             3, 12, (double*)omcbcdd, 3, 12, (double*)&omcbcdd[3], 12, (double*)omcdd ) )
      goto failure;
    break;

case G2H_DOMAIN_CURVES_DEG4:
    auxomc = (vector2d*)pkv_GetScratchMem ( 5*hole_k*sizeof(vector2d) );
    if ( !auxomc ) {
      domain->error_code = G2H_ERROR_NO_SCRATCH_MEMORY;
      goto failure;
    }
    if ( !mbs_multiInterp2knHermiteBezd ( hole_k, 2, 4,
             2, 16, (double*)omcbc, 3, 16, (double*)&omcbc[5], 10, (double*)auxomc ) )
      goto failure;
    if ( !mbs_multiBCDegElevd ( hole_k, 2, 10, 4, (double*)auxomc,
                                G2_CROSS00DEG-4, 16, &ndata, (double*)omc ) )
      goto failure;
        /* derivatives of order 2, 3, 4 at 0 are computed "by hand" */
    for ( i = 0; i < hole_k; i++ ) {
      for ( j = 1; j <= 4; j++ )
        for ( l = 4; l >= j; l-- )
          SubtractPoints2d ( &auxomc[5*i+l], &auxomc[5*i+l-1], &auxomc[5*i+l] );
      MultVector2d ( 12.0, &auxomc[5*i+2], &omcbc[(G2H_OMCDEG+1)*i+2] );
      MultVector2d ( 24.0, &auxomc[5*i+3], &omcbc[(G2H_OMCDEG+1)*i+3] );
      MultVector2d ( 24.0, &auxomc[5*i+4], &omcbc[(G2H_OMCDEG+1)*i+4] );
    }
        /* similarly, first order cross derivatives */
    if ( !mbs_multiInterp2knHermiteBezd ( hole_k, 2, 3,
             1, 14, (double*)omcbcd, 3, 14, (double*)&omcbcd[4], 10, (double*)auxomc ) )
      goto failure;
    if ( !mbs_multiBCDegElevd ( hole_k, 2, 10, 3, (double*)auxomc,
                                G2_CROSS00DEG-1-3, 14, &ndata, (double*)omcd ) )
      goto failure;
    for ( i = 0; i < hole_k; i++ ) {
      for ( j = 1; j <= 3; j++ )
        for ( l = 3; l >= j; l-- )
          SubtractPoints2d ( &auxomc[5*i+l], &auxomc[5*i+l-1], &auxomc[5*i+l] );
      MultVector2d ( 3.0, &auxomc[5*i+1], &omcbcd[G2H_OMCDEG*i+1] );
      MultVector2d ( 6.0, &auxomc[5*i+2], &omcbcd[G2H_OMCDEG*i+2] );
      MultVector2d ( 6.0, &auxomc[5*i+3], &omcbcd[G2H_OMCDEG*i+3] );
    }

        /* second order cross derivatives */
    if ( !mbs_multiInterp2knHermiteBezd ( hole_k, 2, 3,
             1, 12, (double*)omcbcdd, 3, 12, (double*)&omcbcdd[3], 10, (double*)auxomc ) )
      goto failure;
    if ( !mbs_multiBCDegElevd ( hole_k, 2, 10, 3, (double*)auxomc,
                                G2_CROSS00DEG-2-3, 12, &ndata, (double*)omcdd ) )
      goto failure;
    for ( i = 0; i < hole_k; i++ ) {
      for ( j = 1; j <= 2; j++ )
        for ( l = 2; l >= j; l-- )
          SubtractPoints2d ( &auxomc[5*i+l], &auxomc[5*i+l-1], &auxomc[5*i+l] );
      MultVector2d ( 3.0, &auxomc[5*i+1], &omcbcdd[(G2H_OMCDEG-1)*i+1] );
      MultVector2d ( 6.0, &auxomc[5*i+2], &omcbcdd[(G2H_OMCDEG-1)*i+2] );
    }
    break;

default:
    domain->error_code = G2H_ERROR_INVALID_OPTION;
    goto failure;
  }

  pkv_SetScratchMemTop ( sp );
  return true;

failure:
  pkv_SetScratchMemTop ( sp );
  return false;
} /*FindAuxDPatchesd*/

