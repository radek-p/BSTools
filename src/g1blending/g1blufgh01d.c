
/* ///////////////////////////////////////////////////////////////////////// */
/* This file is a part of the BSTools package                                */
/* written by Przemyslaw Kiciak                                              */
/* ///////////////////////////////////////////////////////////////////////// */
/* (C) Copyright by Przemyslaw Kiciak, 2010, 2012                            */
/* this package is distributed under the terms of the                        */
/* Lesser GNU Public License, see the file COPYING.LIB                       */
/* ///////////////////////////////////////////////////////////////////////// */
/* this file was written by Mateusz Markowski                                */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pknum.h"
#include "pkgeom.h"
#include "multibs.h"
#include "g1blendingd.h"

#include "g1blprivated.h"

/* /////////////// obliczanie całego funkcjonału dla JEDNEGO KWADRATU  //////////////////////////
obliczanie całki w jednym kwadracie */
void g1bl_UFuncSQd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv,
                     int pitch, point3d *cp,
                     double tC,
                     int isq, int jsq,
                     double *ftab )
{
  int      fcpn; /* first control point number */
  vector3d pder[5], lapP;                       /* parameterization derivatives */
  double   Gstar[3], detG;  /* first form and related data */
  double   Bstar[3], tH;       /* second form and related data */
  double   W, X, f, g; /* f funkcjonał podstawowy, g - funkcjonał jednoznaczności parametryzacji*/
  int      i, j;
  double   sum0, sum1;

  fcpn = isq*pitch + jsq; /* który kwadrat*/
  sum0 = 0.0;
  for ( i = 0; i < nkn; i++ ) {
    sum1 = 0.0;
    for ( j = 0; j < nkn; j++ ) {
        /* compute the parameterization derivatives */
      _g1bl_UCompPDerd ( nkn, Nitab, pitch, cp, fcpn, i, j, pder );
        /* compute the first fundamental form and its derivatives */
      _g1bl_UCompGStard ( pder, Gstar );
        /* compute the second fundamental form and its derivatives */
      _g1bl_UCompBStard ( pder, Bstar );
        /* the first functional */
      
      detG = g11*g22 - g12*g12;
      X = 1.0 / (detG * detG *sqrt(detG));
      tH = g11*tb22 + tb11*g22 - 2.0*g12*tb12;
      W = (tH*tH)/4.0;
      f = X*W;
      SetVector3d(&lapP, 0.0, 0.0, 0.0);
      AddVector3d(&pvv, &puu, &lapP);
      
      /* testowanie krzywizny średniej */
      /*{
	if ( i==1 && j == 2){
	  double H, m;
	  m = 2*detG * sqrt(detG);
	  H = tH / m;
	  printf(" Krzywizna liczona przez g1bl_UFuncSQd %f \n", H);
	}
      }*/

        /* the second functional */
          /* scalar product of laplacian parametrization p */
      g = DotProduct3d ( &lapP, &lapP);
          /* second norm of the projection of the Laplacian (scalar)*/
          /* on the normal vector */
      g -= (tb11 + tb22)*(tb11 + tb22)/detG;

        /* add the quadrature term */
      sum1 += (f + tC*g)*qcoeff[j];
    }
    sum0 += sum1*qcoeff[i];
  }
        /* store the integral over the square in the array */
  ftab[0] = sum0;
} /*g2bl_UFuncSQd*/


/* oblicza wartość funkcjonału Q(p) = int ||P(delta p)||_2^2*/
void g1bl_QFuncSQd ( int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv,
                     int pitch, point3d *cp,
                     double tC,
                     int isq, int jsq,
                     double *ftab )
{
  int      fcpn; /* first control point number */
  vector3d pder[5], lapP;                       /* parameterization derivatives */
  double   Gstar[3], detG;  /* first form and related data */
  double   Bstar[3];       /* second form and related data */
  double   g; /* f funkcjonał podstawowy, g - funkcjonał jednoznaczności parametryzacji*/
  int      i, j;
  double   sum0, sum1;

  fcpn = isq*pitch + jsq; /* który kwadrat*/
  sum0 = 0.0;
  for ( i = 0; i < nkn; i++ ) {
    sum1 = 0.0;
    for ( j = 0; j < nkn; j++ ) {
        /* compute the parameterization derivatives */
      _g1bl_UCompPDerd ( nkn, Nitab, pitch, cp, fcpn, i, j, pder );
        /* compute the first fundamental form and its derivatives */
      _g1bl_UCompGStard ( pder, Gstar );
        /* compute the second fundamental form and its derivatives */
      _g1bl_UCompBStard ( pder, Bstar );
        /* the first functional */
      
      detG = g11*g22 - g12*g12;
      SetVector3d(&lapP, 0.0, 0.0, 0.0);
      AddVector3d(&pvv, &puu, &lapP);

        /* the second functional */
          /* scalar product of laplacian parametrization p */
      g = DotProduct3d ( &lapP, &lapP);
          /* second norm of the projection of the Laplacian (scalar)*/
          /* on the normal vector */
      g -= (tb11 + tb22)*(tb11 + tb22)/detG;

        /* add the quadrature term */
      sum1 += g*qcoeff[j];
    }
    sum0 += sum1*qcoeff[i];
  }
        /* store the integral over the square in the array */
  ftab[0] = sum0;
} /*g1bl_QFuncSQd*/

/*oblicza  wartość funkcjonału B(p) = int || delta p||_2^2 w kwadracie (isq jsq)*/
void g1bl_biharmFuncSQd(int nkn, const double *qcoeff, double *Nitab,
                     int lastknotu, int lastknotv,
                     int pitch, point3d *cp,
                     double tC,
                     int isq, int jsq,
                     double *ftab){
    int      fcpn; /* first control point number */
  vector3d pder[5], lapP;                       /* parameterization derivatives */
  double   g; /* f funkcjonał podstawowy, g - funkcjonał jednoznaczności parametryzacji*/
  int      i, j;
  double   sum0, sum1;

  fcpn = isq*pitch + jsq; /* który kwadrat*/
  sum0 = 0.0;
  for ( i = 0; i < nkn; i++ ) {
    sum1 = 0.0;
    for ( j = 0; j < nkn; j++ ) {
        /* compute the parameterization derivatives */
      _g1bl_UCompPDerd ( nkn, Nitab, pitch, cp, fcpn, i, j, pder );
        /* compute the first fundamental form and its derivatives */
      /*_g1bl_UCompGStard ( pder, Gstar );*/
        /* compute the second fundamental form and its derivatives */
      /*_g1bl_UCompBStard ( pder, Bstar );*/
      
      SetVector3d(&lapP, 0.0, 0.0, 0.0);
      AddVector3d(&pvv, &puu, &lapP);
        /* the second functional */
          /* scalar product of laplacian parametrization p */
      g = DotProduct3d ( &lapP, &lapP);
        /* add the quadrature term */
      sum1 +=g*qcoeff[j];
    }
    sum0 += sum1*qcoeff[i];
  }
        /* store the integral over the square in the array */
  ftab[0] = sum0;
} /* g1bl_biharmFuncSQd */


/* //////////////////////////////// funkcjonał + gradient////////////////
Omega = (2,N-2)x(2,M-2)
lastknotu - N
lastknotv - M
isq - numer kwadratu w kierunku osi u 
jsq - numer kwadratu w kierunku osi v
ip0 - ograniczenia na punkty brzegowe żeby nie liczyć poochodnych po brzegu
ip1 -
jp0 - 
jp1 - 
///////////////////////// */
void g1bl_UFuncGradSQd ( int nkn, const double *qcoeff, double *Nitab,
                         int lastknotu, int lastknotv,
                         int pitch, point3d *cp,
                         double tC,
                         int isq, int jsq, int ip0, int ip1, int jp0, int jp1,
                         double *ftab, double *gtab )
{
  int      fcpn;
  vector3d pder[5], lapP;
  double   Gstar[3], DGstar[3*3*9], Bstar[3], DBstar[3*3*9];
  double   detG, tH;
  double   DdetG[3], DtH[3];
  double   W, X, Y, Z;
  double   DW[3], DX[3], DY[3], DZ[3], DlapP[3];
  double   *DGk, *DBk;
  double   lapNi, *Ni;
  double   sum0, sum1, f, g;
  double   Df[3*9], Dg[3*9], sumD1[3*9];
  int      i, j, ii, jj, kk, kk3, ip, jp/*, kp*/;

  fcpn = isq*pitch + jsq;
  memset ( Df, 0, SQUAREGDS*sizeof(double) );
  memset ( Dg, 0, SQUAREGDS*sizeof(double) );
  sum0 = 0.0;
  memset ( gtab, 0, SQUAREGDS*sizeof(double) );
  for ( i = 0; i < nkn; i++ ) {
    sum1 = 0.0;
    memset ( sumD1, 0, SQUAREGDS*sizeof(double) );
    for ( j = 0; j < nkn; j++ ) {
        /* compute the parameterization derivatives */
      _g1bl_UCompPDerd ( nkn, Nitab, pitch, cp, fcpn, i, j, pder );
        /* compute the first fundamental form and its derivatives */
      _g1bl_UCompGStard ( pder, Gstar );
      _g1bl_UCompDGStard ( nkn, Nitab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                           isq, jsq, i, j, pder, DGstar );
        /* compute the second fundamental form and its derivatives */
      _g1bl_UCompBStard ( pder, Bstar );
      _g1bl_UCompDBStard ( nkn, Nitab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                           isq, jsq, i, j, pder, DBstar );
        /* the first functional - liczony jest funkcjonał - kod z pierwszej funkcji */
	
      detG = g11*g22 - g12*g12;
      X = 1.0 / (detG * detG *sqrt(detG));
      tH = g11*tb22 + tb11*g22 - 2.0*g12*tb12;
      W = (tH*tH)/4.0;
      f = X*W;
      SetVector3d(&lapP, 0,0,0);
      AddVector3d(&pvv, &puu, &lapP);

        /* the second functional */
          /* scalar product of laplacian parametrization p */
      g = DotProduct3d ( &lapP, &lapP);
          /* second norm of the projection of the Laplacian (scalar)*/
          /* on the normal vector */
      Y = 1/detG;
      Z = (tb11 + tb22)*(tb11 + tb22);
      g -= Z*Y;
      
        /* add the quadrature term */
      sum1 += (f + tC*g)*qcoeff[j];

        /* compute the functional derivatives */
      for ( ii = 0, ip = isq;  ii <  3;  ii++, ip++ )
        if ( ip >= ip0 && ip < ip1 )
          for ( jj = 0, jp = jsq;  jj < 3;  jj++, jp++ )
            if ( jp >= jp0 && jp < jp1 ) {
              kk = 3*ii + jj;
              kk3 = 3*kk;
              /*kp = ip*pitch + jp;*/
          /* compute derivatives of W, X with respect */
          /* to the coordinates of the point cp[kp] */
              DGk = &DGstar[3*3*kk];
              DBk = &DBstar[3*3*kk];
            /* derivatives of det G */
              pkn_MVectorLinCombd ( 3, 3, DdetG,
                      &DGk[0], g22, &DGk[6], g11, &DGk[3], -2.0*g12 );
            /* derivatives of tH */
              pkn_MVectorLinCombd ( 6, 3, DtH,
                      &DGk[0], tb22, &DBk[6], g11, &DBk[0], g22,
                      &DGk[6], tb11, &DGk[3], -2.0*tb12, &DBk[3], -2.0*g12 );
            
            /* derivatives of X */
              MultVector3d ( -2.5*X/detG, (vector3d*)DdetG, (vector3d *) &DX );
            /* derivatives of W */
              MultVector3d ( 0.5*tH, (vector3d*)DtH, (vector3d *) &DW );
          /* compute derivatives of f */
              pkn_MVectorLinCombd ( 2, 3, &Df[kk3], &DX, W, &DW, X );
	      
          /* compute derivatives of g */
	  /* derivatives of Y */
	    MultVector3d(-Y*Y, (vector3d *) DdetG, (vector3d *) DY);
	    /* derivatives of Z */
	    AddVector3d((vector3d*)&DBk[0], (vector3d*)&DBk[6], (vector3d *)DZ);
	    MultVector3d(2*(tb11 + tb22) , (vector3d *)DZ, (vector3d *)DZ);
	    
            /* derivative of square of the second norm of Laplacian */
              Ni = &Nitab[((kk*nkn + i)*nkn + j)*5];
              lapNi = Ni20 + Ni02;
	      MultVector3d(2*lapNi, &lapP, (vector3d *) DlapP);
              
            /* derivative of the projection on the surface normal */
              AddVector3Md ( (vector3d*)DlapP, (vector3d*)DY, -Z, (vector3d*)&Dg[kk3] );
              AddVector3Md ( (vector3d*)&Dg[kk3], (vector3d*)DZ, -Y,
                             (vector3d*)&Dg[kk3] );

            }
        /* integrate the functional derivatives */
      pkn_MatrixLinCombd ( 1, 3*9, 0, Df, 1.0, 0, Dg, tC, 0, Df );
      pkn_AddMatrixMd ( 1, 3*9, 0, sumD1, 0, Df, qcoeff[j], 0, sumD1 );
    }
    sum0 += sum1*qcoeff[i];
    pkn_AddMatrixMd ( 1, 3*9, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
  }
  *ftab = sum0;
} /*g1bl_UFuncGradSQd*/

/* //////////////////////////////funkcjonał + gradient + hesjan w JEDNYM KWADRACIE JEDNOSTKOWYM!!! /////////////////////// */
void g1bl_UFuncGradHessianSQd ( int nkn, const double *qcoeff, double *Nitab, 
                                double *Nijtab, double *Mijtab,
                                int lastknotu, int lastknotv,
                                int pitch, point3d *cp,
                                double tC,
                                int isq, int jsq, int ip0, int ip1, int jp0, int jp1,
                                double *ftab, double *gtab, double *htab)
{
    int      fcpn;
  vector3d pder[5], lapP;
  double   Gstar[3], DGstar[3*3*9], DDGstar[9*45];
  double   Bstar[3], DBstar[3*3*9], DDBstar[9*45];
  /* u mnie 9 *10 /2 - tyle jest pochodnych 2giego rzedu w tablicy jest diagonala i to co pod - hesjan jest symetryczny*/
  double   detG, DdetG[3*9], DDdetG[9];
  double   tH, DtH[3*9], DDtH[9];
  double   W, X, Y, Z, tb1122, detGpow25, detGpow35, detGpow45, detGpow3, detGpow2;
  double   DW[3*9], DDW[9], DX[3*9], DDX[9], DY[3*9], DDY[9], DZ[3*9], DDZ[9], DlapP[3*9];
  double   lapNi, *Ni;
  double   sum0, sum1, f, g, a;
  double   Df[3*9], Dg[3*9], sumD1[3*9], aux[3];
  int      i, j, ddi;
  
  double   DDf[9], DDg[9], sumDD1[9*5*9], DDaux[9];
  double DBa1122[3], DBb1122[3];
  
  double   *DGa, *DBa, *DGb, *DBb, *gstar, *bstar;
  double   aa, bb;
  double   *Na, *Nb;
  int      ia, ja, ka, ka3, iap, jap/*, kap*/;
  int      ib, jb, kb, kb3, ibp, jbp/*, kbp*/;

  fcpn = isq*pitch + jsq;
  memset ( Df, 0, SQUAREGDS*sizeof(double) );
  memset ( Dg, 0, SQUAREGDS*sizeof(double) );

  sum0 = 0.0;
  memset ( gtab, 0, SQUAREGDS*sizeof(double) );
  memset ( htab, 0, SQUAREHDS*sizeof(double) );
  for ( i = 0; i < nkn; i++ ) {
    sum1 = 0.0;
    memset ( sumD1, 0, SQUAREGDS*sizeof(double) );
    memset ( sumDD1, 0, SQUAREHDS*sizeof(double) );
    for ( j = 0; j < nkn; j++ ) {
        /* compute the parameterization derivatives */
      _g1bl_UCompPDerd ( nkn, Nitab, pitch, cp, fcpn, i, j, pder );
        /* compute the first fundamental form and its derivatives */
      _g1bl_UCompGStard ( pder, Gstar );
      _g1bl_UCompDGStard ( nkn, Nitab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                           isq, jsq, i, j, pder, DGstar );
      _g1bl_UCompDDGStard ( nkn, Nijtab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                            isq, jsq, i, j, DDGstar );
        /* compute the second fundamental form and its derivatives */
      _g1bl_UCompBStard ( pder, Bstar );
      _g1bl_UCompDBStard ( nkn, Nitab, lastknotu, lastknotv, ip0, ip1, jp0, jp1,
                           isq, jsq, i, j, pder, DBstar );
			   
      _g1bl_UCompDDBStard (nkn, Mijtab, lastknotu, lastknotv, ip0, ip1, jp0, jp1, isq, jsq, i, j,
                           pder, DDBstar);
      
        /* the first functional */
	
	        /* the first functional - liczony jest funkcjonał - kod z pierwszej funkcji */
	
      detG = g11*g22 - g12*g12;
      detGpow2 = detG*detG;
      detGpow25 =detGpow2 *sqrt(detG); 
      detGpow3 = detGpow2*detG;
      detGpow35 = detGpow25* detG;
      detGpow45 = detGpow35* detG; 
      
      X = 1.0 /detGpow25;
      tH = g11*tb22 + tb11*g22 - 2.0*g12*tb12;
      W = (tH*tH)/4.0;
      
      f = W/detGpow25;
      SetVector3d(&lapP, 0,0,0);
      AddVector3d(&pvv, &puu, &lapP);

        /* the second functional */
          /* scalar product of laplacian parametrization p */
      g = DotProduct3d ( &lapP, &lapP);
          /* second norm of the projection of the Laplacian (scalar)*/
          /* on the normal vector */
      Y = 1/detG;
      tb1122 = tb11 + tb22;
      Z = tb1122 * tb1122;
      g -= Z/detG;
      
        /* add the quadrature term */
      sum1 += (f + tC*g)*qcoeff[j];
/*  tu kończy się wyliczanie wartości funkcjonału*/    

        /* compute the functional derivatives */
	
	for ( ia = 0, iap = isq;  ia < 3;  ia++, iap++ )
        if ( iap >= ip0 && iap < ip1 )
          for ( ja = 0, jap = jsq;  ja < 3;  ja++, jap++ )
            if ( jap >= jp0 && jap < jp1 ) {
              ka = 3*ia + ja;
              ka3 = 3*ka;
              /*kap = iap*pitch + jap;*/
          /* compute derivatives of W, X, Y, Z with respect */
          /* to the coordinates of the point cp[kap] */
              DGa = &DGstar[3*ka3];
              DBa = &DBstar[3*ka3];
            /* derivatives of det G */
              pkn_MVectorLinCombd ( 3, 3, &DdetG[ka3],
                      &DGa[0], g22, &DGa[6], g11, &DGa[3], -2.0*g12 );

	    /* derivatives of tH */
              pkn_MVectorLinCombd ( 6, 3, &DtH[ka3],
                      &DGa[0], tb22, &DBa[6], g11, &DBa[0], g22,
                      &DGa[6], tb11, &DGa[3], -2.0*tb12, &DBa[3], -2.0*g12 );

            /* derivatives of X */
              MultVector3d ( -2.5/detGpow35, (vector3d*)&DdetG[ka3], (vector3d *) &DX[ka3] ); 
            /* derivatives of W */
              MultVector3d ( 0.5*tH, (vector3d*)&DtH[ka3], (vector3d *) &DW[ka3] );
          /* compute derivatives of f */
              pkn_MVectorLinCombd ( 2, 3, &Df[ka3], &DX[ka3], W, &DW[ka3], X );
	      /*SetVector3d((vector3d*)&Df[ka3], DX[ka3], DX[ka3+1], DX[ka3+2]);*/
	      
          /* COMPUTATION of g DERIVATIVES */
	  /* derivatives of Y */
	    MultVector3d(-1.0/detGpow2, (vector3d *)&DdetG[ka3], (vector3d *) &DY[ka3]);
	    /* derivatives of Z */
	    AddVector3d((vector3d*)&DBa[0], (vector3d*)&DBa[6], (vector3d *)&DZ[ka3]);
	    MultVector3d(2*(tb1122) , (vector3d *)&DZ[ka3], (vector3d *)&DZ[ka3]);
	    
            /* derivative of square of the second norm of Laplacian */
              Ni = &Nitab[((ka*nkn + i)*nkn + j)*5];
              lapNi = Ni20 + Ni02;
	      MultVector3d(2*lapNi, &lapP, (vector3d *) &DlapP[ka3]);
              
            /* derivative of the projection on the surface normal */
	      AddVector3Md ( (vector3d*)&DlapP[ka3], (vector3d*)&DY[ka3], -Z, (vector3d*)&Dg[ka3] );
              AddVector3Md ( (vector3d*)&Dg[ka3], (vector3d*)&DZ[ka3], -Y, (vector3d*)&Dg[ka3] );
	    
	      /*AddVector3Md ( (vector3d*)&DlapP[ka3], (vector3d*)&DZ[ka3], -1.0, (vector3d*)&Dg[ka3] );
	     MultVector3d ( 1.0, (vector3d*)&DlapP[ka3], (vector3d*)&Dg[ka3] );*/
            }
        /* integrate the functional derivatives */
      pkn_MatrixLinCombd ( 1, 3*9, 0, Df, 1.0, 0, Dg, tC, 0, Df );
      pkn_AddMatrixMd ( 1, 3*9, 0, sumD1, 0, Df, qcoeff[j], 0, sumD1 );
	
        /* compute the second order functional derivatives HESJAN */
      for ( ia = 0, iap = isq;  ia < 3;  ia++, iap++ )
        if ( iap >= ip0 && iap < ip1 )
          for ( ja = 0, jap = jsq;  ja < 3;  ja++, jap++ )
            if ( jap >= jp0 && jap < jp1 ) {
              ka = 3*ia + ja;
              ka3 = 3*ka;
              /*kap = iap*pitch + jap;*/
              Na = &Nitab[((ka*nkn + i)*nkn + j)*5];
              DGa = &DGstar[3*ka3];
              DBa = &DBstar[3*ka3];
              for ( ib = 0, ibp = isq;  ib < 3;  ib++, ibp++ )
                if ( ibp >= ip0 && ibp < ip1 )
                  for ( jb = 0, jbp = jsq;  jb < 3;  jb++, jbp++ ) {
                    kb = 3*ib + jb;
                    if ( jbp >= jp0 && jbp < jp1 && kb <= ka ) {
                      kb3 = 3*kb;
                      /*kbp = ibp*pitch + jbp;*/
                      Nb = &Nitab[((kb*nkn + i)*nkn + j)*5];
                      DGb = &DGstar[3*kb3];
                      DBb = &DBstar[3*kb3];

                      gstar = &DDGstar[3*pkn_SymMatIndex(ka,kb)];
                      bstar = &DDBstar[3*3*pkn_SymMatIndex(ka,kb)];

		      /* hesjan detG */
                      DDdetG[0] = DDdetG[4] = DDdetG[8] =
                        g22*gstar[0] + g11*gstar[2] - 2.0*g12*gstar[1];
                      DDdetG[0] +=
                        DGb[0]*DGa[6] + DGb[6]*DGa[0] - 2.0*DGb[3]*DGa[3];
                      DDdetG[1] =
                        DGb[0]*DGa[7] + DGb[6]*DGa[1] - 2.0*DGb[3]*DGa[4];
                      DDdetG[2] =
                        DGb[0]*DGa[8] + DGb[6]*DGa[2] - 2.0*DGb[3]*DGa[5];
                      DDdetG[3] =
                        DGb[1]*DGa[6] + DGb[7]*DGa[0] - 2.0*DGb[4]*DGa[3];
                      DDdetG[4] +=
                        DGb[1]*DGa[7] + DGb[7]*DGa[1] - 2.0*DGb[4]*DGa[4];
                      DDdetG[5] =
                        DGb[1]*DGa[8] + DGb[7]*DGa[2] - 2.0*DGb[4]*DGa[5];
                      DDdetG[6] =
                        DGb[2]*DGa[6] + DGb[8]*DGa[0] - 2.0*DGb[5]*DGa[3];
                      DDdetG[7] =
                        DGb[2]*DGa[7] + DGb[8]*DGa[1] - 2.0*DGb[5]*DGa[4];
                      DDdetG[8] +=
                        DGb[2]*DGa[8] + DGb[8]*DGa[2] - 2.0*DGb[5]*DGa[5];


		      /* hesjan tH */
                      DDtH[0] = DDtH[4] = DDtH[8] =
                        tb22*gstar[0] + tb11*gstar[2] - 2.0*tb12*gstar[1];
                      a = g11*bstar[8] + g22*bstar[2] - 2.0*g12*bstar[5]; 
                      DDtH[1] = a;   DDtH[3] = -a;
                      a = g11*bstar[7] + g22*bstar[1] - 2.0*g12*bstar[4];
                      DDtH[2] = -a;  DDtH[6] = a;
                      a = g11*bstar[6] + g22*bstar[0] - 2.0*g12*bstar[3];
                      DDtH[5] = a;   DDtH[7] = -a;
                      DDtH[0] +=
                        DGb[0]*DBa[6] + DGb[6]*DBa[0] + DBb[0]*DGa[6] + DBb[6]*DGa[0] -
                        2.0*(DGb[3]*DBa[3] + DBb[3]*DGa[3]);
                      DDtH[1] +=
                        DGb[0]*DBa[7] + DGb[6]*DBa[1] + DBb[0]*DGa[7] + DBb[6]*DGa[1] -
                        2.0*(DGb[3]*DBa[4] + DBb[3]*DGa[4]);
                      DDtH[2] +=
                        DGb[0]*DBa[8] + DGb[6]*DBa[2] + DBb[0]*DGa[8] + DBb[6]*DGa[2] -
                        2.0*(DGb[3]*DBa[5] + DBb[3]*DGa[5]);
                      DDtH[3] +=
                        DGb[1]*DBa[6] + DGb[7]*DBa[0] + DBb[1]*DGa[6] + DBb[7]*DGa[0] -
                        2.0*(DGb[4]*DBa[3] + DBb[4]*DGa[3]);
                      DDtH[4] +=
                        DGb[1]*DBa[7] + DGb[7]*DBa[1] + DBb[1]*DGa[7] + DBb[7]*DGa[1] -
                        2.0*(DGb[4]*DBa[4] + DBb[4]*DGa[4]);
                      DDtH[5] +=
                        DGb[1]*DBa[8] + DGb[7]*DBa[2] + DBb[1]*DGa[8] + DBb[7]*DGa[2] -
                        2.0*(DGb[4]*DBa[5] + DBb[4]*DGa[5]);
                      DDtH[6] +=
                        DGb[2]*DBa[6] + DGb[8]*DBa[0] + DBb[2]*DGa[6] + DBb[8]*DGa[0] -
                        2.0*(DGb[5]*DBa[3] + DBb[5]*DGa[3]);
                      DDtH[7] +=
                        DGb[2]*DBa[7] + DGb[8]*DBa[1] + DBb[2]*DGa[7] + DBb[8]*DGa[1] -
                        2.0*(DGb[5]*DBa[4] + DBb[5]*DGa[4]);
                      DDtH[8] +=
                        DGb[2]*DBa[8] + DGb[8]*DBa[2] + DBb[2]*DGa[8] + DBb[8]*DGa[2] -
                        2.0*(DGb[5]*DBa[5] + DBb[5]*DGa[5]);

          /* the second order derivatives of W
	  ka3 = d_i, kb3 = d_j*/
		      /* D_di detG * D_dj^T detG*/
		      memset ( DDf, 0, 9*sizeof(double) );
                      pkn_MultMatrixd ( 3, 1, 1, &DtH[kb3],
                                        3, 0, &DtH[ka3], 3, DDW );
                      pkn_MatrixLinCombd ( 1, 9, 0, DDW, 0.5,
                                           0, DDtH, 0.5 * tH, 0, DDW );
                      
          /* the second order derivatives of X */
		      pkn_MultMatrixd ( 3, 1, 1, &DdetG[kb3],
                                        3, 0, &DdetG[ka3], 3, DDX );
		      pkn_MatrixLinCombd ( 1, 9, 0, DDX, 8.75/detGpow45,
                                           0, DDdetG, -2.5/detGpow35, 0, DDX );
					   
        

          /* the second order derivatives of f */
		      pkn_MatrixLinCombd(1, 9, 0, DDX, W, 0, DDW, X, 0, DDf);
		      pkn_MultMatrixAddd(3, 1, 1, &DX[kb3], 3, 0, &DW[ka3], 3, DDf);
		      pkn_MultMatrixAddd(3, 1, 1, &DW[kb3], 3, 0, &DX[ka3], 3, DDf);
 		      /*pkn_MultMatrixNumd(1,9, 0, DDX, 1.0, 0, DDf);*/
		      
	      /* for historic reasons this is the transposition of the */
              /* necessary matrix, to be fixed later, now just transpose it */
                      pkv_TransposeMatrixd ( 3, 3, 3, DDf, 3, DDaux );
                      memcpy ( DDf, DDaux, 9*sizeof(double) );
		      
          /* compute the second order derivatives of g */
            /* second order derivatives of second norm of the Laplacian */
                      memset ( DDg, 0, 9*sizeof(double) );
                      aa = Na[2] + Na[4];
                      bb = Nb[2] + Nb[4];
                      DDg[0] = DDg[4] = DDg[8] = 2*aa*bb;
		      
            /* second orer derivatives of projection on the surface normal */
	    
	    /* second order derivatives of Y */
		      pkn_MultMatrixd ( 3, 1, 1, &DdetG[kb3],
                                        3, 0, &DdetG[ka3], 3, DDY );
		      pkn_MatrixLinCombd ( 1, 9, 0, DDY, 2.0/detGpow3,
                                           0, DDdetG, -1.0/detGpow2,
                                           0, DDY );
	/* second order derivatives of Z */
		      
                      pkn_AddMatrixd ( 1, 3, 0, &bstar[0], 0, &bstar[2*3], 0, aux );
		      pkn_MultMatrixNumd(1, 3, 0, aux, 2*tb1122, 0, aux);
		      pkn_AddMatrixd( 1,3,0, &DBa[0], 0, &DBa[2*3], 0, DBa1122);
		      pkn_AddMatrixd( 1,3,0, &DBb[0], 0, &DBb[2*3], 0, DBb1122);
		      pkn_MultMatrixd(3, 1, 1, DBb1122, 3, 0, DBa1122, 3, DDZ);
		      pkn_MultMatrixNumd(1, 9, 0, DDZ, 2, 0, DDZ);
		      
		      DDZ[1] += aux[2];  DDZ[3] -= aux[2];
                      DDZ[2] -= aux[1];  DDZ[6] += aux[1];
                      DDZ[5] += aux[0];  DDZ[7] -= aux[0];
	/* second order derivative of g */
	
		     pkn_AddMatrixMd(1, 9, 0, DDg, 0, DDY, -Z, 0, DDg);
		      pkn_AddMatrixMd(1, 9, 0, DDg, 0, DDZ, -Y, 0, DDg);
		      pkn_MultMatrixSubd(3, 1, 1, &DY[kb3], 3, 3, &DZ[ka3], 3, DDg);
		      pkn_MultMatrixSubd(3, 1, 1, &DZ[kb3], 3, 3, &DY[ka3], 3, DDg); 
		      
		     /*pkn_AddMatrixMd(1, 9, 0, DDg, 0, DDZ, -1.0, 0, DDg);
		      pkn_MultMatrixNumd(1, 9, 0, DDY, 1.0, 0, DDg);*/
		      
              /* for historic reasons this is the transposition of the */
              /* necessary matrix, to be fixed later, now just transpose it */
                      pkv_TransposeMatrixd ( 3, 3, 3, DDg, 3, DDaux );
                      memcpy ( DDg, DDaux, 9*sizeof(double) );
		      
            /* integrate */
                      pkn_MatrixLinCombd ( 1, 9, 0, DDf, 1.0, 0, DDg, tC, 0, DDf );
                      ddi = 9*pkn_SymMatIndex ( ka, kb );
                      pkn_AddMatrixMd ( 1, 9, 0, &sumDD1[ddi],
                                        0, DDf, qcoeff[j], 0, &sumDD1[ddi] );
                    }
                  }
            }
    }
    sum0 += sum1*qcoeff[i];
    pkn_AddMatrixMd ( 1, SQUAREGDS, 0, gtab, 0, sumD1, qcoeff[i], 0, gtab );
    pkn_AddMatrixMd ( 1, SQUAREHDS, 0, htab, 0, sumDD1, qcoeff[i], 0, htab );
  }
  ftab[0] = sum0;
} /*g1bl_UFuncGradHessianSQd*/

