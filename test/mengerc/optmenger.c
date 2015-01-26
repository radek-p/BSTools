
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>

#include "pkvaria.h"
#include "pkvthreads.h"
#include "pkgeom.h"
#include "pknum.h"
#include "multibs.h"
#include "camera.h"
#include "bsfile.h"
#include "mengerc.h"

#define MAXLKNOT 1000
#define MAXITER  200


boolean      read_in;
int          deg;
int          lkn;
double       kn[MAXLKNOT+1];
point3d      cp[MAXLKNOT];
unsigned int mkcp[MAXLKNOT];
int          ncpu;

char    infilename[100], outfilename[100], logfilename[100],
        parfilename[100];
int     opt;

/* default values */
double w = 4.0;
double penalty[MENGERC_NPPARAM] =
  {1.420666814e+06, 1.106542879e+07, 1.063609655e+04,
   1.219877706e+05, 5.771308492e+03};

/* ////////////////////////////////////////////////////////////////////////// */
FILE *out;

void ReadPenaltyParameters ( char *fn, double *penalty )
{
  FILE *f;
  int  i;

  f = fopen ( fn, "r+" );
  if ( !f )
    exit ( 1 );
  for ( i = 0; i < MENGERC_NPPARAM; i++ )
    fscanf ( f, "%lf", &penalty[i] );
  fclose ( f );
} /*ReadPenaltyParameters*/

void PrintHelpMsg ( void )
{
  printf ( "Usage:\n\n" );
  printf ( "  optmenger -i input_file.bs [-e exponent] [-p param_file.txt] [-o?] [-l log_file.txt] [-o output_file.bs]\n\n" );
  printf ( "exponent is a real number greater than 3\n");
  printf ( "-o1 -o2 -o3 specifies method of choosing penalty parameters\n\n" );
  exit ( 0 );
} /*PrintHelpMsg*/

void ReadCMDParameters ( int argc, char **argv )
{
  boolean inf, outf, logf, expo, paramf;
  int     i;

  inf = outf = logf = expo = paramf = false;
  opt = MENGERC_OPT_NONE;
  for ( i = 0; i < argc; i++ ) {
    if ( !strcmp ( argv[i], "-help" ) ) {
      PrintHelpMsg ();
    }
    else if ( !strcmp ( argv[i], "-i" ) ) {  /* input file name */
      i ++;
      strncpy ( infilename, argv[i], 99 );
      inf = true;
    }
    else if ( !strcmp ( argv[i], "-o" ) ) { /* output file name */
      i ++;
      strncpy ( outfilename, argv[i], 99 );
      outf = true;
    }
    else if ( !strcmp ( argv[i], "-l" ) ) { /* log file name */
      i ++;
      strncpy ( logfilename, argv[i], 99 );
      logf = true;
    }
    else if ( !strcmp ( argv[i], "-p" ) ) { /* penalty parameters file name */
      i ++;
      strncpy ( parfilename, argv[i], 99 );
      paramf = true;
      ReadPenaltyParameters ( parfilename, penalty );
    }
    else if ( !strcmp (argv[i], "-e" ) ) { /* Menger curvature exponent */
      i ++;
      if ( sscanf ( argv[i], "%lf", &w ) != 1 )
        exit ( 1 );
      expo = true;
    }
    else if ( !strcmp ( argv[i], "-o1" ) )
      opt = MENGERC_OPT_FULL1;
    else if ( !strcmp ( argv[i], "-o2" ) )
      opt = MENGERC_OPT_FULL2;
    else if ( !strcmp ( argv[i], "-o3" ) )
      opt = MENGERC_OPT_PART;
  }
  if ( !inf ) {
    PrintHelpMsg ();
  }
  if ( !outf ) {
    sprintf ( outfilename, "wyn%s", infilename );
  }
  if ( !logf ) {
    sprintf ( logfilename, "wyn%sg", infilename );
    i = strlen ( logfilename ) - 3;
    logfilename[i] = 'l';
    logfilename[i+1] = 'o';
  }
} /*ReadCMDParameters*/

boolean WriteMKCP ( void *usrdata )
{
  mengerc_data *md;

  md = (mengerc_data*)usrdata;
  memset ( mkcp, 0, (lkn-deg)*sizeof(unsigned int) );
  mkcp[md->mdi] = 9;
  if ( md->mdi < deg )
    mkcp[md->mdi+lkn-2*deg] = 9;
  bsf_WritePointsMK ( lkn-deg, mkcp );
  return true;
} /*WriteMKCP*/

void ReadTheCurve ( void *usrdata, const char *name, int ident,
                    int degree, int lastknot,
                    const double *knots, boolean closed,
                    const point4d *cpoints,
                    int spdimen, boolean rational )
{
  int i, ncp;

  if ( !read_in ) {
    if ( degree < 3 ) {
      printf ( "the curve degree must be at least 3\n" );
      return;
    }
    deg = degree;
    lkn = lastknot;
    memcpy ( kn, knots, (lkn+1)*sizeof(double) );
    ncp = lkn-deg;
    for ( i = 0; i < ncp; i++ )
      Point4to3d ( &cpoints[i], &cp[i] );
    read_in = true;
  }
} /*ReadTheCurve*/

typedef struct {
    char    *outfilename;
    char    *logfilename;
    double  *penalty_param;
    int     deg, lkn;
    double  *knots;
    point3d *cpoints;
  } mydata;
  
boolean WriteCAttr ( void *usrdata )
{
  mydata *md;

  md = usrdata;
  bsf_WritePointsMK ( md->lkn-md->deg, mkcp );
  return true;
} /*WriteCAttr*/

void OutIter ( void *usrdata, boolean ppopt, int mdi,
               int it, int itres, double f, double gn )
{
  mydata *md;
  char   name[100], code;
  FILE   *out;
  int    i;

  if ( !usrdata )
    return;

  md = (mydata*)usrdata;
  sprintf ( name, "c%02d", it );
  bsf_OpenOutputFile ( md->outfilename, it > 0 );
  for ( i = 0; i < md->lkn-md->deg; i++ )
    mkcp[i] = 1;
  mkcp[mdi] = 9;
  bsf_WriteBSplineCurved ( 3, 3, false,
                           md->deg, md->lkn, md->knots, true, &md->cpoints[0].x,
                           name, it, WriteCAttr, usrdata );
  bsf_CloseOutputFile ();
  switch ( itres ) {
case PKN_LMT_ERROR:            code = 'e';  break;
case PKN_LMT_CROSSED_LIMIT:    code = 'l';  break;
case PKN_LMT_CONTINUE_N:       code = '+';  break;
case PKN_LMT_CONTINUE_LM:      code = '-';  break;
case PKN_LMT_CONTINUE_LM_P:    code = '#';  break;
case PKN_LMT_FOUND_MINIMUM:    code = '.';  break;
case PKN_LMT_FOUND_ZEROGRAD:   code = 'z';  break;
case PKN_LMT_FOUND_ZEROGRAD_P: code = 'Z';  break;
case PKN_LMT_FOUND_BARRIER:    code = 'b';  break;
case PKN_LMT_NO_PROGRESS:      code = 's';  break;
default:                       code = ' ';  break;
  }
  out = fopen ( md->logfilename, "a+" );
  if ( out ) {
    if ( ppopt ) {
      fprintf ( out, "penalty parameters:\n" );
      for ( i = 0; i < MENGERC_NPPARAM; i++ )
        fprintf ( out, "%14.8e ", md->penalty_param[i] );
      fprintf ( out, "\n" );
    }
    fprintf ( out, "%3d (%c): %14.8g %14.8g\n", it, code, f, gn );
    fclose ( out );
  }
  if ( ppopt ) {
    printf ( "penalty parameters:\n" );
    for ( i = 0; i < MENGERC_NPPARAM; i++ )
      printf ( "%14.8e ", md->penalty_param[i] );
    printf ( "\n" );
  }
  printf ( "%3d (%c): %14.8g %14.8g\n", it, code, f, gn );
} /*OutIter*/

boolean _OptimizeMengerCurvature (
                      int deg, int lkn, double *knots, point3d *cpoints,
                      double w, double penalty_param[MENGERC_NPPARAM],
                      int nqkn, int npthr, int opt, int maxit,
                      void (*outiter)(void *usrdata,
                                      boolean ppopt, int mdi,
                                      int itres, int it, double f, double g),
                      void *usrdata )
{
  FILE         *out;
  mengerc_data md;
  int          i;
  boolean      finished;
  double       f;

  if ( !mengerc_InitMCOptimization ( deg, lkn, knots, cpoints, w,
                                     penalty_param, nqkn, npthr, opt, &md ) )
    return false;
  if ( outiter )
    outiter ( usrdata, md.ppopt, md.mdi, 0, -2, md.lastf, md.gn );
  finished = false;
  for ( i = 0;  i < maxit;  i++ ) {  
    if ( !mengerc_IterMCOptimization ( &md, &finished ) )
      goto failure;
    if ( outiter )
      outiter ( usrdata, md.ppopt, md.mdi, i+1, md.itres, md.lastf, md.gn );
    if ( finished )
      break;
  }
  if ( finished ) {
    mengerc_IntegralMengerf ( 3*(lkn-deg), &md, (double*)cpoints, &f );
    out = fopen ( logfilename, "a+" );
    fprintf ( out, "terms:\n" );
    fprintf ( out, "Kp = %14.8g, Gp = %14.8g, L = %14.8g\n",
                md.ffkM, md.ffkMe, md.lgt );
    for ( i = 0; i < 5; i++ )
      fprintf ( out, "%14.8g\n", md.ffR[i] );
    fclose ( out );
  }
  mengerc_UntieTheCurve ( &md );
  return true;  

failure:
  mengerc_UntieTheCurve ( &md );   
  return false;
} /*_OptimizeMengerCurvature*/

int main ( int argc, char **argv )
{
  char            hostname[100];
  bsf_UserReaders rr;
  char            cwd[256];
  mydata          md;
  int             time;
  int             i;

  pkv_InitScratchMem ( 160*1048576 );
  ncpu = pkv_FindNCPU ();
  if ( !pkv_InitPThreads ( ncpu ) ) {
    printf ( "cannot initialize pthreads\n" );
    exit ( 1 );
  }
  ReadCMDParameters ( argc, argv );

  hostname[0] = 0;
  getcwd ( cwd, 255 );
  if ( !chdir ( "/etc" ) ) {
    out = fopen ( "HOSTNAME", "r" );
    if ( out ) {
      fscanf ( out, "%s", hostname );
      fclose ( out );
      for ( i = 0; hostname[i]; i++ )
        if ( hostname[i] == '.' ) {
          hostname[i] = 0;
          break;
        }
      printf ( "host %s: ncpu = %d\n", hostname, ncpu );
    }
    chdir ( cwd );
  }

  bsf_ClearReaders ( &rr );
  read_in = false;
  bsf_BSC4ReadFuncd ( &rr, ReadTheCurve, 10, MAXLKNOT );
  printf ( "%s\n", infilename );
  if ( !bsf_ReadBSFiled ( infilename, &rr ) ) {
    printf ( "Something wrong with this file\n" );
    exit ( 1 );
  }
  if ( !read_in ) {
    printf ( "no cure in it\n" );
    exit ( 1 );
  }

  out = fopen ( logfilename, "w+" );
  if ( !out )
    exit ( 1 );
  setvbuf ( out, NULL, _IONBF, 0 );
  if ( hostname[0] )
    fprintf ( out, "host %s: ", hostname );
  fprintf ( out, "ncpu = %d\n", ncpu );
  fprintf ( out, "input file: %s\n", infilename );
  fprintf ( out, "exponent = %5.2f\n", w );
  fclose ( out );

  md.outfilename = outfilename;
  md.logfilename = logfilename;
  md.penalty_param = penalty;
  md.deg = deg;
  md.lkn = lkn;
  md.knots = kn;
  md.cpoints = cp;
  pkv_Tic ( NULL );
  _OptimizeMengerCurvature ( deg, lkn, kn, cp, w, penalty, 3,
                             ncpu, opt, MAXITER, OutIter, (void*)&md );
  time = pkv_Toc ( NULL );
  printf ( "%s\n", outfilename );
  if ( (out = fopen ( logfilename, "a+" )) ) {
    fprintf ( out, "time = %6.2f\n", pkv_Seconds ( time ) );
    fprintf ( out, "penalty parameters:\n" );
    for ( i = 0; i < MENGERC_NPPARAM; i++ )
      fprintf ( out, "%10.6e\n", penalty[i] );
    fclose ( out );
  }

  pkv_DestroyPThreads ();
  pkv_DestroyScratchMem ();
  printf ( "\n" );
  exit ( 0 );
} /*main*/

