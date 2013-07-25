
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_PATH_LGT     1024
#define MAX_PATH_SHRT      63
#define MAX_FILENAME_LGT   64


void RysujOkno ( xge_widget *er, boolean onscreen );
void ResetCurveMapping ( void );
void FindCurveMapping ( void );

void RysujWOkno ( xge_widget *er, boolean onscreen );
void RedrawCurveAndKnotsArea ( void );
boolean WOknoMsg ( xge_widget *er, int msg, int key, short x, short y );

void init_edwin ( void );
void destroy_edwin ( void );
int CallBack ( xge_widget *er, int msg, int key, short x, short y );

