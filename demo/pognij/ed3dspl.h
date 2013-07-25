
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

#define MAX_PATH_LGT     1024
#define MAX_PATH_SHRT      63
#define MAX_FILENAME_LGT   64

void SetupParProj ( void );
void ResetCPos ( void );
void FindCPos ( void );
void InitPerspProj ( void );
void UpdatePerspProj ( void );
void ResizeWinContents ( void );
void RedrawGeomWindows ( void );

void RysujOkno ( xge_widget *er, boolean onscreen );
void PanParWindows ( xge_widget *er, short x, short y );
void SetParZoom ( float zf );
void ResetCurveMapping ( void );
void FindCurveMapping ( void );
boolean OknoMsg ( xge_widget *er, int msg, int key, short x, short y );

void RysujPOkno ( xge_widget *er, boolean onscreen );
void SetZoom ( float zf );
boolean POknoMsg ( xge_widget *er, int msg, int key, short x, short y );

void RysujWOkno ( xge_widget *er, boolean onscreen );
boolean WOknoMsg ( xge_widget *er, int msg, int key, short x, short y );

void init_edwin ( void );
void destroy_edwin ( void );

int CallBack ( xge_widget *er, int msg, int key, short x, short y );

