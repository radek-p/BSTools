
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

void SetupRefBox ( float x0, float x1, float y0, float y1, float z0, float z1 );
void SetupParProj ( float aspect );
void ResetCPos ( void );
void InitPerspProj ( float aspect );
void UpdatePerspProj ( void );
void ResizeWinContents ( void );
void RedrawGeomWindows ( void );
void DisplayHoleK ( xge_widget *er, boolean onscreen );

void RysujOkno ( xge_widget *er, boolean onscreen );
void SetParZoom ( float zf );
boolean OknoMsg ( xge_widget *er, int msg, int key, short x, short y );

void RysujPOkno ( xge_widget *er, boolean onscreen );
void SetZoom ( float zf );
boolean POknoMsg ( xge_widget *er, int msg, int key, short x, short y );

void init_edwin ( void );
void destroy_edwin ( void );

int CallBack ( xge_widget *er, int msg, int key, short x, short y );

