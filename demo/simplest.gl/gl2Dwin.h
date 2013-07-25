
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

void RysujOkno ( xge_widget *er, boolean onscreen );
int NearestPointFound ( CameraRecf *CPos, short x, short y );
void MoveCurrentPoint ( xge_widget *er, short x, short y );
void SelectPoints ( CameraRecf *CPos, short x0, short x1, short y0, short y1 );
void UnSelectPoints ( CameraRecf *CPos, short x0, short x1, short y0, short y1 );
void TransformPoints ( void );
int CallBack ( xge_widget *er, int msg, int key, short x, short y );
void init_edwin ( void );
void destroy_edwin ( void );

