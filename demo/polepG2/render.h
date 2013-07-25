
extern boolean RenderingIsOn;
extern XImage  *theimage;

/* rendering switches */
extern boolean rswGaussian, rswMean, rswSections,
               rswVDepRefl, rswVIndRefl, rswChess,
               eswSections, eswLines1, eswLines2, eswLight,
               eswEdRendering;

extern point3f RendPoints[10];

void RedrawRendPoints ( int id );
int FindNearestRendPoint ( int id, int x, int y );
void SetRendPoint ( int id, int np, int x, int y );
void ResetRendLines ( void );

boolean InitRenderer ( void );
void DestroyRenderer ( void );
boolean ResetRenderer ( void );
boolean RendEnterBezierPatch ( int n, int m, point3f *cp );
boolean BeginRendering ( ed_rect *er );
void StopRendering ( void );
void RenderLine ( void );

