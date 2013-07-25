
/* //////////////////////////////////////////////////// */
/* This file is a part of the BSTools procedure package */
/* written by Przemyslaw Kiciak.                        */
/* //////////////////////////////////////////////////// */

extern pid_t parent_pid, child_pid;
extern int pipe_in[2], pipe_out[2];

/* the procedures of the parent process */
boolean MakeTheChild ( const char *name, int magic );
void CallTheChild ( void );
void SignalTheChild ( void );
void ProcessChildCommand ( void );

/* the procedures of the child process */
void SignalTheParent ( void );
void ChildLoop ( void );
boolean ChildInitCommunication ( int argc, char **argv, int magic );
