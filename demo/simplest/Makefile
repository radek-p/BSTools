
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2014
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib
BINDIR = ../bin

include ../../options.mak

LDFLAGS =

LIBS = -lxgedit -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria \
  -lX11 -lXext -lpthread -lm

GLIBS = -lxgedit -lmultibs -lcamera -lpkgeom -lpknum -lpkvaria \
  -lGL -lX11 -lXext -lpthread -lm

default: $(BINDIR)/simplest $(BINDIR)/widgets $(BINDIR)/ncolours \
  $(BINDIR)/2Dwin $(BINDIR)/3Dwin $(BINDIR)/2Dwind $(BINDIR)/3Dwind \
  $(BINDIR)/scrollwidg $(BINDIR)/ipc $(BINDIR)/ipc_proc \
  $(BINDIR)/qbtest

.c.o:
	$(CC) $(CFLAGS) -c $< -o $*.o

$(BINDIR)/simplest: $(SRCDIR)/simplest.c $(SRCDIR)/simplest.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/simplest.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/simplest

$(BINDIR)/widgets: $(SRCDIR)/widgets.c $(SRCDIR)/widgets.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/widgets.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/widgets
	rm -f $(SRCDIR)/widgets
	ln -s $(BINDIR)/widgets $(SRCDIR)/widgets

$(BINDIR)/ncolours: $(SRCDIR)/ncolours.c $(SRCDIR)/ncolours.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/ncolours.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/ncolours
	rm -f $(SRCDIR)/ncolours
	ln -s $(BINDIR)/ncolours $(SRCDIR)/ncolours

$(BINDIR)/2Dwin: $(SRCDIR)/2Dwin.c $(SRCDIR)/2Dwin.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/2Dwin.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/2Dwin

$(BINDIR)/2Dwind: $(SRCDIR)/2Dwind.c $(SRCDIR)/2Dwind.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/2Dwind.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/2Dwind

$(BINDIR)/3Dwin: $(SRCDIR)/3Dwin.c $(SRCDIR)/3Dwin.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/3Dwin.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/3Dwin

$(BINDIR)/3Dwind: $(SRCDIR)/3Dwind.c $(SRCDIR)/3Dwind.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/3Dwind.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/3Dwind

$(BINDIR)/scrollwidg: $(SRCDIR)/scrollwidg.c $(SRCDIR)/scrollwidg.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/scrollwidg.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/scrollwidg

$(BINDIR)/glsimplest: $(SRCDIR)/glsimplest.c $(SRCDIR)/glsimplest.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/glsimplest.c \
  $(LIBDIR) $(GLIBS) -o $(BINDIR)/glsimplest

$(BINDIR)/ipc: $(SRCDIR)/ipc.c $(SRCDIR)/ipc.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/ipc.c \
  $(LIBDIR) $(GLIBS) -o $(BINDIR)/ipc

$(BINDIR)/ipc_proc: $(SRCDIR)/ipc_proc.c $(SRCDIR)/ipc_proc.h
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/ipc_proc.c \
  $(LIBDIR) $(GLIBS) -o $(BINDIR)/ipc_proc

$(BINDIR)/qbtest: $(SRCDIR)/qbtest.c
	$(CC) $(LDFLAGS) $(CFLAGS) -I$(INCDIR) $(SRCDIR)/qbtest.c \
  $(LIBDIR) $(LIBS) -o $(BINDIR)/qbtest
	rm -f $(SRCDIR)/qbtest
	ln -s $(BINDIR)/qbtest $(SRCDIR)/qbtest

clean:
	rm -f *.o

mrproper:
	rm -f *.o
	rm -f $(BINDIR)/simplest $(BINDIR)/widgets $(BINDIR)/ncolours \
  $(BINDIR)/2Dwin $(BINDIR)/3Dwin $(BINDIR)/2Dwind $(BINDIR)/3Dwind \
  $(BINDIR)/scrollwidg $(BINDIR)/glsimplest \
  $(BINDIR)/ipc $(BINDIR)/ipc_proc $(BINDIR)/qbtest

