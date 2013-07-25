
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2009
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME    = /home/przemek
SRCDIR  = .
INCDIR  = ../../include
LIBDIR  = -L/usr/X11R6/lib -L../../lib
BINDIR = ../bin

include ../../options.mak

OBJECTS = $(SRCDIR)/oldxgedit.o $(SRCDIR)/g1ekernel.o $(SRCDIR)/edg1hole.o \
  $(SRCDIR)/datagenf.o $(SRCDIR)/edg1hole1.o $(SRCDIR)/edg1hole2.o \
  $(SRCDIR)/edg1hole3.o $(SRCDIR)/drawbez.o $(SRCDIR)/g1info.o \
  $(SRCDIR)/render.o

LIBS = -leghole -lraybez -lmultibs -lcamera -lpsout -lpkgeom -lpknum \
  -lpkvaria -lm -lX11

default: $(BINDIR)/polepG1

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

$(SRCDIR)/oldxgedit.o: $(SRCDIR)/oldxgedit.c $(SRCDIR)/oldxgedit.h

$(SRCDIR)/datagenf.o: $(SRCDIR)/datagenf.c $(SRCDIR)/datagenf.h

$(SRCDIR)/g1ekernel.o: $(SRCDIR)/g1ekernel.c $(SRCDIR)/g1ekernel.h \
  $(SRCDIR)/drawbez.h $(SRCDIR)/render.h

$(SRCDIR)/edg1hole.o: $(SRCDIR)/edg1hole.c $(SRCDIR)/edg1hole.h \
  $(SRCDIR)/g1ekernel.h $(SRCDIR)/render.h

$(SRCDIR)/edg1hole1.o: $(SRCDIR)/edg1hole1.c $(SRCDIR)/edg1hole.h \
  $(SRCDIR)/g1ekernel.h $(SRCDIR)/render.h

$(SRCDIR)/edg1hole2.o: $(SRCDIR)/edg1hole2.c $(SRCDIR)/edg1hole.h \
  $(SRCDIR)/g1ekernel.h $(SRCDIR)/render.h

$(SRCDIR)/edg1hole3.o: $(SRCDIR)/edg1hole3.c $(SRCDIR)/edg1hole.h \
  $(SRCDIR)/g1ekernel.h $(SRCDIR)/render.h

$(SRCDIR)/drawbez.o: $(SRCDIR)/drawbez.c $(SRCDIR)/drawbez.h

$(SRCDIR)/g1info.o: $(SRCDIR)/g1info.c $(SRCDIR)/g1ekernel.h

$(SRCDIR)/render.o: $(SRCDIR)/render.c $(SRCDIR)/render.h $(SRCDIR)/g1ekernel.h

$(BINDIR)/polepG1: $(OBJECTS)
	$(CC) $(CFLAGS) -I$(INCDIR) $(OBJECTS) $(LIBDIR) $(LIBS) \
  -o $(BINDIR)/polepG1

clean:
	rm -f $(SRCDIR)/*.o

mrproper:
	rm -f $(SRCDIR)/*.o
	rm -f $(BINDIR)/polepG1

