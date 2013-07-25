
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2009
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME   = /home/przemek
SRCDIR = .
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib
BINDIR = ../bin

include ../../options.mak

OBJECTS = $(SRCDIR)/oldxgedit.o $(SRCDIR)/g2ekernel.o $(SRCDIR)/edg2hole.o \
  $(SRCDIR)/datagenf.o $(SRCDIR)/edg2hole1.o $(SRCDIR)/edg2hole2.o \
  $(SRCDIR)/edg2hole3.o $(SRCDIR)/drawbez.o $(SRCDIR)/g2info.o \
  $(SRCDIR)/render.o

LIBS = -leghole -lraybez -lmultibs -lcamera -lpsout -lpkgeom -lpknum \
  -lpkvaria -lm -lX11

default: $(BINDIR)/polepG2

.c.o:
	$(CC) -c $(CFLAGS) -I$(INCDIR) $< -o $*.o

$(SRCDIR)/datagenf.o: $(SRCDIR)/datagenf.c $(SRCDIR)/datagenf.h

$(SRCDIR)/g2ekernel.o: $(SRCDIR)/g2ekernel.c $(SRCDIR)/g2ekernel.h \
  $(SRCDIR)/drawbez.h $(SRCDIR)/render.h

$(SRCDIR)/edg2hole.o: $(SRCDIR)/edg2hole.c $(SRCDIR)/edg2hole.h \
  $(SRCDIR)/g2ekernel.h $(SRCDIR)/render.h

$(SRCDIR)/edg2hole1.o: $(SRCDIR)/edg2hole1.c $(SRCDIR)/edg2hole.h \
  $(SRCDIR)/g2ekernel.h $(SRCDIR)/render.h

$(SRCDIR)/edg2hole2.o: $(SRCDIR)/edg2hole2.c $(SRCDIR)/edg2hole.h \
  $(SRCDIR)/g2ekernel.h $(SRCDIR)/render.h

$(SRCDIR)/edg2hole3.o: $(SRCDIR)/edg2hole3.c $(SRCDIR)/edg2hole.h \
  $(SRCDIR)/g2ekernel.h $(SRCDIR)/render.h

$(SRCDIR)/drawbez.o: $(SRCDIR)/drawbez.c $(SRCDIR)/drawbez.h

$(SRCDIR)/g2info.o: $(SRCDIR)/g2info.c $(SRCDIR)/g2ekernel.h

$(SRCDIR)/render.o: $(SRCDIR)/render.c $(SRCDIR)/g2ekernel.h \
  $(SRCDIR)/edg2hole.h $(SRCDIR)/render.h

$(BINDIR)/polepG2: $(OBJECTS)
	$(CC) $(CFLAGS) -I$(INCDIR) $(OBJECTS) $(LIBDIR) $(LIBS) \
  -o $(BINDIR)/polepG2

clean:
	rm -f $(SRCDIR)/*.o

mrproper:
	rm -f $(SRCDIR)/*.o
	rm -f $(BINDIR)/polepG2
