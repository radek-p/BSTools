
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
LIBDIR = ../../lib

include ../../options.mak

INCLUDES = $(INCDIR)/pkvaria.h $(INCDIR)/pkvthreads.h \
  $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h \
  $(INCDIR)/multibs.h $(INCDIR)/multibsf.h $(INCDIR)/multibsd.h \
  $(INCDIR)/raybez.h $(INCDIR)/raybezf.h $(INCDIR)/raybezd.h \
  $(SRCDIR)/raybezprivatef.h $(SRCDIR)/raybezprivated.h

OBJECTS = rayboxf.o rayrboxf.o raybezf.o rayrbezf.o raybeztreef.o \
  raybezcoffsetf.o rayrbezcoffsetf.o rayrbeztreef.o rayboxtestf.o \
  rayboxd.o rayrboxd.o raybezd.o rayrbezd.o raybeztreed.o \
  raybezcoffsetd.o rayrbezcoffsetd.o rayrbeztreed.o rayboxtestd.o \
  raypthreads.o rbezintersf.o rbezintersd.o \
  raybezconvhtest2f.o raybezconvhtest2d.o \
  rbshomotopyf.o rbshomotopyd.o \
  raybezphpointsf.o raybezphpointsd.o \
  raybezconvhtest1f.o raybezconvhtest1d.o \
  raybezchpointsf.o raybezchpointsd.o \
  rayrbezchpointsf.o rayrbezchpointsd.o \
  rayrbezphpointsf.o rayrbezphpointsd.o

TARGET = $(LIBDIR)/libraybez.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

raypthreads.o: raypthreads.c $(INCLUDES)

rayboxf.o: rayboxf.c $(INCLUDES)
rayrboxf.o: rayrboxf.c $(INCLUDES)
rayboxtestf.o: rayboxtestf.c $(INCLUDES)
raybezconvhtest1f.o: raybezconvhtest1f.c $(INCLUDES)
raybezconvhtest2f.o: raybezconvhtest2f.c $(INCLUDES)
raybezf.o: raybezf.c $(INCLUDES)
rayrbezf.o: rayrbezf.c $(INCLUDES)
raybeztreef.o: raybeztreef.c $(INCLUDES)
raybezcoffsetf.o: raybezcoffsetf.c $(INCLUDES)
rayrbezcoffsetf.o: rayrbezcoffsetf.c $(INCLUDES)
rayrbeztreef.o: rayrbeztreef.c $(INCLUDES)
rbezintersf.o: rbezintersf.c $(INCLUDES)
rbshomotopyf.o: rbshomotopyf.c $(INCLUDES)
raybezphpointsf.o: raybezphpointsf.c $(INCLUDES)
raybezchpointsf.o: raybezchpointsf.c $(INCLUDES)
rayrbezchpointsf.o: rayrbezchpointsf.c $(INCLUDES)
rayrbezphpointsf.o: rayrbezphpointsf.c $(INCLUDES)

rayboxd.o: rayboxd.c $(INCLUDES)
rayrboxd.o: rayrboxd.c $(INCLUDES)
rayboxtestd.o: rayboxtestd.c $(INCLUDES)
raybezconvhtest1d.o: raybezconvhtest1d.c $(INCLUDES)
raybezconvhtest2d.o: raybezconvhtest2d.c $(INCLUDES)
raybezd.o: raybezd.c $(INCLUDES)
rayrbezd.o: rayrbezd.c $(INCLUDES)
raybeztreed.o: raybeztreed.c $(INCLUDES)
raybezcoffsetd.o: raybezcoffsetd.c $(INCLUDES)
rayrbezcoffsetd.o: rayrbezcoffsetd.c $(INCLUDES)
rayrbeztreed.o: rayrbeztreed.c $(INCLUDES)
rbezintersd.o: rbezintersd.c $(INCLUDES)
rbshomotopyd.o: rbshomotopyd.c $(INCLUDES)
raybezphpointsd.o: raybezphpointsd.c $(INCLUDES)
raybezchpointsd.o: raybezchpointsd.c $(INCLUDES)
rayrbezchpointsd.o: rayrbezchpointsd.c $(INCLUDES)
rayrbezphpointsd.o: rayrbezphpointsd.c $(INCLUDES)

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)
