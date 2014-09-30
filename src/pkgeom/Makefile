
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

OBJECTS = geom2f.o geom3f.o geom4f.o \
  quaternionf.o convhf.o convcone2f.o raylinef.o lbclipf.o \
  geom2d.o geom3d.o geom4d.o \
  quaterniond.o convhd.o convcone2d.o raylined.o lbclipd.o \
  geom3df.o geom3fd.o

TARGET = $(LIBDIR)/libpkgeom.a

default: $(TARGET)

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJECTS)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

geom2f.o: geom2f.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

geom3f.o: geom3f.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

geom3df.o: geom3df.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

geom3fd.o: geom3fd.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

geom4f.o: geom4f.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

quaternionf.o: quaternionf.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

geom2d.o: geom2d.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

geom3d.o: geom3d.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

geom4d.o: geom4d.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

quaterniond.o: quaterniond.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

convhf.o: convhf.c $(INCDIR)/msgpool.h $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h \
  $(INCDIR)/pkgeomd.h $(INCDIR)/convh.h

convhd.o: convhd.c $(INCDIR)/msgpool.h $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h \
  $(INCDIR)/pkgeomd.h $(INCDIR)/convh.h

convcone2f.o: convcone2f.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

convcone2d.o: convcone2d.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

raylinef.o: raylinef.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

raylined.o: raylined.c $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h $(INCDIR)/pkgeomd.h

lbclipf.o: lbclipf.c $(INCDIR)/msgpool.h $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h \
  $(INCDIR)/pkgeomd.h $(INCDIR)/pkgeomclip.h

lbclipd.o: lbclipd.c $(INCDIR)/msgpool.h $(INCDIR)/pkgeom.h $(INCDIR)/pkgeomf.h \
  $(INCDIR)/pkgeomd.h $(INCDIR)/pkgeomclip.h

clean:
	rm -f $(OBJECTS)

mrproper:
	rm -f $(OBJECTS)
	rm -f $(TARGET)
