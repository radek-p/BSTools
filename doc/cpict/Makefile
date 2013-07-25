
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2012
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
INCDIR = ../../include
LIBDIR = ../../lib
G1DIR = ../../src/projects

include ../../options.mak

LIBS = -lraybez -lcamera -lbsfile -lbsmesh -leghole -lmultibs \
  -lpsout -lpkgeom -lpknum -lpkvaria -lm

SOURCES = memory.c memshift.c knotrem.c bsapprox.c oslo.c degel.c \
  splmult.c intbsc.c bspdegel.c bspbez.c nvtest.c \
  chknots.c camera.c stereo.c trimpatch.c psout.c raster.c \
  patchder.c patchpdir.c raybezp.c bsder123.c g1demo.c \
  g2hdomain.c datagen.c degred.c datagen.h datagenf.c datagenf.h g2ps.c g2ps.h \
  addspl.c bsplane.c bsmesh.c

PICTURES = memory.ps memshift.ps knotrem.ps bsapprox.ps oslo.ps degel.ps \
  splmult.ps intbsc.ps bspdegel.ps bspbez.ps nvtest1.ps nvtest2.ps \
  chknots.ps camera.ps stereo.ps trimpatch.ps psout.ps rasterbc.ps \
  rasterbs.ps patchder.ps patchpdir.ps raybezp.ps bsder123.ps g1patches1.ps \
  g1patches2.ps g2hdomain.ps degred1.ps degred2.ps addspl.ps bsplane.ps \
  bsmesh.ps

PROGRAMS = memory memshift knotrem bsapprox oslo degel \
  splmult intbsc bspdegel bspbez nvtest \
  chknots camera stereo trimpatch psout raster \
  patchder patchpdir raybezp bsder123 g1demo \
  g2hdomain degred addspl bsplane bsmesh

.SUFFIXES: .ps

default: $(PICTURES)

.c.ps:
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) $< $(LIBS) -o $*
	./$*
	rm -f $*

# the first picture formally depends on all the source files;
# this is done so in order to execute make once in this directory,
# during the TeX processing of the dounmentation.
memory.ps: $(SOURCES)

nvtest1.ps: nvtest.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) nvtest.c $(LIBS) -o nvtest
	./nvtest
	rm -f ./nvtest

nvtest2.ps: nvtest.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) nvtest.c $(LIBS) -o nvtest
	./nvtest
	rm -f ./nvtest

datagen.o: datagen.c datagen.h
	$(CC) -c $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) datagen.c -o datagen.o

g1patches1.ps: g1demo.c datagen.o
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) g1demo.c datagen.o \
  $(G1DIR)/bookg1holef.o $(LIBS) -o g1demo
	./g1demo
	rm -f ./g1demo

g1patches2.ps: g1demo.c datagen.o
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) g1demo.c datagen.o \
  $(G1DIR)/bookg1holef.o $(LIBS) -o g1demo
	./g1demo
	rm -f ./g1demo

rasterbc.ps: raster.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) raster.c $(LIBS) -o raster
	./raster
	rm -f ./raster

rasterbs.ps: raster.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) raster.c $(LIBS) -o raster
	./raster
	rm -f ./raster

g2ps.o: g2ps.c g2ps.h
	$(CC) -c $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) g2ps.c -o g2ps.o

datagenf.o: datagenf.c datagenf.h
	$(CC) -c $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) datagenf.c -o datagenf.o

g2hdomain.ps: g2hdomain.c datagenf.o g2ps.o
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) g2hdomain.c datagenf.o g2ps.o \
  $(LIBS) -o g2hdomain
	./g2hdomain
	rm -f ./g2hdomain

degred1.ps: degred.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) degred.c $(LIBS) -o degred
	./degred
	rm -f ./degred

degred2.ps: degred.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) degred.c $(LIBS) -o degred
	./degred
	rm -f ./degred

addspl.ps: addspl.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) addspl.c $(LIBS) -o addspl
	./addspl
	rm -f ./addspl

bsplane.ps: bsplane.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bsplane.c $(LIBS) -o bsplane
	./bsplane
	rm -f ./bsplane

bsmesh.ps: bsmesh.c
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bsmesh.c $(LIBS) -o bsmesh
	./bsmesh
	rm -f ./bsmesh

clean:
	rm -f *.o
	rm -f $(PICTURES)
	rm -f $(PROGRAMS)

mrproper:
	rm -f *.o
	rm -f $(PICTURES)
	rm -f *.ps
	rm -f $(PROGRAMS)
