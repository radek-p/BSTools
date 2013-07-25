
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2013
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

#
#first enter the ../src directory and make
#

HOME = /home/przemek
SRCDIR = .
INCDIR = ../include
LIBDIR = ../lib
LIBSFN = $(LIBDIR)/libpkvaria.a $(LIBDIR)/libpkgeom.a $(LIBDIR)/libcamera.a \
  $(LIBDIR)/libmultibs.a $(LIBDIR)/libpsout.a
LIBS = -leghole -lraybez -lcamera -lmultibs -lpsout -lpkgeom -lpknum -lpkvaria \
  -lm

include ../options.mak

TARGET = nvtest knotrem bspdegel bscfe bscfede oslo4 bspbez degel \
  skiplr intbsc raster patchder patchpdir lsq transp bsder \
  addspl raytest trimpatch bsder123f choleskyf clipf rclipf \
  testg2f chknots degredf bscoonsf clcuintf clcuintd \
  nvtestd knotremd bspdegeld bscfed bscfeded oslo4d bspbezd degeld \
  skiplrd intbscd rasterd patchderd patchpdird lsqd transpd bsderd \
  addspld raytestd trimpatchd bsder123d choleskyd clipd rclipd \
  testg2d chknotsd degredd bscoonsd \
  qsort pksort intersf intersd bcoffsetf bcoffsetd \
  lmttestf lmttestd


default: $(TARGET)

.c.o:
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

nvtest.o: nvtest.c $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/multibs.h \
  $(INCDIR)/psout.h

nvtestd.o: nvtestd.c $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/multibs.h \
  $(INCDIR)/psout.h

knotrem.o: knotrem.c $(INCDIR)/pkvaria.h $(INCDIR)/pkgeom.h $(INCDIR)/multibs.h \
  $(INCDIR)/psout.h

knotremd.o: knotremd.c $(INCDIR)/pkvaria.h $(INCDIR)/pkgeom.h $(INCDIR)/multibs.h \
  $(INCDIR)/psout.h

bspdegel.o: bspdegel.c $(INCDIR)/pkvaria.h $(INCDIR)/pkgeom.h $(INCDIR)/camera.h \
  $(INCDIR)/multibs.h $(INCDIR)/psout.h

bspdegeld.o: bspdegeld.c $(INCDIR)/pkvaria.h $(INCDIR)/pkgeom.h $(INCDIR)/camera.h \
  $(INCDIR)/multibs.h $(INCDIR)/psout.h

chknots: chknots.c $(INCDIR)/pkvaria.h $(INCDIR)/pkgeom.h $(INCDIR)/camera.h \
  $(INCDIR)/multibs.h $(INCDIR)/psout.h
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) chknots.c $(LIBS) -o chknots
	./chknots

chknotsd: chknotsd.c $(INCDIR)/pkvaria.h $(INCDIR)/pkgeom.h $(INCDIR)/camera.h \
  $(INCDIR)/multibs.h $(INCDIR)/psout.h
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) chknotsd.c $(LIBS) -o chknotsd
	./chknotsd

raytest: raytest.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) raytest.c $(LIBS) -o raytest
	./raytest

raytestd: raytestd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) raytestd.c $(LIBS) -o raytestd
	./raytestd

nvtest: nvtest.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) nvtest.c $(LIBS) -o nvtest
	./nvtest

nvtestd: nvtestd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) nvtestd.c $(LIBS) -o nvtestd
	./nvtestd

knotrem: knotrem.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) knotrem.c $(LIBS) -o knotrem
	./knotrem

knotremd: knotremd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) knotremd.c $(LIBS) -o knotremd
	./knotremd

bspdegel: bspdegel.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bspdegel.c $(LIBS) -o bspdegel
	./bspdegel

bspdegeld: bspdegeld.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bspdegeld.c $(LIBS) -o bspdegeld
	./bspdegeld

degel: degel.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) degel.c $(LIBS) -o degel
	./degel

degeld: degeld.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) degeld.c $(LIBS) -o degeld
	./degeld

bscfe: bscfe.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bscfe.c $(LIBS) -o bscfe
	./bscfe

bscfed: bscfed.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bscfed.c $(LIBS) -o bscfed
	./bscfed

bscfede: bscfede.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bscfede.c $(LIBS) -o bscfede
	./bscfede

bscfeded: bscfeded.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bscfeded.c $(LIBS) -o bscfeded
	./bscfeded

oslo4: oslo4.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) oslo4.c $(LIBS) -o oslo4
	./oslo4

oslo4d: oslo4d.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) oslo4d.c $(LIBS) -o oslo4d
	./oslo4d

bspbez: bspbez.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bspbez.c $(LIBS) -o bspbez
	./bspbez

bspbezd: bspbezd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bspbezd.c $(LIBS) -o bspbezd
	./bspbezd

skiplr: skiplr.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) skiplr.c $(LIBS) -o skiplr
	./skiplr

skiplrd: skiplrd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) skiplrd.c $(LIBS) -o skiplrd
	./skiplrd

intbsc: intbsc.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) intbsc.c $(LIBS) -o intbsc
	./intbsc

intbscd: intbscd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) intbscd.c $(LIBS) -o intbscd
	./intbscd

raster: raster.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) raster.c $(LIBS) -o raster
	./raster

rasterd: rasterd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) rasterd.c $(LIBS) -o rasterd
	./rasterd

patchder: patchder.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) patchder.c $(LIBS) -o patchder
	./patchder

patchderd: patchderd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) patchderd.c $(LIBS) -o patchderd
	./patchderd

patchpdir: patchpdir.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) patchpdir.c $(LIBS) -o patchpdir
	./patchpdir

patchpdird: patchpdird.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) patchpdird.c $(LIBS) -o patchpdird
	./patchpdird

lsq: lsq.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) lsq.c $(LIBS) -o lsq
	./lsq

lsqd: lsqd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) lsqd.c $(LIBS) -o lsqd
	./lsqd

transp: transp.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) transp.c $(LIBS) -o transp
	./transp

transpd: transpd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) transpd.c $(LIBS) -o transpd
	./transpd

bsder: bsder.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bsder.c $(LIBS) -o bsder
	./bsder

bsderd: bsderd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bsderd.c $(LIBS) -o bsderd
	./bsderd

addspl: addspl.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) addspl.c $(LIBS) -o addspl
	./addspl

addspld: addspld.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) addspld.c $(LIBS) -o addspld
	./addspld

bsder123f: bsder123f.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bsder123f.c $(LIBS) -o bsder123f
	./bsder123f

bsder123d: bsder123d.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bsder123d.c $(LIBS) -o bsder123d
	./bsder123d

choleskyf: choleskyf.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) choleskyf.c $(LIBS) -o choleskyf
	./choleskyf

choleskyd: choleskyd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) choleskyd.c $(LIBS) -o choleskyd
	./choleskyd

clipf: clipf.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) clipf.c $(LIBS) -o clipf
	./clipf

clipd: clipd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) clipd.c $(LIBS) -o clipd
	./clipd

rclipf: rclipf.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) rclipf.c $(LIBS) -o rclipf
	./rclipf

rclipd: rclipd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) rclipd.c $(LIBS) -o rclipd
	./rclipd

clcuintf: clcuintf.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) clcuintf.c $(LIBS) -o clcuintf
	./clcuintf

clcuintd: clcuintd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) clcuintd.c $(LIBS) -o clcuintd
	./clcuintd

degredf: degredf.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) degredf.c $(LIBS) -o degredf
	./degredf

degredd: degredd.c $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) degredd.c $(LIBS) -o degredd
	./degredd

raytest.o: raytest.c $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/multibs.h \
  $(INCDIR)/raybez.h $(INCDIR)/psout.h

raytestd.o: raytestd.c $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/multibs.h \
  $(INCDIR)/raybez.h $(INCDIR)/psout.h

trimpatch.o: trimpatch.c $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/multibs.h \
  $(INCDIR)/psout.h

trimpatchd.o: trimpatchd.c $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/multibs.h \
  $(INCDIR)/psout.h

trimpatch: trimpatch.o $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) trimpatch.o $(LIBS) -o trimpatch
	./trimpatch

trimpatchd: trimpatchd.o $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) trimpatchd.o $(LIBS) -o trimpatchd
	./trimpatchd

datagenf.o: datagenf.c datagenf.h $(INCDIR)/pkvaria.h $(INCDIR)/pknum.h \
  $(INCDIR)/pkgeom.h $(INCDIR)/multibs.h $(INCDIR)/egholef.h

testg2f: testg2f.c datagenf.o $(INCDIR)/pkvaria.h $(INCDIR)/pknum.h \
  $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/psout.h \
  $(INCDIR)/multibs.h $(INCDIR)/egholef.h $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) testg2f.c datagenf.o $(LIBS) \
  -o testg2f
	./testg2f

datagend.o: datagend.c datagend.h $(INCDIR)/pkvaria.h $(INCDIR)/pknum.h \
  $(INCDIR)/pkgeom.h $(INCDIR)/multibs.h $(INCDIR)/egholed.h

testg2d: testg2d.c datagend.o $(INCDIR)/pkvaria.h $(INCDIR)/pknum.h \
  $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/psout.h \
  $(INCDIR)/multibs.h $(INCDIR)/egholed.h $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) testg2d.c datagend.o $(LIBS) \
  -o testg2d
	./testg2d

bscoonsf: bscoonsf.c $(INCDIR)/pkvaria.h $(INCDIR)/pknum.h \
  $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/psout.h \
  $(INCDIR)/multibs.h $(INCDIR)/egholed.h $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bscoonsf.c $(LIBS) \
  -o bscoonsf
	./bscoonsf

bscoonsd: bscoonsd.c $(INCDIR)/pkvaria.h $(INCDIR)/pknum.h \
  $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/psout.h \
  $(INCDIR)/multibs.h $(INCDIR)/egholed.h $(LIBSFN)
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bscoonsd.c $(LIBS) \
  -o bscoonsd
	./bscoonsd

bcoffsetf: bcoffsetf.c $(INCDIR)/pkvaria.h $(INCDIR)/pknum.h \
  $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/psout.h \
  $(INCDIR)/multibs.h
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bcoffsetf.c $(LIBS) \
  -o bcoffsetf
	./bcoffsetf

bcoffsetd: bcoffsetd.c $(INCDIR)/pkvaria.h $(INCDIR)/pknum.h \
  $(INCDIR)/pkgeom.h $(INCDIR)/camera.h $(INCDIR)/psout.h \
  $(INCDIR)/multibs.h
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) bcoffsetd.c $(LIBS) \
  -o bcoffsetd
	./bcoffsetd

qsort: qsort.c $(INCDIR)/pkvaria.h
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) qsort.c $(LIBS) -o qsort
	./qsort

pksort: pksort.c $(INCDIR)/pkvaria.h
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) pksort.c $(LIBS) -o pksort
	./pksort

intersf: intersf.c
	$(CC) $(CFLAGS) -I$(INCDIR) intersf.c -L$(LIBDIR) $(LIBS) -o intersf
	./intersf

intersd: intersd.c
	$(CC) $(CFLAGS) -I$(INCDIR) intersd.c -L$(LIBDIR) $(LIBS) -o intersd
	./intersd

lmttestf: lmttestf.c
	$(CC) $(CFLAGS) -I$(INCDIR) lmttestf.c -L$(LIBDIR) $(LIBS) -o lmttestf
	./lmttestf

lmttestd: lmttestd.c
	$(CC) $(CFLAGS) -I$(INCDIR) lmttestd.c -L$(LIBDIR) $(LIBS) -o lmttestd
	./lmttestd

alltests:
	(cd testg1;make)
	(cd testg2;make)

clean:
	rm -f $(TARGET)
	rm -f *.o
	(cd testg1;make clean)
	(cd testg2;make clean)

mrproper:
	rm -f $(TARGET)
	rm -f *.o
	rm -f *.ps
	(cd testg1;make mrproper)
	(cd testg2;make mrproper)

