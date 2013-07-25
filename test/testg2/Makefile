
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2007, 2008
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

HOME = /home/przemek
INCDIR = ../../include
LIBDIR = -L/usr/X11R6/lib -L../../lib

include ../../options.mak

OBJECTSF = datagenf.o testgraphf.o

OBJECTSD = datagend.o testgraphd.o

LIBS = -leghole -lmultibs -lcamera -lpsout -lpkgeom -lpknum -lpkvaria -lm

default: test0f test1f test2f test3f test4f test5f test6f \
         test7f test8f test9f \
         test0d test1d test2d test3d test4d test5d test6d \
         test7d test9d test9d

.c.o:
	$(CC) -c $(CFLAGS) -I$(INCDIR) -c $< -o $*.o

datagenf.o: datagenf.c datagenf.h

datagend.o: datagend.c datagend.h

bslapf.o: bslapf.c bslapf.h

bslapd.o: bslapd.c bslapd.h

drawitf.o: drawitf.c drawitf.h

drawitd.o: drawitd.c drawitd.h

testgraphf.o: testgraphf.c testgraphf.h

testgraphd.o: testgraphd.c testgraphd.h

readdatd.o: readdatd.c readdatd.h

writdatf.o: writdatf.c writdatf.h

writdatd.o: writdatd.c writdatd.h

test0f: test0f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test0f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test0f
	./test0f

test0d: test0d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test0d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test0d
	./test0d

test1f: test1f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test1f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test1f
	./test1f

test1d: test1d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test1d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test1d
	./test1d

test2f: test2f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test2f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test2f
	./test2f

test2d: test2d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test2d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test2d
	./test2d

test3f: test3f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test3f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test3f
	./test3f

test3d: test3d.c datagend.o testgraphd.o \
  readdatd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test3d.c $(OBJECTSD) readdatd.o \
  $(LIBDIR) $(LIBS) -o test3d
	./test3d

test4f: test4f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test4f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test4f
	./test4f

test4d: test4d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test4d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test4d
	./test4d

test5f: test5f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test5f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test5f
	./test5f

test5d: test5d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test5d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test5d
	./test5d

test6f: test6f.c datagenf.o testgraphf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test6f.c $(OBJECTSF) \
  $(LIBDIR) $(LIBS) -o test6f
	./test6f

test6d: test6d.c datagend.o testgraphd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test6d.c $(OBJECTSD) \
  $(LIBDIR) $(LIBS) -o test6d
	./test6d

test7f: test7f.c datagenf.o testgraphf.o drawitf.o bslapf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test7f.c $(OBJECTSF) drawitf.o bslapf.o \
  $(LIBDIR) $(LIBS) -o test7f
	./test7f

test7d: test7d.c datagend.o testgraphd.o drawitd.o bslapd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test7d.c $(OBJECTSD) drawitd.o bslapd.o \
  $(LIBDIR) $(LIBS) -o test7d
	./test7d

test8f: test8f.c datagenf.o testgraphf.o drawitf.o writdatf.o
	$(CC) $(CFLAGS) -I$(INCDIR) test8f.c $(OBJECTSF) \
  drawitf.o bslapf.o writdatf.o $(LIBDIR) $(LIBS) -o test8f
	./test8f

test8d: test8d.c datagend.o testgraphd.o drawitd.o writdatd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test8d.c $(OBJECTSD) \
  drawitd.o bslapd.o writdatd.o $(LIBDIR) $(LIBS) -o test8d
	./test8d

test9f: test9f.c datagenf.o testgraphf.o drawitf.o writdatf.o 
	$(CC) $(CFLAGS) -I$(INCDIR) test9f.c $(OBJECTSF) \
  drawitf.o bslapf.o writdatf.o $(LIBDIR) $(LIBS) -o test9f
	./test9f

test9d: test9d.c datagend.o testgraphd.o drawitd.o  writdatd.o
	$(CC) $(CFLAGS) -I$(INCDIR) test9d.c $(OBJECTSD) \
  drawitd.o bslapd.o writdatd.o $(LIBDIR) $(LIBS) -o test9d
	./test9d

clean:
	rm -f *.o
	rm -f *.ps
	rm -f test?f test?d

mrproper:
	rm -f *.o
	rm -f *.ps
	rm -f test?f test?d

