
# /////////////////////////////////////////////////////////////////////
# This file is a part of the BSTools package
# written by Przemyslaw Kiciak
# /////////////////////////////////////////////////////////////////////
# (C) Copyright by Przemyslaw Kiciak, 2005, 2012
# this package is distributed under the terms of the
# Lesser GNU Public License, see the file COPYING.LIB
# /////////////////////////////////////////////////////////////////////

default: options.mak
	(cd src;make)
	(cd test;make)
	(cd doc;make)
	(cd demo;make)

# all tests include the tests of libg1hole and libg2hole,
# somewhat time consuming and generating PostScript pictures
# in rather large files
alltests:
	(make)
	(cd test;make alltests)

options.mak:
	rm -f options.mak
	ln -s gcc.mak options.mak

gcc:
	rm -f options.mak
	ln -s gcc.mak options.mak
	(cd src;make)
	(cd test;make)
	(cd doc;make)
	(cd demo;make)

32:
	rm -f options.mak
	ln -s gcc.32.mak options.mak
	(cd src;make)
	(cd test;make)
	(cd doc;make)
	(cd demo;make)

64:
	rm -f options.mak
	ln -s gcc.64.mak options.mak
	(cd src;make)
	(cd test;make)
	(cd doc;make)
	(cd demo;make)

native:
	rm -f options.mak
	ln -s gcc.native.mak options.mak
	(cd src;make)
	(cd test;make)
	(cd doc;make)
	(cd demo;make)

debug:
	rm -f options.mak
	ln -s gcc.debug.mak options.mak
	(cd src;make)
	(cd test;make)
	(cd doc;make)
	(cd demo;make)

clean: options.mak
	(cd src;make clean)
	(cd test;make clean)
	(cd doc;make clean)
	(cd demo;make clean)

mrproper: options.mak
	(cd src;make mrproper)
	(cd test;make mrproper)
	(cd doc;make mrproper)
	(cd demo;make mrproper)
	rm -f options.mak

