# Makefile for samtools, utilities for the Sequence Alignment/Map format.
#
#    Copyright (C) 2008-2021 Genome Research Ltd.
#    Portions copyright (C) 2010-2012 Broad Institute.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

CC       = gcc
CXX	= g++
AR       = ar
AWK      = awk
CFLAGS   = -g -Wall -O2 -pedantic -std=c99 -D_XOPEN_SOURCE=600
CPPFLAGS  = -std=c++11 -g -Wall -O3
LDFLAGS  = 
LIBS     =

OBJ = ParsingBam.o Util.o PhasingProcess.o Phasing.o PhasingGraph.o main.o

PROGRAMS = longphase

all: $(PROGRAMS) 

ALL_CPPFLAGS = -I. -Ihtslib $(CPPFLAGS)
ALL_LDFLAGS  = $(LDFLAGS)

# Usually config.mk and config.h are generated by running configure
# or config.status, but if those aren't used create defaults here.

config.mk:
	@sed -e '/^prefix/,/^LIBS/d;s/@Hsource@//;s/@Hinstall@/#/;s#@HTSDIR@#htslib#g;s/@HTSLIB_CPPFLAGS@/-I$$(HTSDIR)/g;s/@CURSES_LIB@/-lcurses/g' config.mk.in > $@

config.h:
	echo '/* Basic config.h generated by Makefile */' > $@
	echo '#define HAVE_CURSES' >> $@
	echo '#define HAVE_CURSES_H' >> $@

include config.mk


$(PROGRAMS): $(OBJ) $(HTSLIB)
	$(CXX) $(ALL_CPPFLAGS) -o $@ $^ $(HTSLIB_LIB)

%.o: %.cpp
	$(CXX) $(ALL_CPPFLAGS) -o $@ -c $^ 


mostlyclean: testclean
	-rm -f *.o 

clean: mostlyclean
	-rm -f $(PROGRAMS) 

distclean: clean
	-rm -f config.cache config.h config.log config.mk config.status
	-rm -f TAGS
	-rm -rf autom4te.cache

clean-all: clean clean-htslib

distclean-all: distclean distclean-htslib

mostlyclean-all: mostlyclean mostlyclean-htslib

testclean-all: testclean testclean-htslib

tags:
	ctags -f TAGS *.[ch] misc/*.[ch]


force:


.PHONY: all check check-all clean clean-all distclean distclean-all force
.PHONY: install lib mostlyclean mostlyclean-all print-version tags
.PHONY: test test-all testclean testclean-all
