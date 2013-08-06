CC=gcc
LD=gcc
AR=ar

prefix := /usr/local

CFLAGS=-std=gnu99 -Wall -O2 -I$(FITSIO_BASE)/include
ARFLAGS=rcs

SRCDIR=./src

# note order
TEST_LINKFLAGS=-L$(SRCDIR) -L$(FITSIO_BASE)/lib -lpsfex -lcfitsio -lm
REC_LINKFLAGS=-L$(SRCDIR) -L$(FITSIO_BASE)/lib -lpsfex -lcfitsio -lm

LIB_SOURCES = $(SRCDIR)/psfex.c $(SRCDIR)/psfex_fits.c $(SRCDIR)/poly.c

TEST_SOURCES = $(SRCDIR)/test.c
REC_SOURCES = $(SRCDIR)/psfex-rec.c

ALL_SOURCES = $(LIB_SOURCES) \
			  $(TEST_SOURCES) \
			  $(REC_SOURCES)

LIB_OBJECTS=$(patsubst %.c,%.o,$(LIB_SOURCES)) 
TEST_OBJECTS=$(patsubst %.c,%.o,$(TEST_SOURCES)) 
REC_OBJECTS=$(patsubst %.c,%.o,$(REC_SOURCES)) 

# these installed
LIB_BASE=libpsfex.a
LIB = $(SRCDIR)/$(LIB_BASE)
HEADER = $(SRCDIR)/psfex.h

TEST_PROG = $(SRCDIR)/test

REC_BASE = psfex-rec
REC_PROG = $(SRCDIR)/$(REC_BASE)

PROGS=$(TEST_PROG) $(REC_PROG)

DEPFILE=$(SRCDIR)/.depend

default: all

depend: $(DEPFILE)

$(DEPFILE): $(ALL_SOURCES)
	$(CC) $(CFLAGS) -MM $^ > $(DEPFILE);

-include $(DEPFILE)


install: $(LIB)
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/include
	mkdir -p $(prefix)/bin
	cp $(LIB) $(prefix)/lib/
	cp $(HEADER) $(prefix)/include/
	cp $(REC_PROG) $(prefix)/bin
	chmod a+x $(prefix)/bin/$(REC_BASE)

all: $(TEST_PROG) $(REC_PROG)

lib: $(LIB)

$(LIB): $(LIB_OBJECTS)
	$(AR) $(ARFLAGS) $(LIB) $(LIB_OBJECTS)

$(TEST_PROG): $(LIB) $(TEST_OBJECTS)
	$(LD) -o $@  $(TEST_OBJECTS) $(TEST_LINKFLAGS)

$(REC_PROG): $(LIB) $(REC_OBJECTS)
	$(LD) -o $@  $(REC_OBJECTS) $(REC_LINKFLAGS)

clean:
	rm -f $(SRCDIR)/*.o $(LIB) $(PROGS) $(DEPFILE)
