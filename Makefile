CC=gcc
LD=gcc
AR=ar

prefix := /usr/local

CFLAGS=-std=gnu99 -Wall -Werror -O2
ARFLAGS=rcs

SRCDIR=./src

# note order
TEST_LINKFLAGS=-L$(SRCDIR) -lpsfex -lcfitsio -lm

LIB_SOURCES = $(SRCDIR)/psfex.c $(SRCDIR)/psfex_fits.c

TEST_SOURCES = $(SRCDIR)/test.c

ALL_SOURCES = $(LIB_SOURCES) \
			  $(TEST_SOURCES)

LIB_OBJECTS=$(patsubst %.c,%.o,$(LIB_SOURCES)) 
TEST_OBJECTS=$(patsubst %.c,%.o,$(TEST_SOURCES)) 

# these installed
LIB_BASE=libpsfex.a
LIB = $(SRCDIR)/$(LIB_BASE)
HEADER = $(SRCDIR)/psfex.h

# just for tests
TEST_PROG = $(SRCDIR)/test

PROGS=$(TEST_PROG)
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

all: $(TEST_PROG)

lib: $(LIB)
	
$(LIB): $(LIB_OBJECTS)
	$(AR) $(ARFLAGS) $(LIB) $(LIB_OBJECTS)

$(TEST_PROG): $(LIB) $(TEST_OBJECTS)
	$(LD) -o $@  $(TEST_OBJECTS) $(TEST_LINKFLAGS)

clean:
	rm -f $(SRCDIR)/*.o $(LIB) $(PROGS) $(DEPFILE)
