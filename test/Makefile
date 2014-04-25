####### Select target computer
SYSTYPE = 'PUC-geryon'
#SYSTYPE = 'i686'
#SYSTYPE = 'i386'

####### Compilation options
#OPT = -DDEBUG

####### Files
EXEC  = rdtest
OBJS  = rdtest.o allvars_test.o read_rscale_table.o \
        ../lib/nrsrc/nrutil.o ../lib/nrsrc/locate.o \
        ../lib/nrsrc/dpolint.o ../lib/nrsrc/dsubmatrix.o \
        ../lib/nrsrc/d4tensor.o ../lib/nrsrc/polin2.o

INCL  = allvars_test.h proto_test.h ../lib/nrsrc/nr.h \
        ../lib/nrsrc/nrsag.h ../lib/nrsrc/nrutil.h

.KEEP_STATE:

####### Compiler settings for different systems
CC       = gcc       # Default values
OPTIMIZE = -O2

LIBS = -lm

CFLAGS = $(OPTIMIZE) $(OPT)

ifeq ($(SYSTYPE),'PUC-geryon')
  CFLAGS += -march=nocona -pipe -I/usr/local/include
endif

ifeq ($(SYSTYPE),'i686')
  CFLAGS += -march=native -Wno-unused-result
endif

ifeq ($(SYSTYPE),'i386')
  CFLAGS += -Wall -pedantic -std=c99
endif

rdtables: $(OBJS) Makefile
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LIBS)

clean:
	rm -f $(OBJS) $(EXEC) *~ ../lib/nrsrc/*~

tidy:
	rm -f $(OBJS) *~ ../lib/nrsrc/*~
