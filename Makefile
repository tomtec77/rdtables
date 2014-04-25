####### Select target computer
SYSTYPE = 'PUC-geryon'
#SYSTYPE = 'i686'
#SYSTYPE = 'i386'

####### Compilation options
OPT += -DTESTMODE    # To run a single calculation in test mode
#OPT += -DDEBUG

####### Files
EXEC  = rdtables
OBJS  = rdtables.o allvars.o readparameters.o init.o output_times.o \
        concentration.o axes.o generate_rscale_table.o rditer2.o\
        lib/myutils/myutils.o lib/myutils/validate.o lib/myutils/indexrm.o \
        lib/myutils/gethostname.o \
        lib/nrsrc/nrutil.o lib/nrsrc/locate.o lib/nrsrc/dpolint.o \
        lib/nrsrc/gaulag.o lib/nrsrc/gammln.o lib/nrsrc/bessel.o

INCL  = allvars.h proto.h readparameters.h lib/myutils/myutils.h

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
  CFLAGS += -Wall -pedantic -std=c99 -D_USE_BSD
endif

all: rdtables

rdtables: $(OBJS) Makefile
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) $(LIBS)

clean:
	rm -f $(OBJS) $(EXEC) *~ lib/myutils/*~ lib/nrsrc/*~

tidy:
	rm -f $(OBJS) *~ lib/myutils/*~ lib/nrsrc/*~
