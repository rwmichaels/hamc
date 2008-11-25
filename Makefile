# hamc = Hall A Monte Carlo
# R. Michaels, May 2008
# this is a test code

# Choose the compiler.
GCC=g++
GLD=g++

MAKENODICTIONARY=1

export OSNAME := $(shell uname)

ifeq ($(OSNAME),SunOS)

   ROOTCFLAGS    = $(shell root-config --cflags)
   ROOTLIBS      = $(shell root-config --libs)
   ROOTGLIBS     = $(shell root-config --glibs)
   CXX           = $(GCC)
   CXXFLAGS      = -KPIC -DSUNVERS -I$(ROOTSYS)/include -I$(MAINDIR)
   CXXFLAGS     += $(ROOTCFLAGS)
   LD            = $(GLD)
   LDFLAGS       = -g -D
   SOFLAGS       = -G
   GLIB =  -lm -lc -lgen -lw -lnsl -ldl
   SLIB = -L/opt/SUNWspro/SC4.2/lib -lF77 -lM77 -lsunmath
   ET_AC_FLAGS = -D_REENTRANT -D_POSIX_THREAD_SEMANTICS
   ET_CFLAGS = -mt -fast -xO5 -KPIC $(ET_AC_FLAGS) -DSUNVERS
   ONLIBS = -lposix4 -lnsl -lsocket -lresolv
   LIBS = $(GLIB)

endif

# Linux with egcs

ifeq ($(OSNAME),Linux)
   ROOTLIBS      = $(shell root-config --libs)
   ROOTGLIBS     = $(shell root-config --glibs)
   INCLUDES      = -I$(ROOTSYS)/include
   CXX           = $(GCC)
   CXXFLAGS      = -Wall -fno-exceptions -fPIC $(INCLUDES) 
   LD            = $(GLD)
   LDFLAGS       = 
   SOFLAGS       = -shared 
   GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11
   ET_AC_FLAGS = -D_REENTRANT -D_POSIX_PTHREAD_SEMANTICS
   ET_CFLAGS = -02 -fPIC -I. $(ET_AC_FLAGS) -DLINUXVERS
   ONLIBS = -lieee -lpthread -ldl -lresolv
   LIBS = $(GLIBS) $(ROOTLIBS) $(ROOTGLIBS)

endif

MAKEDEPEND    = $(GCC)

ALL_LIBS = $(LIBS) 

DCDIR=../codaclass
INCLUDES += -I$(DCDIR)

ifdef PROFILE
   CXXFLAGS += -pg
endif

ifdef OPTIMIZE
   CXXFLAGS += -O
else
   CXXFLAGS += -g
endif

SRC = hamcExpt.C hamcSingles.C hamcExptPREX.C \
      hamcPhysics.C hamcPhyPREX.C \
      hamcEvent.C hamcSpecHRS.C \
      hamcTrans.C hamcTransMat.C \
      hamcTransLerWarmSeptum.C \
      hamcTransLerColdSeptum.C \
      hamcTransLerHRS.C \
      hamcTarget.C hamcTgtPREX.C \
      hamcRad.C hamcKine.C \
      hamcTrack.C hamcBeam.C hamcTrackOut.C \
      hamcInout.C THaString.C

DEPS = $(SRC:.C=.d)
DEP  = $(SRC:.C=.d)
HEAD = $(SRC:.C=.h) 

PROGS = hamc
HAMCLIBS = libhamc.a
HAMCLIBS_NODICT = libhamc_NODICT.a

# Make the dictionary
ifdef MAKENODICTIONARY
  OBJS = $(SRC:.C=_NODICT.o)
else
  OBJS = $(SRC:.C=.o)
  OBJS += hamcDict.o
endif
OBJS += prex_forward.o monte_trans_hrs.o R6_forward.o ls_6d_forward.o

install: all

all: $(PROGS) $(HAMCLIBS) $(HAMCLIBS_NODICT) libhamc.so

$(PROGS): main.o $(OBJS) $(SRC) $(HEAD) 
	rm -f $@
	$(LD) $(CXXFLAGS) -o $@ main.o $(OBJS) $(ALL_LIBS)

$(HAMCLIBS_NODICT): $(LOBJS_NODICT) $(LSRC) $(HEAD)
	rm -f $@
	ar cr $@ $(LOBJS)

$(HAMCLIBS): $(LOBJS) $(LSRC) $(HEAD)
	rm -f $@
	ar cr $@ $(LOBJS)

prex_forward.o: prex_forward.f
	rm -f $@
	g77 -c prex_forward.f

main.o: main.C
	$(CXX) -c $(INCLUDES) $(CXXFLAGS) $<	

libhamc.so: $(OBJS) $(HEAD)
	$(CXX) $(SOFLAGS) -O -o libhamc.so $(OBJS) $(ALL_LIBS)

# Dictionary

hamcDict.C: $(OBJS) hamcLinkDef.h
	@echo "Generating Decoder Dictionary..."
	$(ROOTSYS)/bin/rootcint -f hamcDict.C -c -p $(HEAD) hamcLinkDef.h

# create a tar file of ./$(VERS)/* (all code)

tarfile: clean version
	tar cvf $(VERS).tar ./$(VERS)

clean:
	rm -f *.o core hamcDict* $(PROGS) $(HAMCLIBS) $(HAMCLIBS_NODICT) libhamc.so

realclean:  clean
	rm -f *.d *.tar *~


%.o:	%.C
	$(CXX) $(CXXFLAGS) -c $<

%_NODICT.o:	%.C
	$(CXX) $(CXXFLAGS) -c -DNODICT -o $*_NODICT.o $<

%.d:	%.C 
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(MAKEDEPEND) -MM $(INCLUDES) -c $< \
                | sed '\''s%^.*\.o%$*\.o%g'\'' \
                | sed '\''s%\($*\)\.o[ :]*%\1.o $@ : %g'\'' > $@; \
                [ -s $@ ] || rm -f $@'

-include $(DEPS)
