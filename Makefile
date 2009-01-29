# hamc = Hall A Monte Carlo
# R. Michaels, May 2008


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
   LIBS = $(GLIBS) $(ROOTLIBS) $(ROOTGLIBS) -lg2c

endif

MAKEDEPEND    = $(GCC)

ALL_LIBS = $(LIBS) 

SRCDIR=./src
INCLUDES += -I$(SRCDIR)

ifdef PROFILE
   CXXFLAGS += -pg
endif

ifdef OPTIMIZE
   CXXFLAGS += -O
else
   CXXFLAGS += -g
endif

SRC = $(SRCDIR)/hamcExpt.C $(SRCDIR)/hamcSingles.C $(SRCDIR)/hamcPhysics.C \
      $(SRCDIR)/hamcEvent.C $(SRCDIR)/hamcSpecHRS.C \
      $(SRCDIR)/hamcTrans.C $(SRCDIR)/hamcTransMat.C \
      $(SRCDIR)/hamcAccAvg.C \
      $(SRCDIR)/hamcTransLerWarmSeptum.C \
      $(SRCDIR)/hamcTransLerColdSeptum.C \
      $(SRCDIR)/hamcTransLerHRS.C \
      $(SRCDIR)/hamcTarget.C \
      $(SRCDIR)/hamcEloss.C $(SRCDIR)/hamcKine.C \
      $(SRCDIR)/hamcTrack.C $(SRCDIR)/hamcBeam.C $(SRCDIR)/hamcTrackOut.C \
      $(SRCDIR)/hamcInout.C $(SRCDIR)/THaString.C

DEPS = $(SRC:.C=.d)
DEP  = $(SRC:.C=.d)
HEAD = $(SRC:.C=.h) 

PROGS = prex happex
HAMCLIBS = libhamc.a
HAMCLIBS_NODICT = libhamc_NODICT.a

# Make the dictionary
ifdef MAKENODICTIONARY
  OBJS = $(SRC:.C=_NODICT.o)
else
  OBJS = $(SRC:.C=.o)
  OBJS += hamcDict.o
endif
OBJS += $(SRCDIR)/prex_forward.o $(SRCDIR)/monte_trans_hrs.o $(SRCDIR)/R6_forward.o $(SRCDIR)/ls_6d_forward.o

# PREX experiment
PREX_SRC = ./PREX/main_PREX.C ./PREX/hamcExptPREX.C ./PREX/hamcPhyPREX.C ./PREX/hamcTgtPREX.C
# Make the dictionary
ifdef MAKENODICTIONARY
  PREX_OBJS = $(PREX_SRC:.C=_NODICT.o)
else
  PREX_OBJS = $(PREX_SRC:.C=.o)
endif

# HAPPEX experiment
HAPPEX_SRC = ./HAPPEX/main_HAPPEX.C ./HAPPEX/hamcExptHAPPEX.C ./HAPPEX/hamcPhyHAPPEX.C ./HAPPEX/hamcTgtHAPPEX.C
# Make the dictionary
ifdef MAKENODICTIONARY
  HAPPEX_OBJS = $(HAPPEX_SRC:.C=_NODICT.o)
else
  HAPPEX_OBJS = $(HAPPEX_SRC:.C=.o)
endif

install: all

all: $(PROGS) $(HAMCLIBS) $(HAMCLIBS_NODICT) libhamc.so 

prex: $(PREX_OBJS) $(PREX_HEAD) $(OBJS) $(SRC) $(HEAD) 
	rm -f $@
	$(LD) $(CXXFLAGS) -o $@ $(OBJS) $(PREX_OBJS) $(ALL_LIBS)

happex: $(HAPPEX_OBJS) $(HAPPEX_HEAD) $(OBJS) $(SRC)  $(HEAD) 
	rm -f $@
	$(LD) $(CXXFLAGS) -o $@ $(OBJS) $(HAPPEX_OBJS) $(ALL_LIBS)

#Add PVDIS here when it exists


$(HAMCLIBS_NODICT): $(LOBJS_NODICT) $(LSRC) $(HEAD)
	rm -f $@
	ar cr $@ $(LOBJS)

$(HAMCLIBS): $(LOBJS) $(LSRC) $(HEAD)
	rm -f $@
	ar cr $@ $(LOBJS)

$(SRCDIR)/prex_forward.o: $(SRCDIR)/prex_forward.f
	rm -f $@
	cd $(SRCDIR) ; g77 -c prex_forward.f ; cd ../

$(SRCDIR)/monte_trans_hrs.o: $(SRCDIR)/monte_trans_hrs.f
	rm -f $@
	cd $(SRCDIR) ; g77 -c monte_trans_hrs.f ; cd ../

$(SRCDIR)/ls_6d_forward.o: $(SRCDIR)/ls_6d_forward.f
	rm -f $@
	cd $(SRCDIR) ; g77 -c ls_6d_forward.f ; cd ../

$(SRCDIR)/R6_forward.o: $(SRCDIR)/R6_forward.f
	rm -f $@
	cd $(SRCDIR) ; g77 -c R6_forward.f ; cd ../

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
	rm -f $(SRCDIR)/*.o core hamcDict* $(PROGS) $(HAMCLIBS) $(HAMCLIBS_NODICT) libhamc.so ./PREX/*.o ./HAPPEX/*.o 

realclean:  clean
	rm -f *.d *.tar  *~


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
