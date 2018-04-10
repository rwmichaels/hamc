# hamc = Hall A Monte Carlo
# R. Michaels, May 2008
# Update  Jan 2018  with USEFORTRAN option

# Choose the compiler.
GCC=g++
GLD=g++

MAKENODICTIONARY=1

# Using fortran-dependent code (=1) or not (=0), soon to be obsolete.
# If not, then rely on HRSTRans (libhrstrans) for transport.
USEFORTRAN=1

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
   LIBS = $(GLIB) /usr/lib64/libg2c.so.0

endif

# Linux with egcs

ifeq ($(OSNAME),Linux)
   ROOTLIBS      = $(shell root-config --libs)
   ROOTGLIBS     = $(shell root-config --glibs)
   INCLUDES      = -I$(ROOTSYS)/include
   CXX           = $(GCC)
   CXXFLAGS      = -fno-exceptions -fpermissive  -std=c++11 -fPIC $(INCLUDES)
   LD            = $(GLD)
   LDFLAGS       = 
   SOFLAGS       = -shared 
   GLIBS         = $(ROOTGLIBS) -L/usr/X11R6/lib -lXpm -lX11 
   ET_AC_FLAGS = -D_REENTRANT -D_POSIX_PTHREAD_SEMANTICS
   ET_CFLAGS = -02 -fPIC -I. $(ET_AC_FLAGS) -DLINUXVERS
   ONLIBS = -lieee -lpthread -ldl -lresolv 
   LIBS = $(GLIBS) $(ROOTLIBS) $(ROOTGLIBS) /usr/lib64/libg2c.so.0

endif

HRSTRANSLIB   = /home/rom/hrstrans1/lib

MAKEDEPEND    = $(GCC)

ALL_LIBS = $(LIBS) $(HRSTRANSLIB)/libhrstrans.so

SRCDIR=./src
INCLUDES += -I$(SRCDIR)
INCLUDES += -I$(HRSTRANSLIB)

ifdef PROFILE
   CXXFLAGS += -pg
endif

ifdef OPTIMIZE
   CXXFLAGS += -O
else
   CXXFLAGS += -g -ggdb
endif

SRC = $(SRCDIR)/hamcExpt.C $(SRCDIR)/hamcSingles.C $(SRCDIR)/hamcPhysics.C \
      $(SRCDIR)/hamcEvent.C $(SRCDIR)/hamcSpecHRS.C \
      $(SRCDIR)/hamcTrans.C $(SRCDIR)/hamcTransMat.C \
      $(SRCDIR)/hamcAccAvg.C \
      $(SRCDIR)/hamcTransGuido.C \
      $(SRCDIR)/hamcTHRSTrans.C \
      $(SRCDIR)/hamcTarget.C \
      $(SRCDIR)/hamcEloss.C $(SRCDIR)/hamcKine.C \
      $(SRCDIR)/hamcTrack.C $(SRCDIR)/hamcBeam.C $(SRCDIR)/hamcTrackOut.C \
      $(SRCDIR)/hamcInout.C $(SRCDIR)/THaString.C $(SRCDIR)/hamcMultScatt.C

# Fortran-dependent sources
# which of these you need depends on what you are compiling.
# too tired to think much now, but
#  Ler4deg was for CREX
#  WarmSeptum for CREX-I,  ColdSeptum is largely irrelevant.
#  LerHRS is Std. HRS, so HAPPEX-III etc.
#FORSRC=$(SRCDIR)/hamcTransLer4deg.C \
#      $(SRCDIR)/hamcTransLerWarmSeptum.C \
#      $(SRCDIR)/hamcTransLerColdSeptum.C \
#      $(SRCDIR)/hamcTransLerHRS.C 

FORSRC=$(SRCDIR)/hamcTransLerHRS.C $(SRCDIR)/hamcTransLerWarmSeptum.C $(SRCDIR)/hamcTransLerColdSeptum.C $(SRCDIR)/hamcTransLer4deg.C

ifdef USEFORTRAN 
   SRC+=$(FORSRC)
   CXXFLAGS += -DUSEFORTRAN -lgfortran
endif

DEPS = $(SRC:.C=.d)
DEP  = $(SRC:.C=.d)
HEAD = $(SRC:.C=.h) 

PROGS = prex happex 
ifdef DOPVDIS
  PROGS += pvdis
endif
HAMCLIBS = libhamc.a
HAMCLIBS_NODICT = libhamc_NODICT.a

# Make the dictionary
ifdef MAKENODICTIONARY
  OBJS = $(SRC:.C=_NODICT.o)
else
  OBJS = $(SRC:.C=.o)
  OBJS += hamcDict.o
endif

ifdef USEFORTRAN
# PREX-I or C-REX Optics HRS+septum
#OBJS += $(SRCDIR)/crex_4degr.o $(SRCDIR)/prex_forward.o $(SRCDIR)/monte_trans_hrs.o $(SRCDIR)/R6_forward.o $(SRCDIR)/ls_6d_forward.o
# Standard HRS Optics, i.e. HRS w/out Septum
#OBJS += $(SRCDIR)/prex_retune_for.o $(SRCDIR)/monte_trans_hrs.o $(SRCDIR)/R6_forward.o $(SRCDIR)/ls_6d_forward.o
# Try a combination
OBJS += $(SRCDIR)/crex_4degr.o $(SRCDIR)/prex_retune_for.o $(SRCDIR)/monte_trans_hrs.o $(SRCDIR)/R6_forward.o $(SRCDIR)/ls_6d_forward.o
endif

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

# PVDIS experiment
PVDIS_SRC = ./PVDIS/main_PVDIS.C ./PVDIS/hamcExptPVDIS.C ./PVDIS/hamcPhyPVDIS.C ./PVDIS/hamcTgtPVDIS.C 
ifdef MAKENODICTIONARY
  PVDIS_OBJS = $(PVDIS_SRC:.C=_NODICT.o)
else
  PVDIS_OBJS = $(PVDIS_SRC:.C=.o)
endif

ifdef DOPVDIS
  PVDIS_OBJS += ./PVDIS/getpdf_mrst2003c.o ./PVDIS/mrst2003c.o ./PVDIS/NextUn.o ./PVDIS/r1998.o ./PVDIS/readpdf_single.o ./PVDIS/xsec.o
endif

install: all

all: $(PROGS) $(HAMCLIBS) $(HAMCLIBS_NODICT) libhamc.so 

prex: $(PREX_OBJS) $(PREX_HEAD) $(OBJS) $(SRC) $(HEAD) 
	rm -f $@
	$(LD) $(CXXFLAGS) -o $@ $(OBJS) $(PREX_OBJS) $(ALL_LIBS)

happex: $(HAPPEX_OBJS) $(HAPPEX_HEAD) $(OBJS) $(SRC)  $(HEAD) 
	rm -f $@
	$(LD) $(CXXFLAGS) -o $@ $(OBJS) $(HAPPEX_OBJS) $(ALL_LIBS)

# Note, PVDIS relies heavily on Fortran.  This needs to be fixed if
# we care about it (noted, Jan 2018).

ifdef DOPVDIS
pvdis: $(PVDIS_OBJS) $(PVDIS_HEAD) $(OBJS) $(SRC) $(HEAD)
	rm -f $@
	$(LD) $(CXXFLAGS) -o $@ $(OBJS) $(PVDIS_OBJS) $(ALL_LIBS)
endif

$(HAMCLIBS_NODICT): $(LOBJS_NODICT) $(LSRC) $(HEAD)
	rm -f $@
	ar cr $@ $(LOBJS)

$(HAMCLIBS): $(LOBJS) $(LSRC) $(HEAD)
	rm -f $@
	ar cr $@ $(LOBJS)

# These are FORTRAN-dependent codes.
# They are obsolescent (soon to be obsolete), to be replaced
# by hrstranslib.
ifdef USEFORTRAN
# 4 degree Septum (April 2013)
$(SRCDIR)/crex_4degr.o: $(SRCDIR)/crex_4degr.f
	rm -f $@
	cd $(SRCDIR) ; g77 -c -fPIC crex_4degr.f ; cd ../

# Tune B (what was used during production)
$(SRCDIR)/prex_forward.o: $(SRCDIR)/prex_forward.f
	rm -f $@
	cd $(SRCDIR) ; g77 -c -fPIC prex_forward.f ; cd ../
#Replace 3 lines below with 3 lines above for standard HRS optics
# Tune Y:
$(SRCDIR)/prex_retune_for.o: $(SRCDIR)/prex_retune_for.f
	rm -f $@
	cd $(SRCDIR) ; g77 -c -fPIC prex_retune_for.f ; cd ../

$(SRCDIR)/monte_trans_hrs.o: $(SRCDIR)/monte_trans_hrs.f
	rm -f $@
	cd $(SRCDIR) ; g77 -c -fPIC monte_trans_hrs.f ; cd ../

$(SRCDIR)/ls_6d_forward.o: $(SRCDIR)/ls_6d_forward.f
	rm -f $@
	cd $(SRCDIR) ; g77 -c -fPIC ls_6d_forward.f ; cd ../

$(SRCDIR)/R6_forward.o: $(SRCDIR)/R6_forward.f
	rm -f $@
	cd $(SRCDIR) ; g77 -c -fPIC R6_forward.f ; cd ../
endif

#############################   PVDIS OBJS   ###############

# The PVDIS code wont work until we fix or replace fortran, I believe.
ifdef USEFORTRAN
./PVDIS/getpdf_mrst2003c.o: ./PVDIS/getpdf_mrst2003c.f ./PVDIS/mrst2003c.o ./PVDIS/NextUn.o ./PVDIS/r1998.o ./PVDIS/readpdf_single.o
	rm -f $@
	cd ./PVDIS/; g77 -c -fPIC getpdf_mrst2003c.f mrst2003c.o NextUn.o r1998.o readpdf_single.o; cd ../

./PVDIS/mrst2003c.o: ./PVDIS/mrst2003c.f
	rm -f $@
	cd ./PVDIS/; g77 -c -fPIC mrst2003c.f; cd ../

./PVDIS/NextUn.o: ./PVDIS/NextUn.f
	rm -f $@
	cd ./PVDIS/; g77 -c -fPIC NextUn.f; cd ../

./PVDIS/r1998.o: ./PVDIS/r1998.f
	rm -f $@
	cd ./PVDIS/; g77 -c -fPIC r1998.f; cd ../

./PVDIS/readpdf_single.o: ./PVDIS/readpdf_single.f
	rm -f $@
	cd ./PVDIS/; g77 -c -fPIC readpdf_single.f; cd ../

./PVDIS/xsec.o: ./PVDIS/xsec.f
	rm -f $@
	cd ./PVDIS/; g77 -c -fPIC xsec.f; cd ../
endif

################################################################

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
	rm -f $(SRCDIR)/*.o core hamcDict* $(PROGS) $(HAMCLIBS) $(HAMCLIBS_NODICT) libhamc.so ./PREX/*.o ./HAPPEX/*.o ./PVDIS/*.o

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
