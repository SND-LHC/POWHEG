#-*- Makefile -*-
## Choose compiler: gfortran,ifort (g77 not supported, F90 constructs in use!)
COMPILER=gfortran
FC=$(COMPILER)
#CXX=g++
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf
#Choose Analysis: dummy, process specific
## default analysis may require FASTJET package, that has to be installed separately (see below)
# trijet     : all-round analysis
# minlo      : an analysis to test the minlo option with a pt cut on the hardest jet and 
# various pt cuts on the second hardest
# minlo-ycut : the same as the previous, but this with an additional rapidity cut on all jets 
# lhc        : Analyses originally used by various experiments at LHC
# tev        : Analysis originally used by various experiments at TeVatron
#ANALYSIS=lh-sh
ANALYSIS=paper-Rdep
#ANALYSIS=paper
## For static linking uncomment the following
#STATIC= -static
#

OBJ=obj-$(COMPILER)

ifeq ("$(COMPILER)","gfortran")	
F77=gfortran -fno-automatic
CXX=g++
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
#FFLAGS= -Wall -Wimplicit-interface -fbounds-check
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O2
## For debugging uncomment the following
#DEBUG= -ggdb 
ifdef DEBUG
OPT=-O0
FPE=-ffpe-trap=invalid,zero,overflow
endif
endif

ifeq ("$(COMPILER)","ifort")
F77 = ifort -save  -extend_source
CXX=icpc
#CXX = g++
#LIBS = -limf
#FFLAGS =  -checkm
## For floating point exception trapping  uncomment the following 
#FPE = -fpe0
OPT = -O1 #-fast -check -ftrapuv 
## For debugging uncomment the following
#DEBUG=  -g
ifdef DEBUG
OPT=-O0 
FPE = -fpe0
endif
endif

CXXFLAGS=$(OPT)

PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))
VPATH= ./:../:$(OBJ)

INCLUDE0=$(PWD)
INCLUDE1=$(PWD)/include
#INCLUDE1.2=$(shell dirname $(PWD))/progress/include
INCLUDE2=$(shell dirname $(PWD))/include
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1)  -I$(INCLUDE2)

INCLUDE =$(wildcard ../include/*.h *.h include/*.h)

ifeq ("$(PDF)","lhapdf")
LHAPDF_CONFIG=lhapdf-config
PDFPACK=lhapdf6if.o lhapdf6ifcc.o
LHCXXFLAGS+= $(shell $(LHAPDF_CONFIG) --cxxflags)
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF -lstdc++
ifeq  ("$(STATIC)","-static") 
## If LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
## libgfortran as well. The same holds for libstdc++. 
## One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
## does perform this inclusion. The path has to be set by the user. 
# LIBGFORTRANPATH= #/usr/lib/gcc/x86_64-redhat-linux/4.1.2
# LIBSTDCPP=/lib64
LIBSLHAPDF+=  -L$(LIBGFORTRANPATH)  -lgfortranbegin -lgfortran -L$(LIBSTDCPP) -lstdc++
endif
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif

# PYTHIA 8
PYTHIA8LOCATION=$(shell pythia8-config --prefix)
PY8CXXFLAGS+=-I$(PYTHIA8LOCATION)/include/ -I$(PYTHIA8LOCATION)/include/Pythia8
LIBPYTHIA8=-L$(PYTHIA8LOCATION)/lib/archive -lpythia8 -llhapdfdummy -lstdc++

PWHGANAL=pwhg_bookhist-new.o pwhg_analysis-dummy.o
ifeq ("$(ANALYSIS)","trijet")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o  pwhg_analysis.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetfortran.o
PWHGANAL+=  fastjetfortran.o fastjetfortran-BBUY.o
endif
ifeq ("$(ANALYSIS)","minlo")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o  pwhg_analysis-minlo.o  observables.o 
PWHGANAL+=  fastjetfortran.o
endif
ifeq ("$(ANALYSIS)","minlo-ycut")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o observables.o pwhg_analysis-minlo-ycut.o
PWHGANAL+=  fastjetfortran.o
endif
ifeq ("$(ANALYSIS)","paper")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o observables.o pwhg_analysis-minlo-ycut-paper.o
PWHGANAL+=  fastjetfortran.o
endif
ifeq ("$(ANALYSIS)","paper-Rdep")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o observables.o pwhg_analysis-minlo-ycut-paper-Rdep.o
PWHGANAL+=  fastjetfortran.o
endif
ifeq ("$(ANALYSIS)","lhc")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o  pwhg_analysis-lhc.o
PWHGANAL+= fastjetfortran.o fastjet-addon.o  splitandfilt.o \
	   observables.o
endif
ifeq ("$(ANALYSIS)","tev")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o  pwhg_analysis-tev.o
PWHGANAL+=  splitandfilt.o observables.o fastjetD0RunIIConewrap.o
endif


ifeq ("$(ANALYSIS)","lh-sh")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o observables.o pwhg_analysis-lh-sh.o
PWHGANAL+=  fastjetfortran.o
endif


%.o: %.f $(INCLUDE) | $(OBJ)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.f90 $(INCLUDE) | $(OBJ)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.c | $(OBJ)
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^ 

%.o: %.cc | $(OBJ)
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS) $(CXXFLAGS)


LIBS+=-lz

USER= init_processes_born_real.o init_processes.o  init_couplings.o \
     Born_phsp.o Born.o virtual.o real.o  $(PWHGANAL)


PWHG= $(USER) pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_analysis_driver.o checkmomzero.o		                \
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	reshufflemoms.o \
	sigcollsoft.o sigvirtual.o  ubprojections.o                     \
	pwhgreweight.o setlocalscales.o mint_upb.o opencount.o          \
	merging.o validflav.o $(PDFPACK) $(FPEOBJ) lhefread.o pwhg_io_interface.o rwl_weightlists.o rwl_setup_param_weights.o

LIBDIRMG=$(OBJ)
LINKMGLIBS =  -L$(LIBDIRMG)  -lmadgraph -lmodel -ldhelas3 
MADLIBS=libdhelas3.a libmadgraph.a libmodel.a

LINKNLOJETLIB = -lnlojet
NLOJETLIB=libnlojet.a

# target to generate LHEF output

pwhg_main: $(PWHG) $(MADLIBS) $(NLOJETLIB)
	$(FF) $(patsubst %,$(OBJ)/%,$(PWHG)) $(LINKMGLIBS) $(LINKNLOJETLIB) $(LIBS) $(LIBSFASTJET) $(STATIC)  -lstdc++ -o  $@


LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o		\
     opencount.o powheginput.o $(PWHGANAL)	\
     lhefread.o pwhg_io_interface.o rwl_weightlists.o newunit.o pwhg_analysis_driver.o $(FPEOBJ)

# target to analyze LHEF output
lhef_analysis:$(LHEF)
	$(FF) $(patsubst %,$(OBJ)/%,$(LHEF)) $(LIBS) $(LIBSFASTJET) $(STATIC)  -o $@ 



# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig.o boostrot.o bra_ket_subroutines.o	\
	powheginput.o $(PWHGANAL) lhefread.o pwhg_io_interface.o rwl_weightlists.o newunit.o pdfdummies.o changescalup.o	\
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o $(FPEOBJ) 

main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,$(OBJ)/%,$(HERWIG))  $(LIBSFASTJET)  $(STATIC) -o $@

# target to read event file, shower events with HERWIG + analysis
HERWIGNOHAD=main-HERWIG-nohad.o setup-HERWIG-lhef.o herwig.o boostrot.o bra_ket_subroutines.o	\
	powheginput.o $(PWHGANAL) lhefread.o pwhg_io_interface.o rwl_weightlists.o newunit.o pdfdummies.o changescalup.o	\
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o $(FPEOBJ) 

main-HERWIGNOHAD-lhef: $(HERWIGNOHAD)
	$(FF) $(patsubst %,$(OBJ)/%,$(HERWIGNOHAD))  $(LIBSFASTJET)  $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o \
	$(PWHGANAL) lhefread.o pwhg_io_interface.o rwl_weightlists.o newunit.o bra_ket_subroutines.o	\
        reshufflemoms.o changescalup.o \
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	$(FPEOBJ)

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA)) $(LIBS) $(LIBSFASTJET)  $(STATIC) -o $@

# target to read event file, shower events with PYTHIA8 + analysis
PYTHIA8=main-PYTHIA8.o pythia8F77.o boostrot.o powheginput.o \
	$(PWHGANAL) opencount.o lhefread.o pwhg_io_interface.o rwl_weightlists.o newunit.o pdfdummies.o \
	random.o cernroutines.o bra_ket_subroutines.o changescalup.o \
	$(FPEOBJ)

main-PYTHIA8-lhef: $(PYTHIA8)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA8)) $(LIBSFASTJET) $(LIBPYTHIA8) $(LIBS) $(STATIC) -o $@

# target to cleanup
.PHONY: clean libdhelas3.a libmadgraph.a libmodel.a libnlojet.a

XFFLAGS=$(DEBUG) $(OPT)

libdhelas3.a:
	cd DHELAS ; make -j FC="$(F77)" F77="$(F77)" XFFLAGS="$(XFFLAGS)" OBJ="$(OBJ)"; wait

libmadgraph.a:
	cd Madlib ; make -j FC="$(F77)" F77="$(F77)" XFFLAGS"=$(XFFLAGS)" OBJ="$(OBJ)"; wait

libnlojet.a:
	cd nlojet ; make -j CC="$(CXX)" CDEBUG="$(DEBUG)"  CXXFLAGS="$(CXXFLAGS)" OBJ="$(OBJ)"; wait


ifeq ("$(COMPILER)","gfortran")
XFFLAGS +=-ffixed-line-length-132
else
XFFLAGS +=-extend-source
endif

libmodel.a:
	cd MODEL ; make FC="$(F77)" F77="$(F77)" XFFLAGS="$(XFFLAGS)" OBJ="$(OBJ)"; wait


clean:
	rm -f $(patsubst %,$(OBJ)/%,$(USER) $(PWHG) $(LHEF) $(HERWIG) $(PYTHIA)) pwhg_main lhef_analysis main-HERWIG-lhef	\
	main-PYTHIA-lhef main-PYTHIA8-lhef main-HERWIGNOHAD-lhef  $(OBJ)/libnlojet.a


veryclean:
	rm -f $(OBJ)/*.o $(OBJ)/*.mod $(OBJ)/*.a pwhg_main lhef_analysis main-HERWIG-lhef	\
	main-PYTHIA-lhef main-PYTHIA8-lhef main-HERWIGNOHAD-lhef

# target to generate object directory if it does not exist
$(OBJ):
	mkdir -p $(OBJ)

##########################################################################


## Specific rules

lhapdf6ifcc.o: lhapdf6ifcc.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^  $(LHCXXFLAGS)

fastjet-addon.o: fastjet-addon.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^  $(FJCXXFLAGS)

fastjetD0RunIIConewrap.o: fastjetD0RunIIConewrap.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^  $(FJCXXFLAGS)

fastjetfortran-BBUY.o: fastjetfortran-BBUY.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^  $(FJCXXFLAGS)

fastjetfortran.o: fastjetfortran.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^  $(FJCXXFLAGS)

splitandfilt.o: splitandfilt.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^  $(FJCXXFLAGS)

pythia8F77.o: pythia8F77.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^  $(PY8CXXFLAGS)



