# ----------------------------------------------------------------------------
#
# Makefile for MODEL directory
# Jan 25 2006
#
# ----------------------------------------------------------------------------


#F77           = gfortran
FFLAGS        = $(XFFLAGS) -ffixed-line-length-132
LIBDIR        = ../$(OBJ)
LIBRARY	      = $(LIBDIR)/libmodel.a
MODEL         = printout.o couplings.o

VPATH=./:../$(OBJ)

%.o:%.f
	$(F77) $(FFLAGS) -c $< -o../$(OBJ)/$@


$(LIBRARY): $(MODEL)
	ar cru $(LIBRARY) $(patsubst %,../$(OBJ)/%,$(MODEL))
	ranlib $(LIBRARY)

smcalc: 
	$(F77) $(FFLAGS) SMCalc.f -o smcalc 
