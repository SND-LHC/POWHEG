
#F77           = gfortran
FFLAGS        = $(XFFLAGS)
LIBDIR        = ../$(OBJ)
LIBRARY	      = $(LIBDIR)/libmadgraph.a
CODE          =  sborn_proc.o $(patsubst %.f,%.o,$(wildcard sborn_sf_???.f)) $(patsubst %.f,%.o,$(wildcard sborn_cl_???.f)) $(patsubst %.f,%.o,$(wildcard b_sf_???_???.f)) $(patsubst %.f,%.o,$(wildcard b_cl_???_???.f)) nexper.o

VPATH=./:../$(OBJ)

%.o: %.f
	$(F77) $(FFLAGS) -c $*.f -o../$(OBJ)/$@

$(LIBRARY): $(CODE)
	ar cru $(LIBRARY) $(patsubst %,../$(OBJ)/%,$(CODE))
	ranlib $(LIBRARY)


