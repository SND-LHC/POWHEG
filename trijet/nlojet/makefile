ifdef CDEBUG
DEBUG= -g -O0
endif

LIBDIR        = ../$(OBJ)
LIBRARY	      = $(LIBDIR)/libnlojet.a
CODE          =  $(patsubst %.cc,%.o,$(wildcard *.cc))


VPATH=./:../$(OBJ)

INCLUDE1=$(PWD)/bits
INCLUDE2=$(PWD)/include

%.o: %.cc | 
	$(CC) $(DEBUG)  -c  -o  ../$(OBJ)/$@ $^ $(CXXFLAGS) -I. -I$(INCLUDE1) -I$(INCLUDE2)


all: $(LIBRARY)

$(LIBRARY): $(CODE)
	 ar cru $(LIBRARY) $(patsubst %,../$(OBJ)/%, $(CODE)); ranlib $(LIBRARY)

