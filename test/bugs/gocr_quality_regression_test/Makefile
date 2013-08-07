CXX		:= g++
LD		:= g++

CXXFLAGS	:= -g3 -O2
CPPFLAGS	:= -I/usr/include/gocr
LDFLAGS		:= -L/usr/lib
LIBS		:= -lPgm2asc -lnetpbm

OBJ		= osra_gocr.o

.PHONY: all clean

.SUFFIXES: .c .cpp

all: test

test: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

clean:
	$(RM) -f *.o test
