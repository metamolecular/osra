CXX		:= g++
LD		:= g++

CXXFLAGS	:= -g3 -ggdb -O0 -DOSRA_LIB
CPPFLAGS	:= -I../../src
LDFLAGS		:= -L../../src -L/usr/local/lib

LIBS		:= -losra

.PHONY: all clean

.SUFFIXES: .c .cpp

OBJ			= lib_sample.o

all: lib_sample

#	LD_LIBRARY_PATH=../../src ./lib_sample /usr/local/lib

lib_sample: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o lib_sample
