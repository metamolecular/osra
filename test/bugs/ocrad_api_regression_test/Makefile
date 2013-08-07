CXX		:= g++
LD		:= g++

CXXFLAGS	:= -g3 -O2
CPPFLAGS	:= -I/usr/include/ocrad
LDFLAGS		:= -L/usr/lib
LIBS		:= -locrad

OBJ		= osra_ocr.o

.PHONY: all clean

.SUFFIXES: .c .cpp

all: test

test: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

clean:
	$(RM) -f *.o test
