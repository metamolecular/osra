CXX		:= g++
LD		:= g++
RM		:= /bin/rm

CXXFLAGS	:= -g3 -O2
CPPFLAGS	:= -I/usr/include/GraphicsMagick
LDFLAGS		:= -L/usr/lib

LIBS		:= -lGraphicsMagick++ 

.PHONY: all clean

.SUFFIXES: .c .cpp

OBJ		= test.o

all: test

test: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

clean:
	$(RM) -f *.o test
