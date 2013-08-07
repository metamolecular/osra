CXX		:= g++
LD		:= g++
RM		:= /bin/rm

CXXFLAGS	:= -g3 -ggdb -O0
CPPFLAGS	:= 
LDFLAGS		:= -L/usr/lib

LIBS		:= -ltesseract_api

.PHONY: all clean

.SUFFIXES: .c .cpp

OBJ		= test.o

all:
	rm -f *.o
	$(MAKE) test_local
	rm -f *.o
	$(MAKE) test_global
	./test_global > 1
	./test_local > 2
	diff -au 1 2

test_local: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

test_global: CXXFLAGS += -DTESS_GLOBAL_INSTANCE
test_global: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

clean:
	$(RM) -f *.o test_*
