#
# This makefile is responsible for building the executable.
#

include ../Makefile.inc
include Makefile.dep

.PHONY: clean_obj

LIB_VERSION	:= $(LIB_MAJOR_VERSION).$(LIB_MINOR_VERSION).$(LIB_PATCH_VERSION)

TARGETS		:= osra$(EXEEXT)

OBJ_LIB		:=  osra_lib.o osra_grayscale.o osra_fragments.o osra_segment.o osra_labels.o osra_thin.o osra_common.o osra_stl.o osra_structure.o osra_anisotropic.o osra_ocr.o osra_openbabel.o mcdlutil.o unpaper.o osra_reaction.o

ifdef TESSERACT_LIB
OBJ_LIB		+= osra_ocr_tesseract.o
endif

OBJ_CLI		:= $(OBJ_LIB) osra.o

ifdef OSRA_LIB
TARGETS		+= libosra.a libosra$(SHAREDEXT)
endif

ifdef OSRA_JAVA
OBJ_JAVA	:= $(OBJ_LIB) osra_java.o
TARGETS		+= libosra_java$(SHAREDEXT)
endif

ifdef OSRA_ANDROID
OBJ_ANDROID	:= $(OBJ_LIB) osra.o
endif

all:
# From here: http://stackoverflow.com/questions/5584872/complex-conditions-check-in-makefile/5586785#5586785
ifneq ($(or $(OSRA_LIB),$(OSRA_JAVA),$(OSRA_ANDROID)),)
	$(RM) -f $(OBJ_CLI)
endif
	$(MAKE) osra$(EXEEXT)
ifdef OSRA_LIB
	$(RM) -f $(OBJ_LIB)
	$(MAKE) libosra.a
	$(RM) -f $(OBJ_LIB)
	$(MAKE) libosra$(SHAREDEXT)
endif
ifdef OSRA_JAVA
	$(RM) -f $(OBJ_JAVA)
	$(MAKE) libosra_java$(SHAREDEXT)
endif
ifdef OSRA_ANDROID
	$(RM) -f $(OBJ_ANDROID)
	$(MAKE) libosra_andriod$(SHAREDEXT)
endif
#	We have to update the timestaps of the targets, otherwise "install" target will try to re-link and will cause missed symbols:
	touch $(TARGETS)

osra$(EXEEXT): $(OBJ_CLI)
	$(LINK.cpp) -o $@ $(OBJ_CLI) $(LIBS)

ifdef OSRA_LIB
libosra.a: CXXFLAGS += -DOSRA_LIB
libosra.a: $(OBJ_LIB)
	$(AR) cru $@ $(OBJ_LIB)
	$(RANLIB) $@

libosra$(SHAREDEXT): CXXFLAGS += -fPIC -DOSRA_LIB
libosra$(SHAREDEXT): $(OBJ_LIB)
	$(LINK.cpp) $(LDSHAREDFLAGS) -o $@ $(OBJ_LIB) $(LIBS)
endif

ifdef OSRA_JAVA
libosra_java$(SHAREDEXT): CXXFLAGS += -fPIC -DOSRA_LIB -DOSRA_JAVA
libosra_java$(SHAREDEXT): $(OBJ_JAVA)
	$(LINK.cpp) $(LDSHAREDFLAGS) -o $@ $(OBJ_JAVA) $(LIBS)
endif

ifdef OSRA_ANDROID
libosra_andriod$(SHAREDEXT): CXXFLAGS += -fPIC -DOSRA_LIB -DOSRA_ANDROID
libosra_andriod$(SHAREDEXT): $(OBJ_ANDROID)
	$(LINK.cpp) $(LDSHAREDFLAGS) -o $@ $(OBJ_ANDROID) $(LIBS)
endif

Makefile.dep: ../Makefile.inc $(patsubst %.o,%.cpp,$(OBJ_CLI))
	$(CXX) $(CPPFLAGS) -MM $^ > Makefile.dep 

# Correct installation for Cygwin/MinGW also needs handling of libosra.dll.a, which is not done here:
install: $(TARGETS)
	$(INSTALL_DIR) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) osra$(EXEEXT) $(DESTDIR)$(bindir)
ifdef OSRA_LIB
	$(INSTALL_DIR) $(DESTDIR)$(libdir) $(DESTDIR)$(includedir) $(DESTDIR)$(libdir)/pkgconfig
	$(INSTALL_PROGRAM) libosra$(SHAREDEXT) $(DESTDIR)$(libdir)/libosra$(SHAREDEXT).$(LIB_VERSION)
	$(LN_S) -f libosra$(SHAREDEXT).$(LIB_VERSION) $(DESTDIR)$(libdir)/libosra$(SHAREDEXT).$(LIB_MAJOR_VERSION)
	$(LN_S) -f libosra$(SHAREDEXT).$(LIB_MAJOR_VERSION) $(DESTDIR)$(libdir)/libosra$(SHAREDEXT)
	$(INSTALL_DATA) libosra.a $(DESTDIR)$(libdir)
	$(INSTALL_DATA) osra_lib.h $(DESTDIR)$(includedir)
	$(INSTALL_DATA) ../package/linux/osra.pc $(DESTDIR)$(libdir)/pkgconfig
endif
ifdef OSRA_JAVA
	$(INSTALL_DIR) $(DESTDIR)$(libdir)
	$(INSTALL_PROGRAM) libosra_java$(SHAREDEXT) $(DESTDIR)$(libdir)/libosra_java$(SHAREDEXT).$(LIB_VERSION)
	$(LN_S) -f libosra_java$(SHAREDEXT).$(LIB_VERSION) $(DESTDIR)$(libdir)/libosra_java$(SHAREDEXT).$(LIB_MAJOR_VERSION)
	$(LN_S) -f libosra_java$(SHAREDEXT).$(LIB_MAJOR_VERSION) $(DESTDIR)$(libdir)/libosra_java$(SHAREDEXT)
endif

# "install" and "make" tools auomatically autodetect the extension for executables, but "rm" needs extension correction.
uninstall:
	-$(RM) -f \
		$(DESTDIR)$(bindir)/osra$(EXEEXT) \
		$(DESTDIR)$(libdir)/libosra$(SHAREDEXT).$(LIB_VERSION) \
		$(DESTDIR)$(libdir)/libosra$(SHAREDEXT).$(LIB_MAJOR_VERSION) \
		$(DESTDIR)$(libdir)/libosra$(SHAREDEXT) \
		$(DESTDIR)$(libdir)/libosra.a \
		$(DESTDIR)$(libdir)/libosra_java$(SHAREDEXT).$(LIB_VERSION) \
		$(DESTDIR)$(libdir)/libosra_java$(SHAREDEXT).$(LIB_MAJOR_VERSION) \
		$(DESTDIR)$(libdir)/libosra_java$(SHAREDEXT)
		$(DESTDIR)$(includedir)/osra_lib.h \
		$(DESTDIR)$(libdir)/pkgconfig/osra.pc

clean:
	-$(RM) -f *.o osra$(EXEEXT) libosra*.*

distclean: clean
	-$(RM) -f config.h Makefile.dep

../Makefile.inc: ../Makefile.inc.in ../config.status
	cd .. && ./config.status
