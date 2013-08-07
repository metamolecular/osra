dnl
dnl Autoconf configuration for OSRA.
dnl
dnl Use "autoconf" to generate a configure script.
dnl
dnl Follow this link for complete autoconf macro list:
dnl http://www.gnu.org/software/autoconf/manual/html_node/Autoconf-Macro-Index.html
dnl

AC_INIT([osra], [2.0.0])
AC_PREREQ([2.65])

dnl
dnl Library version. If library ABI changes, increase the major number.
dnl

dnl Increment if interfaces have been added, removed or changed:
AC_SUBST([LIB_MAJOR_VERSION], [1])
dnl Increment if source code has changed:
AC_SUBST([LIB_MINOR_VERSION], [0])

dnl
dnl Platform detection
dnl

AC_CANONICAL_BUILD

AS_CASE(["$build_os"], [linux*], [build_os="linux"], [darwin*], [build_os="darwin"], [mingw*], [build_os="mingw"])

AC_LANG([C++])
AX_GNU_LD
AC_CXX_HAVE_STL

AS_IF([test "${ac_cv_cxx_have_stl}" != "yes"], [
	AC_MSG_ERROR([STL was not found; make sure you have installed libstdc++-dev package or check config.log])
])

dnl
dnl Checks for general purpose tools
dnl

AC_PROG_CXX
AC_PATH_PROG(RM, rm)
AC_PROG_LN_S
AC_PROG_INSTALL
AC_PATH_PROG(AR, ar)
AC_PROG_RANLIB

dnl
dnl Checks for programs needed to build a documentation
dnl

AC_ARG_ENABLE(
	[docs],
	[AC_HELP_STRING([--enable-docs], [Enable generating of documentation (man pages, API reference, etc)])],
	[
		dnl NaturalDocs executable is present as "/usr/bin/naturaldocs" in Debian and "/usr/bin/NaturalDocs" in Fedora: 
		AC_PATH_PROG(NATURALDOCS, naturaldocs)
		AS_IF([test -z "${NATURALDOCS}"], [
			AC_PATH_PROG(NATURALDOCS, NaturalDocs)
			dnl AS_IF([test -z "${NATURALDOCS}"], [
			dnl 	AC_MSG_ERROR([NaturalDocs executable was not found; make sure you have installed NaturalDocs package])
			dnl ])
		])
		
		AC_PATH_PROG(XSLTPROC, xsltproc)
		AS_IF([test -z "${XSLTPROC}"], [
			AC_MSG_ERROR([xsltproc executable was not found; make sure you have installed xsltproc package])
		])
	]
)

dnl
dnl Include "/usr/local" tree into the list of default locations
dnl

AS_IF([test -d "/usr/local/include"], [
	CPPFLAGS="-I/usr/local/include ${CPPFLAGS}"
])

AS_IF([test -d "/usr/local/lib"], [
	LDFLAGS="-L/usr/local/lib ${LDFLAGS}"
])

AS_IF([test -d "/opt/local/lib"], [
	MACPORTSLIBS="/opt/local/lib"
])

dnl
dnl OpenMP support
dnl

AC_ARG_ENABLE(
	[openmp],
	[AC_HELP_STRING([--enable-openmp], [Enable OpenMP extensions for compilers that support it])]
)

dnl Enable OpenMP on Linux by default (unless it is disabled):
AS_IF([test -z "${enable_openmp}" -a "$build_os" == "linux"], [enable_openmp=yes])

AS_IF([test "${enable_openmp}" == "yes"], [
	AC_OPENMP()
])

dnl
dnl Compilation/linking configuration options
dnl

AC_ARG_ENABLE(
	[static-linking],
	[AC_HELP_STRING([--enable-static-linking], [Try to link all libraries statically (fallback to dynamic linking if statis was not found) (GCC-only)])],
	[
		LDFLAGS="${LDFLAGS} -static-libgcc"
		dnl GNU LD-specific option:
		AS_IF([test "${ac_gnu_ld}" == "yes"], [LDFLAGS="${LDFLAGS} -Wl,-static"])
	]
)

AC_ARG_ENABLE(
	[profiling],
	[AC_HELP_STRING([--enable-profiling], [Include verbose profiling/debugging information into the binary])],
	[
		LDFLAGS="-pg ${LDFLAGS}"
		CXXFLAGS="-Wall -g3 -pg ${CXXFLAGS}"
	],
	[INSTALL_PROGRAM="${INSTALL_PROGRAM} -s"]
)

AC_ARG_WITH(
	[macports],
	[AC_HELP_STRING([--with-macports], [Alternative location for MacPorts libraries on OS X (default: /opt/local/lib)])],
	[MACPORTSLIBS="$withval"]
)

dnl Library compilation

AC_ARG_ENABLE(
	[lib],
	[AC_HELP_STRING([--enable-lib], [Include building of OSRA dynamic (.so) and static (.a) library into build cycle.])],
	[
		AC_SUBST([OSRA_LIB], [yes])
	]
)

dnl
dnl Java bridge support
dnl

AC_ARG_ENABLE(
	[java],
	[AC_HELP_STRING([--enable-java], [Include Java JNI bridge compilation into build cycle. Make sure you have JDK installed and JAVA_HOME environment variable pointing to its location or jni.h is in compiler search path.])]
)

AS_IF([test "${enable_java}" == "yes"], [
	AS_IF([test -n "${JAVA_HOME}"], [
		CPPFLAGS="-I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux -I${JAVA_HOME}/include/win32 ${CPPFLAGS}"
	])

	dnl Fix for jlong definition in jni.h on Windows
	AS_CASE(["$build_os"],
		[cygwin | mingw], [
			AC_DEFINE([__int64], [long long])
		]
	)

	AC_CHECK_HEADER([jni.h], [
		AC_TRY_COMPILE([
			#include <stddef.h>
			
			#include <jni.h>
		], [
			JNIEnv *j_env;
			jobjectArray j_obj_arr;
			jbyteArray j_byte_arr;
			jstring j_str;
			jboolean j_output_confidence;
			
			j_env->GetArrayLength(j_byte_arr);
			j_env->GetObjectArrayElement(j_obj_arr, 0);
			j_env->GetByteArrayElements(j_byte_arr, NULL);
			j_env->ReleaseByteArrayElements(j_byte_arr, (jbyte*) NULL, JNI_ABORT);
			j_env->GetStringUTFChars(j_str, NULL);
			j_env->ReleaseStringUTFChars(j_str, (const char *) NULL);
		], [
			AC_SUBST([OSRA_JAVA], [yes])
		], [
			AC_MSG_ERROR([Java JNI bridge support was requested, but failed because JNI API is not supported])
		])
	], [
		AC_MSG_ERROR([Java JNI bridge support was requested, but failed because "jni.h" was not found (search path: ${CPPFLAGS})])
	])
])

dnl
dnl Define includes paths and libraries paths
dnl

dnl tclap library

AX_PROBE_OBLIGATORY_HEADER(
	[tclap],
	[tclap/CmdLine.h],
	[auto],
	[custom location of tclap/CmdLine.h from TCLAP (Templatized C++ Command Line Parser Library), see http://tclap.sourceforge.net/]
)

dnl potrace library

AX_PROBE_OBLIGATORY_LIBRARY(
	[potrace],
	[potracelib.h],
	[auto],
	[custom location of potracelib.h from POTRACE (vector tracing library), see http://potrace.sourceforge.net/]
)

dnl gocr library

AX_PROBE_OBLIGATORY_LIBRARY(
	[gocr],
	[pgm2asc.h],
	[/usr/include/gocr /usr/local/include/gocr auto],
	[custom location of gocr/pgm2asc.h from GOCR (optical character recognition library), see http://jocr.sourceforge.net/]
)

dnl ocrad library

AX_PROBE_OBLIGATORY_LIBRARY(
	[ocrad],
	[ocradlib.h],
	[auto],
	[custom location of ocradlib.h from OCRAD (optical character recognition program), see http://www.gnu.org/software/ocrad/ocrad.html]
)

dnl tesseract library (optional)

AC_ARG_WITH(
	[tesseract-include],
	[AC_HELP_STRING([--with-tesseract-include], [Location of tesseract/baseapi.h header file])],
	[CPPFLAGS="-I${withval} ${CPPFLAGS}"]
)

AX_PROBE_OPTIONAL_LIBRARY(
	[tesseract],
	[tesseract/baseapi.h],
	[auto],
	[Enable optional support for Tesseract (OCR Engine developed at HP Labs), see http://code.google.com/p/tesseract-ocr/]
)

dnl cuneiform library (optional)

dnl AC_ARG_WITH(
dnl	[cuneiform-include],
dnl	[AC_HELP_STRING([--with-cuneiform-include], [Location of cuneiform.h header file])],
dnl	[CPPFLAGS="-I${withval} ${CPPFLAGS}"]
dnl)

AX_PROBE_OPTIONAL_LIBRARY(
	[cuneiform],
	[cuneiform.h],
	[auto],
	[Enable optional support for Cuneiform (OCR system developed by Cognitive technologies), see https://launchpad.net/cuneiform-linux]
)

dnl openbabel library

AX_PROBE_OBLIGATORY_LIBRARY(
	[openbabel],
	[openbabel/mol.h openbabel/obconversion.h openbabel/builder.h openbabel/alias.h],
	[/usr/include/openbabel-2.0 /usr/local/include/openbabel-2.0],
	[custom location of openbabel/builder.h from OpenBabel (open source chemistry toolbox), see http://openbabel.sourceforge.net/wiki/Main_Page]
)

dnl
dnl Checks for header files, libraries, types and functions
dnl

dnl limits.h header

osra_header_limits=yes

AC_CHECK_HEADER([limits.h], [
	AC_TRY_COMPILE([#include <limits.h>], [ int a = INT_MAX; ], [], [osra_header_limits=no])
], [osra_header_limits=no])

AS_IF([test "${osra_header_limits}" != "yes"], [
	AC_MSG_ERROR([limits.h was not found; make sure you have installed libc-dev package or check config.log])
])

dnl float.h header

osra_header_float=yes

AC_CHECK_HEADER([float.h], [
	AC_TRY_COMPILE([#include <float.h>], [ float f = FLT_MAX; ], [], [osra_header_float=no])
], [osra_header_float=no])

AS_IF([test "${osra_header_float}" != "yes"], [
	AC_MSG_ERROR([float.h was not found; make sure you have installed libc-dev package or check config.log])
])

dnl libgen.h header

osra_header_libgen=yes

AC_CHECK_HEADER([libgen.h], [
	AC_TRY_COMPILE([#include <libgen.h>], [ char *a = dirname(""); ], [], [osra_header_libgen=no])
], [osra_header_libgen=no])

AS_IF([test "${osra_header_libgen}" != "yes"], [
	AC_MSG_ERROR([libgen.h was not found; make sure you have installed libc-dev package or check config.log])
])

dnl math library

osra_lib_math=yes

AC_CHECK_HEADER([math.h], [
	dnl "-lm" should be checked only for non-OSX platforms:
	AS_IF([test "$build_os" != "darwin"], [
		AX_TRY_LINK([m], [#include <math.h>], [ sqrt(1); fabs(1); ], [], [osra_lib_math=no])
	])
], [osra_lib_math=no])

AS_IF([test "${osra_lib_math}" != "yes"], [
	AC_MSG_ERROR([math.h or libm was not found; make sure you have installed libc-dev package or check config.log])
])

dnl
dnl TCLAP library
dnl

osra_lib_tclap=yes

AC_TRY_COMPILE([
	#include <stdio.h> // Defines the EOF constant; needed for MinGW
	#include <stddef.h>
	#include <tclap/CmdLine.h>
], [
	TCLAP::CmdLine cmd("");
	TCLAP::ValueArg<int> a1("", "", "", false, 0, "");
	cmd.add(a1);
	TCLAP::SwitchArg a2("", "", "", false);
	cmd.add(a2);
	cmd.parse(0, NULL);
], [], [osra_lib_tclap=no])

AS_IF([test "${osra_lib_tclap}" != "yes"], [
	AC_MSG_ERROR([TCLAP API check failed; make sure you have installed libtclap-dev package or check config.log])
])

dnl
dnl Potrace library
dnl

osra_lib_potrace=yes

AX_TRY_LINK([potrace], [
	extern "C" {
	#include <potracelib.h>
	}
], [
	potrace_state_t state;
	potrace_param_t param;
	potrace_bitmap_t bitmap;
	potrace_path_t path;
	potrace_dpoint_t dpoint;

	potrace_trace(&param, &bitmap);
	potrace_state_free(&state);
	potrace_param_default();
], [], [osra_lib_potrace=no])

AS_IF([test "${osra_lib_potrace}" != "yes"], [
	AC_MSG_ERROR([POTRACE API check failed; make sure you have installed libpotrace-dev package or check config.log])
])

dnl
dnl Ocrad library
dnl

osra_lib_ocrad=yes

AX_TRY_LINK([ocrad], [
	#include <ocradlib.h>
], [
	struct OCRAD_Pixmap opix;
	
	OCRAD_Descriptor * const ocrdes = OCRAD_open();

	OCRAD_get_errno(ocrdes);
	OCRAD_set_image(ocrdes, &opix, 0);
	OCRAD_scale(ocrdes, 2);
	OCRAD_recognize(ocrdes, 0);
	OCRAD_result_first_character(ocrdes);
	OCRAD_close(ocrdes);
], [], [osra_lib_ocrad=no])

AS_IF([test "${osra_lib_ocrad}" != "yes"], [
	AC_MSG_ERROR([OCRAD API check failed; make sure you have installed libocrad-dev package or check config.log])
])

dnl
dnl Tesseract library (optional)
dnl

AS_IF([test "${ac_lib_tesseract}" == "yes"], [
	dnl Note: This hack to link statically against Tesseract is not portable: the list of libraries should be discovered via pkg-config.
	AS_IF([test "${enable_static_linking:+set}" == "set"], [
		LIBS="-ljpeg -lz ${LIBS}"
		AC_CHECK_LIB([pthread], [sem_init], [], [ac_lib_tesseract=no])
		AC_CHECK_LIB([tiff], [TIFFGetField], [], [ac_lib_tesseract=no])
	])

	AX_TRY_LINK([tesseract], [
		#include <stddef.h>
		#include <tesseract/baseapi.h>
		using namespace tesseract;
	], [
		TessBaseAPI tess;
		tess.Init(NULL, "");
		char* text = tess.TesseractRect((const unsigned char*) NULL, 1, 0, 0, 0, 0, 0);
		tess.End();
	], [
		AC_DEFINE([HAVE_TESSERACT_LIB])
		AC_SUBST([TESSERACT_LIB], [yes])
	], [ac_lib_tesseract=no])
])

AS_IF([test "${with_tesseract}" != "no" -a "${ac_lib_tesseract}" == "no"], [
	AC_MSG_WARN([Tesseract library support was requested, but disabled due to API check failure (see above messages for the reason or check config.log)])
])

dnl
dnl Cuneiform library (optional)
dnl

AS_IF([test "${ac_lib_cuneiform}" == "yes"], [
	AX_TRY_LINK([cuneiform], [
		#include <stddef.h>
		#include <cuneiform.h>
	], [
		PUMA_XOpen(NULL, "");
		PUMA_XFinalRecognition();
		PUMA_XClose();
		PUMA_Init(0, NULL);
		PUMA_SetImportData(PUMA_Word32_Language, NULL);
		PUMA_SaveToMemory(NULL, PUMA_TOTEXT, PUMA_CODE_ASCII, NULL, 1);
		PUMA_Done();
	], [AC_DEFINE([HAVE_CUNEIFORM_LIB])], [ac_lib_cuneiform=no])
])

AS_IF([test "${with_cuneiform}" != "no" -a "${ac_lib_cuneiform}" == "no"], [
	AC_MSG_WARN([Cineiform library support was requested, but disabled due to API check failure (see above messages for the reason or check config.log)])
])

dnl
dnl OpenBabel library
dnl

osra_lib_openbabel=yes

AS_IF([test "$build_os" == "mingw"], [LIBS="-lws2_32 ${LIBS}"])

dnl These libraries are needed for static linking with openbabel:
AS_IF([test "${enable_static_linking:+set}" == "set"], [
	AS_IF([test "$build_os" == "darwin"], [
		dnl Special case for OSX: the libraries should be at fixed locations:
		AC_CHECK_LIB([inchi], [GetINCHI])
		LIBS="${MACPORTSLIBS}/libxml2.a ${MACPORTSLIBS}/libiconv.a ${MACPORTSLIBS}/libz.a ${LIBS}"
	], [
		dnl Note: This hack to link statically against OpenBabel is not portable: the list of libraries should be discovered via pkg-config.
		AC_CHECK_LIB([inchi], [GetINCHI])
		AC_CHECK_LIB([z], [inflateEnd])
		AC_CHECK_LIB([xml2], [xmlParseFile])
	])
])

AX_TRY_LINK([openbabel], [
	#include <openbabel/mol.h>
	#include <openbabel/obconversion.h>
	#include <openbabel/builder.h>
	#include <openbabel/alias.h>
], [
	OpenBabel::OBPairData label;
	OpenBabel::AliasData aliasData;
	OpenBabel::OBAtom atom;
	OpenBabel::OBMol mol;
	OpenBabel::OBBond bond;
	OpenBabel::OBRing ring;
	OpenBabel::OBConversion conversion;

	label.SetAttribute("");
	label.SetValue("");
	
	aliasData.SetAlias("");
	aliasData.SetOrigin(OpenBabel::external);
	aliasData.Expand(mol, 1);

	atom.IsHeteroatom();
	atom.IsCarbon();
	atom.IsNitrogen();
	atom.IsOxygen();
	atom.IsSulfur();
	atom.IsChiral();
	atom.GetAtomicNum();
	atom.SetData(&aliasData);
	atom.GetFormalCharge();
	atom.SetVector(0, 0, 0);
	atom.SetAtomicNum(0);
	atom.SetFormalCharge(0);
	atom.SetIdx(0);
	atom.Clear();

	mol.NumAtoms();
	mol.NumBonds();
	mol.NumRotors();
	mol.GetAtom(0);
	mol.GetBond(0);
	mol.GetFirstAtom();
	mol.GetSSSR();
	mol.NewAtom();
	mol.AddBond(0, 0, 0, 0);
	mol.AssignSpinMultiplicity(false);
	mol.SetDimension(0);
	mol.SetData(&label);
	std::vector< std::vector<int> > v;
	mol.ContigFragList(v);
	mol.ConnectTheDots();
	mol.Separate();
	mol.StripSalts(0);
	mol.BeginModify();
	mol.EndModify();
	mol.FindRingAtomsAndBonds();

	bond.GetBondOrder();
	bond.GetFlags();
	bond.SetBondOrder(0);
	bond.IsInRing();
	bond.IsDouble();
	bond.IsTriple();
	bond.UnsetUp();
	bond.UnsetDown();
	bond.UnsetAromatic();
	bond.UnsetHash();
	bond.UnsetWedge();
	bond.GetBeginAtom();

	ring.IsAromatic();

	conversion.SetInFormat("");
	conversion.SetOutFormat("");
	conversion.ReadString(&mol, "");
	conversion.Read(&mol);
	conversion.WriteString(&mol, false);

], [], [osra_lib_openbabel=no])

AS_IF([test "${osra_lib_openbabel}" != "yes"], [
	AC_MSG_ERROR([OpenBabel API check failed; make sure you have installed libopenbabel-dev package or check config.log])
])

dnl
dnl GraphicsMagick library
dnl

AC_ARG_ENABLE(
	[graphicsmagick-config],
	[AC_HELP_STRING([--disable-graphicsmagick-config], [Disable use of GraphicsMagick++-config for library linking])]
)

dnl The variable ${enable_graphicsmagick_config} has a reverse meaning: when set it means "disable-graphicsmagick-config"
AS_IF([test "${enable_graphicsmagick_config:+set}" != "set"], [
	AC_PATH_PROG(GRAPHICS_MAGICK_CONFIG, GraphicsMagick++-config)

	AS_IF([test "${GRAPHICS_MAGICK_CONFIG}" == ""], [
		AC_MSG_ERROR([GraphicsMagick++-config was not found; make sure you have installed libgraphicsmagick++-dev package])
	])	

	CPPFLAGS="`${GRAPHICS_MAGICK_CONFIG} --cppflags` ${CPPFLAGS}"

	AS_IF([test "${enable_static_linking:+set}" == "set"], [
		dnl Here goes statically linking setup section.
		dnl Note: This hack to link statically against GraphicsMagick is not portable: the list of libraries should be discovered via pkg-config.
		LIBS="`${GRAPHICS_MAGICK_CONFIG} --libs` ${LIBS}"
		AS_IF([test "$build_os" == "cygwin"],
			[LIBS="-lxcb -lXau -lXdmcp ${LIBS}"]
		)
	])
], [
	AS_IF([test "${enable_static_linking:+set}" == "set" -a "$build_os" == "darwin"], [
		dnl Special case for OSX: the libraries should be at fixed locations:
		LIBS="-ljpeg ${MACPORTSLIBS}/liblcms.a ${MACPORTSLIBS}/libtiff.a ${MACPORTSLIBS}/libfreetype.a ${MACPORTSLIBS}/libjasper.a ${MACPORTSLIBS}/libpng.a ${MACPORTSLIBS}/libbz2.a ${MACPORTSLIBS}/libz.a ${LIBS}"
	])
])

_AS_ECHO_LOG([CPPFLAGS="${CPPFLAGS}" LIBS="${LIBS}" after populated with GraphicsMagick++-config])

dnl We cannot move this check to the top, because includes for GraphicsMagick are defined only at this point:

AX_PROBE_OBLIGATORY_LIBRARY(
	[graphicsmagick],
	[Magick++.h],
	[/usr/include/GraphicsMagick /usr/local/include/GraphicsMagick],
	[custom location of Magick++.h from GraphicsMagick (image manipulation library), see http://www.graphicsmagick.org/]
)

osra_lib_graphicsmagick=yes

AX_TRY_LINK([GraphicsMagick++ GraphicsMagick], [
	#include <Magick++.h>
], [
	Magick::Image image;
	Magick::Color color;
	Magick::ColorRGB colorRGB;
	Magick::ColorGray colorGray;

	image.getPixels(0, 0, 0, 0);
	image.pixelColor(0, 0);
	image.columns();
	image.rows();
	image.matte();
	image.backgroundColor(color);
	image.rotate(0);

	color.alpha();

	colorRGB.red();
	colorRGB.green();
	colorRGB.blue();

	colorGray.shade();

	MagickLib::InitializeMagick((char*) 0);
	MagickLib::DestroyMagick();
], [], [osra_lib_graphicsmagick=no])

AS_IF([test "${osra_lib_graphicsmagick}" != "yes"], [
	AC_MSG_ERROR([GraphicsMagick API check failed; make sure you have installed libgraphicsmagick++-dev and all dependent -dev packages or check config.log])
])

dnl netpbm library check (which is used by GOCR)
dnl (copied from GOCR configure.in)

dnl The following symbols are used by GOCR:
dnl for i in /usr/lib64/* /usr/lib/*; do nm $i 2>/dev/null | egrep '(pnm_readpaminit|pnm_allocpamrow|pnm_readpamrow|pm_freerow|pnm_nextimage)' && echo "[$i]"; done
AC_SEARCH_LIBS([pnm_readpnminit], [netpbm pnm])

dnl gocr library
dnl 
dnl TODO: This check should go the last in the row, because of missing global OCR_JOB declaration:
dnl 
dnl pgm2asc.c:2900: In function `add_line_info':
dnl pgm2asc.c:1712: In function `count_subboxes'
dnl pgm2asc.c:2505: In function `divide_vert_glued_boxes'
dnl pgm2asc.c:1764: In function `glue_holes_inside_chars'
dnl pgm2asc.c:1146: In function `distance'
dnl undefined reference to `JOB'
dnl
dnl All other checks after this one will fail. 

osra_lib_gocr=yes

AX_TRY_LINK([Pgm2asc], [
	extern "C" {
	#include <pgm2asc.h>
	#include <gocr.h>
	}
	job_t *OCR_JOB, *JOB;
], [
	job_t job;

	job_init(&job);
	job_init_image(&job);
	job_free_image(&job);
	pgm2asc(&job);
], [], [osra_lib_gocr=no])

AS_IF([test "${osra_lib_gocr}" != "yes"], [
	AC_MSG_ERROR([GOCR API check failed; make sure you have installed libgocr-dev package or check config.log])
])

dnl
dnl Generate Makefile, config.h, etc
dnl

dnl Compute the library patch version from the package version:
AC_SUBST([LIB_PATCH_VERSION], [`echo ${PACKAGE_VERSION} | perl -ne '/(\d+)\.(\d+)\.(\d+)/ || die; printf "%d%02d%02d", $1, $2, $3;'`])

AS_CASE(["$build_os"],
	[darwin], [
		AC_SUBST([LDSHAREDFLAGS], ['-dynamiclib -Wl,-dylib_install_name -Wl,$@.$(LIB_MAJOR_VERSION)'])
		AC_SUBST([SHAREDEXT], [.dylib])
	],
	[cygwin | mingw], [
		AC_SUBST([LDSHAREDFLAGS], ['-shared -Wl,--kill-at -Wl,--out-implib,$@.a'])
		AC_SUBST([SHAREDEXT], [.dll])
	],
	[*], [
		dnl Assume this is Linux/Unix:
		AC_SUBST([LDSHAREDFLAGS], ['-shared -Wl,-soname,$@.$(LIB_MAJOR_VERSION)'])
		AC_SUBST([SHAREDEXT], [.so])
	]
)

AS_IF([test -n "${OPENMP_CXXFLAGS}" ], [
	CPPFLAGS="${OPENMP_CFLAGS} ${CPPFLAGS}"
	CXXFLAGS="${OPENMP_CXXFLAGS} ${CXXFLAGS}"
])

dnl Substitute config.h paths correctly:
test "$prefix" = "NONE" && prefix=$ac_default_prefix

dnl This will resolve 2st-level reference (${datarootdir} and ${PACKAGE_TARNAME}):
datadir=`eval echo "${datadir}"`
docdir=`eval echo "${docdir}"`
dnl This will resolve 1st-level reference (${prefix}):
resolved_datadir=`eval echo "${datadir}"`

dnl C-code needs everything completely resolved:
AC_DEFINE_UNQUOTED(DATA_DIR, ["$resolved_datadir"])
AC_SUBST([resolved_datadir], [$resolved_datadir])

AC_CONFIG_HEADER([src/config.h])
AC_CONFIG_FILES([Makefile Makefile.inc doc/manual.sgml package/linux/osra.pc package/linux/debian/control package/linux/debian/rules package/linux/suse/osra.spec pom.xml])
AC_OUTPUT
AC_MSG_NOTICE([Now you can run "make all install" to compile and install the application.])
