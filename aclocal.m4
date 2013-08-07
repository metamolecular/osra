m4_include([m4/ac_cxx_namespaces.m4])
m4_include([m4/ac_cxx_have_stl.m4])

# SYNOPSIS
#
# AX_GNU_LD()
#
# DESCRIPTION
#
# This macro sets the variable "ac_gnu_ld" to "yes" if GNU linker is present in the system.
# Taken from http://fink.sourceforge.net/files/ltconfig
#
AC_DEFUN([AX_GNU_LD], [
	AC_PATH_PROG(LD, ld)
	AC_MSG_CHECKING([if ld ($LD) is GNU ld])
	AS_IF(["$LD" -v 2>&1 < /dev/null | egrep '(GNU|with BFD)' > /dev/null], [
		AS_ECHO(yes)
		ac_gnu_ld=yes
	], [
		AS_ECHO(no)
		ac_gnu_ld=no
	])
])

# AX_PROBE_OBLIGATORY_HEADER(lib-name, headers ..., default-paths-to-check ..., help-string)
#
# DESCRIPTION
#
# This macro defines the helper argument "--with-${lib-name}" using AC_ARG_WITH and
# probes the given (as 2nd argument) headers first in system-wide locations and
# then for the specified space-separated locations (given as 3rd argument).
# If probing fails, the error is reported.
# See AX_PROBE_LIBRARY for more details about the probing itself.
#
AC_DEFUN([AX_PROBE_OBLIGATORY_HEADER], [
	dnl m4_if() macro on some reason does not work inside AC_HELP_STRING(): 
	AC_ARG_WITH(
		[$1-include],
		[m4_if([$3], [], [AC_HELP_STRING([--with-$1-include], [$4])], [AC_HELP_STRING([--with-$1-include], [$4 (default: "$3")])])],
		[
			with_$1="${with_$1_include}"
			CPPFLAGS="-I${withval} ${CPPFLAGS}"
		],
		[with_$1_include="$3"])

	AS_IF([test "${with_$1_include}" == "no"], [AC_MSG_ERROR([The library $1 is obligatory. You cannot disable it.])])

	dnl Here the value of ${with_$1} is either:
	dnl * if option was given in a command-line, it's value (if empty, then only system paths are checked)
	dnl * if option was omitted, the defaults ($3) are used 
	AX_PROBE_LIBRARY([$1], [$2])

	AS_IF([test "${ac_lib_$1}" != "yes"], [AC_MSG_ERROR([$2 header(s) is missing. Check the default/listed above headers locations.])])
]) # AX_PROBE_OBLIGATORY_HEADER

# SYNOPSIS
#
# AX_PROBE_OBLIGATORY_LIBRARY(lib-name, headers ..., default-paths-to-check ..., help-string)
#
# DESCRIPTION
#
# This macro defines the helper argument "--with-${lib-name}" using AC_ARG_WITH and
# probes the given (as 2nd argument) headers first in system-wide locations and
# then for the specified space-separated locations (given as 3rd argument).
# If probing fails, the error is reported.
# See AX_PROBE_LIBRARY for more details about the probing itself.
#
AC_DEFUN([AX_PROBE_OBLIGATORY_LIBRARY], [
	dnl m4_if() macro on some reason does not work inside AC_HELP_STRING(): 
	AC_ARG_WITH(
		[$1-include],
		[m4_if([$3], [], [AC_HELP_STRING([--with-$1-include], [$4])], [AC_HELP_STRING([--with-$1-include], [$4 (default: "$3")])])],
		[
			with_$1="${with_$1_include}"
			CPPFLAGS="-I${withval} ${CPPFLAGS}"
		],
		[with_$1="$3"])

	AS_IF([test "${with_$1}" == "no"], [AC_MSG_ERROR([The library $1 is obligatory. You cannot disable it.])])

	AC_ARG_WITH(
		[$1-lib],
		[AC_HELP_STRING([--with-$1-lib], [custom location of the library])],
		[LDFLAGS="-L${withval} ${LDFLAGS}"])

	dnl Here the value of ${with_$1} is either:
	dnl * if option was given in a command-line, it's value (if empty, then only system paths are checked)
	dnl * if option was omitted, the defaults ($3) are used 
	AX_PROBE_LIBRARY([$1], [$2])

	AS_IF([test "${ac_lib_$1}" != "yes"], [AC_MSG_ERROR([$2 header(s) is missing. Check the default/listed above headers locations.])])
]) # AX_PROBE_OBLIGATORY_LIBRARY


# SYNOPSIS
#
# AX_PROBE_OPTIONAL_LIBRARY(lib-name, headers ..., default-paths-to-check ..., help-string)
#
# DESCRIPTION
#
# This macro defines the helper argument "--with-${lib-name}" using AC_ARG_WITH and
# probes the given (as 2nd argument) headers first in system-wide locations and
# then for the specified space-separated locations (given as 3rd argument).
# Behaves the same way as AX_PROBE_OBLIGATORY_LIBRARY, but does not report
# the error if the library was not found.
# See AX_PROBE_LIBRARY for more details about the probing itself.
#
AC_DEFUN([AX_PROBE_OPTIONAL_LIBRARY], [
	AC_ARG_WITH(
		[$1],
		[m4_if([$3], [], [AC_HELP_STRING([--with-$1], [$4 (optional)])], [AC_HELP_STRING([--with-$1], [$4 (optional) (default: "$3")])])],
		[],
		[with_$1="no"])

	AS_IF([test "${with_$1}" != "no"], [
		AS_IF([test "${with_$1}" == "" -o "${with_$1}" == "yes"], [
			with_$1="$3"
		])

	AC_ARG_WITH(
		[$1-include],
		[m4_if([$3], [], [AC_HELP_STRING([--with-$1-include], [$4])], [AC_HELP_STRING([--with-$1-include], [$4 (default: "$3")])])],
		[
			with_$1="${with_$1_include}"
			CPPFLAGS="-I${withval} ${CPPFLAGS}"
		],
		[with_$1="$3"])

		AC_ARG_WITH(
			[$1-lib],
			[AC_HELP_STRING([--with-$1-lib], [custom location of the library])],
			[LDFLAGS="-L${withval} ${LDFLAGS}"])

		dnl Here the value of ${with_$1} is:
		dnl * if option was given in a command-line, it's value (if empty, the defaults ($3) are used)
		dnl * if option was omitted (or --without-$1 form was used), this block is not executed  
		AX_PROBE_LIBRARY([$1], [$2])
	])])
]) # AX_PROBE_OPTIONAL_LIBRARY


# SYNOPSIS
#
# AX_PROBE_LIBRARY(lib-name, headers ...)
#
# DESCRIPTION
#
# This macro probes the given (as 2nd argument) headers first in system-wide locations and
# then for the specified space-separated locations (given as ${with_${lib-name}} variable). If
# probing succeeds it adds the headers location to $CPPFLAGS and library locations
# to $LDFLAGS. The special path "auto" means that the library will be autoprobed
# in user's $HOME. If headers have been located this macro defines the variable
# $ac_lib_${lib-name} (which is set to "yes") and also in case of "auto" location
# appends the found dirs with headers to $CPPFLAGS and appends the directories
# with library binaries to $LDFLAGS.
#
# TODO: If directories contain spaces this will cause problems (fixing $IFS will
# cause problems in other places).
#
AC_DEFUN([AX_PROBE_LIBRARY], [
	dnl Testing for default locations, ignoring the optional locations:   
	AC_CHECK_HEADERS([$2], [ac_lib_$1=yes], [ac_lib_$1=no])

	dnl Testing the specified locations:   
	AS_IF([test "${ac_lib_$1}" != "yes" -a "${with_$1}" != ""], [
		AX_RESET_HEADERS_CACHE([$2])

		AS_FOR([], [ac_test_location], [${with_$1}], [
			dnl Probing the library in user's $HOME:
			AS_IF([test "${ac_test_location}" = "auto"], [
				dnl Read the directory entries by mask sorted alphabetically in reverse order:
				AS_FOR([], [ac_location], [`ls -1d $HOME/$1-* 2>/dev/null | tac`], [ 
					AS_IF([test -d "${ac_location}"], [
						dnl Save the current state
						ax_probe_library_save_LDFLAGS=${LDFLAGS}
						ax_probe_library_save_CPPFLAGS=${CPPFLAGS}

						dnl Compose the list of unique locations of headers:
						AS_FOR([], [ac_inc_location], [`find "${ac_location}" -iname '*.h' |
							while read ac_include_location; do dirname "${ac_include_location}"; done |
								sort -u`], [
							CPPFLAGS="-I${ac_inc_location} ${CPPFLAGS}"
						])

						dnl Compose the list of unique locations of libraries (standard library extensions are taken from autoconf/libs.m4:185):
						AS_FOR([], [ac_lib_location], [`find "${ac_location}" -iname '*.so' -o -iname '*.sl' -o -iname '*.dylib' -o -iname '*.a' -o -iname '*.dll' |
							while read ac_library_location; do dirname "${ac_library_location}"; done |
								sort -u`], [
							LDFLAGS="-L${ac_lib_location} ${LDFLAGS}"
						])

						AC_MSG_CHECKING([$1 for $2 in ${ac_location}])
						AS_ECHO()
						_AS_ECHO_LOG([CPPFLAGS="${CPPFLAGS}" and LDFLAGS="${LDFLAGS}" for HOME location check])

						AC_CHECK_HEADERS([$2], [ac_lib_$1=yes], [ac_lib_$1=no])

						dnl We have found the location, leave the loop:
						AS_IF([test "${ac_lib_$1}" = "yes"], [break 2])

						dnl Restore the state to original in case of unsuccessful attempt
						LDFLAGS=${ax_probe_library_save_LDFLAGS}
						CPPFLAGS=${ax_probe_library_save_CPPFLAGS}
						AX_RESET_HEADERS_CACHE([$2])
					])
				])
			], [
				dnl Save the current state
				ax_probe_library_save_CPPFLAGS=${CPPFLAGS}

				CPPFLAGS="-I${ac_test_location} $CPPFLAGS"

				AC_MSG_CHECKING([$1 for $2 in ${ac_test_location}])
				AS_ECHO()
				_AS_ECHO_LOG([CPPFLAGS="${CPPFLAGS}" for custom location check])

				AC_CHECK_HEADERS([$2], [ac_lib_$1=yes], [ac_lib_$1=no])

				dnl We have found the location, leave the loop:
				AS_IF([test "${ac_lib_$1}" = "yes"], [break])

				dnl Restore the state to original in case of unsuccessful attempt
				CPPFLAGS=${ax_probe_library_save_CPPFLAGS}
				AX_RESET_HEADERS_CACHE([$2])
			])
		])
	])
]) # AX_PROBE_LIBRARY

# SYNOPSIS
#
# AX_TRY_LINK(library, includes, function-body [, action-if-true [, action-if-false]])
#
# DESCRIPTION
#
# This macro is a combination of autoconf's AC_TRY_LINK/AC_CHECK_LIB that checks the given given C++ program (3rd argument) successfully compiles.
# If compilation succeeds then:
# * the library (1st argument) is added to the $LIBS list (keeping this list unique).
# * the action-if-true argument is executed (the difference from AC_TRY_LINK/AC_CHECK_LIB is that action-if-true does not replace above step)
AC_DEFUN([AX_TRY_LINK], [
	dnl Below steps is a workaround for the limitation, that variables may not contain symbols
	dnl like "+" or "-" (and library names can). See AC_CHECK_LIB source comments for more information.
	AS_VAR_PUSHDEF([ac_Lib], [ac_cv_lib_$1])

	ax_link_dynamically=no

	AC_CACHE_CHECK(
		[m4_if(m4_index([$1], [ ]), [-1], [for -l$1], [for libs: $1])],
		[ac_Lib],
		[
			AS_UNSET([LIBS_LIST])

			AS_FOR([], [ax_var], [$1], [
				LIBS_LIST="-l${ax_var} ${LIBS_LIST}"
			])

			dnl Save the current state
			AC_LANG_SAVE
			AC_LANG_CPLUSPLUS
			ax_try_link_save_LIBS=${LIBS}
			
			LIBS="${LIBS_LIST} ${LIBS}"

			AC_TRY_LINK([$2], [$3], [AS_VAR_SET([ac_Lib], [yes])], [AS_VAR_SET([ac_Lib], [no])])

			dnl Restore the state to original regardless to the result
			LIBS=${ax_try_link_save_LIBS}

			dnl If the linking failed and "--enable-static-linking" was given, try to link against dynamic library (for the case when library is only available as .so:
			AS_VAR_IF([ac_Lib], [no], [
				dnl Try to disable static linking (if it was enabled) just for this very library. This trick can only work for GNU LD:
				AS_IF([test "${enable_static_linking+set}" == "set" -a "${ac_gnu_ld}" == "yes"], [
					ax_try_link_save_LIBS=${LIBS}

					dnl Note: "-Wl,-static" affects only libraries ("-l"), following this option, see http://stackoverflow.com/questions/809794/use-both-static-and-dynamically-linked-libraries-in-gcc
					LIBS="-Wl,-Bdynamic ${LIBS_LIST} -Wl,-static ${LIBS}"
					
					_AS_ECHO_LOG([LIBS="${LIBS}" for dynamic library presence check])

					AC_TRY_LINK([$2], [$3], [
						AS_VAR_SET([ac_Lib], [yes])
						ax_link_dynamically=yes
					])
	
					dnl Restore the state to original regardless to the result
					LIBS=${ax_try_link_save_LIBS}
				])
			])

			AC_LANG_RESTORE
		]
	)

	dnl If the variable is set, we define a constant and push library to $LIBS by default or execute [4], otherwise execute [5].
	AS_VAR_IF([ac_Lib], [yes],
		[
			$4
			AC_DEFINE_UNQUOTED(AS_TR_CPP(HAVE_LIB$1))

			AS_IF([test "${ax_link_dynamically}" == "yes"], [LIBS="-Wl,-static ${LIBS}"])

			dnl Do not prepend a library, if it is already in the list:
			AS_FOR([], [ax_var], [$1], [
				(echo "${LIBS}" | grep -q -- "-l${ax_var} ") || LIBS="-l${ax_var} ${LIBS}"
			])

			AS_IF([test "${ax_link_dynamically}" == "yes"], [LIBS="-Wl,-Bdynamic ${LIBS}"])

			_AS_ECHO_LOG([LIBS="${LIBS}" after linking check succeeded])
		],
		[$5]
	)
	AS_VAR_POPDEF([ac_Lib])
]) # AX_TRY_LINK

# SYNOPSIS
#
# AX_RESET_HEADERS_CACHE(headers ...)
#
# DESCRIPTION
#
# This macro invalidates the headers cache variables created by previous AC_CHECK_HEADER/AC_CHECK_HEADERS checks.
# Should be used only internally.
#
AC_DEFUN([AX_RESET_HEADERS_CACHE], [
	AS_FOR([], [ax_var], [$1], [
		dnl You can replace "ac_cv_header_" with any prefix from http://www.gnu.org/software/autoconf/manual/html_node/Cache-Variable-Index.html
		AS_VAR_PUSHDEF([ax_Var], [ac_cv_header_${ax_var}])
		AS_UNSET([ax_Var])
		AS_VAR_POPDEF([ax_Var])
	])
]) # AX_RESET_HEADERS_CACHE
