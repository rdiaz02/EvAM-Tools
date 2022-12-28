dnl Initialise Autoconf
AC_PREREQ([2.69])
AC_INIT(
	[CT-CBN],
	[0.1.04b],
	[niko.beerenwinkel@bsse.ethz.ch],
	[ct-cbn],
	[https://www.bsse.ethz.ch/cbg/software/ct-cbn.html])
AC_CONFIG_SRCDIR([src/ct-cbn.c])



dnl =============================================================
dnl Check whether we want to set defaults for CFLAGS and CPPFLAGS
dnl =============================================================
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug], [Modifies various CFLAGS and CPPFLAGS required for debugging])]
)

AS_IF(
	[test "x$enable_debug" = "xyes"],
		[
			enable_cflags_setting=yes
			cflags_test="-O0 -g"
			
			enable_cppflags_setting=yes
			cppflags_test=""
		],
		[
			cflags_test="-O2"
			cppflags_test="-DNDEBUG"
		
			AC_MSG_CHECKING([whether configure should try to set CFLAGS])
			AS_IF(
				[test "x${CFLAGS+set}" = "xset"],
					[enable_cflags_setting=no],
					[enable_cflags_setting=yes]
			)
			AC_MSG_RESULT([$enable_cflags_setting])
                	
			AC_MSG_CHECKING([whether configure should try to set CPPFLAGS])
			AS_IF(
				[test "x${CPPFLAGS+set}" = "xset"],
					[enable_cppflags_setting=no],
					[enable_cppflags_setting=yes]
			)
			AC_MSG_RESULT([$enable_cppflags_setting])
		]
)



dnl =======================================
dnl Check for standard headers and programs
dnl =======================================
AX_COMPILER_VENDOR
AX_COMPILER_VERSION
AC_PROG_SED



dnl =======================
dnl Set CFLAGS and CPPFLAGS
dnl =======================
AS_IF(
	[test "x${enable_cflags_setting}" = "xyes"],
		[
		CFLAGS=""
		AX_APPEND_COMPILE_FLAGS([$cflags_test -Wall], [CFLAGS])
		CFLAGS=$( echo $CFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
		]
)

AS_IF(
	[test "x${enable_cppflags_setting}" = "xyes"],
		[
		CPPFLAGS=""
		AX_APPEND_COMPILE_FLAGS([$cppflags_test], [CPPFLAGS])
		CPPFLAGS=$( echo $CPPFLAGS | $SED -e 's/^ *//' -e 's/ *$//' )
		]
)



dnl ===================
dnl Initialise Automake
dnl ===================
AM_INIT_AUTOMAKE([1.15 foreign dist-bzip2 no-dist-gzip subdir-objects silent-rules])



dnl ===================
dnl Check for libraries
dnl ===================
AC_OPENMP
PKG_CHECK_MODULES([GSL], [gsl])



dnl ========
dnl Finalise
dnl ========
AC_CONFIG_FILES([Makefile])
AC_OUTPUT



dnl ============================
dnl Report configuration to user
dnl ============================
AC_MSG_RESULT([
	$PACKAGE_NAME $VERSION
	Using $ax_cv_c_compiler_vendor, $ax_cv_c_compiler_version

	CC:       $CC
	CFLAGS:   $CFLAGS
	CPPFLAGS: $CPPFLAGS
	LDFLAGS:  $LDFLAGS
])
