dnl Process this file with autoconf to produce a configure script.
AC_INIT([snalys], [0.2], [rahman@physto.se])
dnl Determine the type of architecture.
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE([dist-bzip2])
dnl AC_PACKAGE_TARNAME

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_F77
AC_LANG_FORTRAN77
AC_F77_LIBRARY_LDFLAGS
AC_PATH_PROG(rdocpath,robodoc,'')
AM_CONDITIONAL(ROBODOC, test -x $rdocpath)
AC_PATH_PROG(perlpath,perl,'')
AC_PATH_PROG(latexpath,latex,'')
AC_PATH_PROG(dvipspath,dvips,'')
AM_CONDITIONAL(MAKEPS, test -x $latexpath && test -x $dvipspath)
AC_PROG_INSTALL

AC_ARG_WITH(fitsiodir,
	[  --with-fitsiodir==VAL       Where to find the cfitsio libs],
	[fitsiodir="$withval"],[fitsiodir=""])
dnl AC_CHECK_LIB(cfitsio,ftvers, \
dnl    LIBS=LIBS:" -lcfitsio -lm "; cd src; ln -s fitsio.f fits.f;, \
dnl    cd src; ln -s nofitsio.f fits.f)
if test "$fitsiodir"; then
   AC_MSG_RESULT("Assuming that cfitsio is where you told me: $fitsiodir")
   LIBS="$LIBS -L$fitsiodir/lib -lcfitsio"
   fitsio_include="$fitsiodir/include"
   HAVEFITSIO=yes
   AC_SUBST(fitsio_include)
else
   AC_CHECK_LIB(cfitsio,ftvers,LIBS="$LIBS -lcfitsio";HAVEFITSIO=yes,HAVEFITSIO=no,-lm)
fi
AM_CONDITIONAL(FITSIO, test $HAVEFITSIO = yes)
AC_CHECK_LIB(pgplot,pgopen)

dnl Check for the cernlib
dnl AC_MSG_CHECKING(for the cernlib)
dnl AC_ARG_WITH(cernlib,
dnl   [  --with-cernlib=dir         give the directory to the cernlib],
dnl   cernlibdir="$withval", cernlibdir="/cern/pro")
dnl AC_MSG_RESULT($cernlibdir)
dnl AC_SUBST(cernlibdir)

dnl AC_CHECK_FILE($cernlibdir/libpacklib.a, PACKLIB_FOUND=yes,
dnl 	[ echo "*** Cannot find packlib!" ;
dnl 	  echo "*** Make sure that the path to cernlib is correct";
dnl 	  echo "*** Use --with-cernlib=dir, where dir is the path to cernlib.";
dnl 	  exit 1;] )

dnl Checks for libraries.
dnl AC_CHECK_FILE(lib/libpgplot-${target}.a, \
dnl   cp lib/libpgplot-${target}.a lib/libpgplot.a)

dnl Checks for header files.

dnl Checks for typedefs, structures, and compiler characteristics.

dnl Checks for library functions.

AC_OUTPUT([Makefile src/Makefile doc/Makefile \
		    doc/ascii/Makefile \
		    doc/html/Makefile])
dnl 		    doc/latex/Makefile])
