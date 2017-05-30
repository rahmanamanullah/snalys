A cosmology fitter from magnitudes and redshifts that work with arbitrary
magnitude distributions.  This can be used to account for the effects
introduced by gravitational lensing.

The easiest way to install the software is to enter directory snoc/snalys
and type:

  ./configure
  make

This should produce a binary for your architecture. Currently i386 linux and
alpha OSF are supported but there is no reason why it will not work on other
platforms. The only thing you have to do is to compile the pgplot libraries
(snoc/lib/pgplot*.tar) for that specific architecture. To clean up the
object files type

  make clean

Before you want to compile the software for a different architecture type

  make distclean

and the type ./configure again. Good luck!

To run the program you also have to specify the environment variable
PGPLOT_DIR to point at the snoc/snalys/lib directory.

Typing

  snalys -h

is the closest thing to a documentation today, but I am working on it.


                            Rahman Amanullah 2000-09-21
