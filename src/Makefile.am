if FITSIO
  bin_PROGRAMS = snalys enclosep
  FITS = fitsio.f
  enclosep_SOURCES = enclosep.cpp
  INCLUDES = -I @fitsio_include@
else
  bin_PROGRAMS = snalys
  FITS = nofitsio.f
endif
snalys_SOURCES = main.f init.f analyze.f results.f precalcs.f \
	davidon.f gridsearch.f powell.f genetic.f qsort.f rtsafe.f date.f time.f \
	mlfunc.f mate.f integrate.f cosmo-routines.f \
	cosmo-struct.f gausslin.f christian.f array.f interpolate.f \
	cosmo.inc error.inc io.inc mate.inc operation.inc \
	probability.inc snalys.inc supernova.inc zbin.inc \
	gausslin.inc precalcs.inc $(FITS)
nodist_snalys_SOURCES = gausslin.inc
BUILT_SOURCES = gausslin.inc
MOSTLYCLEANFILES = gausslin.inc
snalys_FFLAGS = -fbounds-check $(AM_FFLAGS)
EXTRA_DIST = probability
dist-hook:
	rm -rf `find $(distdir)/probability -name CVS`

gausslin.inc : probability/sigmav1.out \
	probability/sigmav2.out probability/bv1.out \
	probability/bv2.out probability/m0v1.out \
	probability/m0v2.out probability/partof.pl
	cd probability/; @perlpath@ partof.pl
	mv probability/gausslin.inc .
