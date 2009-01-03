# This is the qpost makefile.
# Use with GNU make.
# Relies on makedepf90 to do dependency lists:
# You can obtain this free program from:
#   http://www.helsinki.fi/~eedelman/makedepf90.html
# New $ARCH entries may be required, depending on your compiler.
# Also requires fpx3 fortran preprocessor available at
#   http://wwwuser.gwdg.de/~jbehren/fpx3.html

SHELL = /bin/sh

EXE = qpost

#ARCH = linux_intel_71
#ARCH = linux_g95
ARCH = linux_intel_81
#ARCH = aix
#ARCH = osx_g95

#--this option only applies for 64-bit mode
#--when using this, it is a good idea to set the evironment variable
#  export OBJECT_MODE=64 in the loadlever script (and perhaps also before
#  "make"ing with this Makefile
q64 = yes

ifeq ($(ARCH),linux_intel_81)
  FC = ifort
  #FFLAGS = -O0 -traceback -g
  #FFLAGS = -O2
  FFLAGS = -O3
  #FFLAGS = -fast
  #FFLAGS = -O3 -ipo
  #FFLAGS = -O3 -mp
  #FFLAGS = -g 
  FFLAGS += -warn all 
  LDFLAGS = -nothreads
  #LDFLAGS = -static -nothreads
  LIBPATH = -L${HOME}/fftw/fftw2/lib
  #LIBS = $(LIBPATH) -lintel81_srfftw -lintel81_sfftw
  LIBS = $(LIBPATH) -lintel81_drfftw -lintel81_dfftw
endif
ifeq ($(ARCH),linux_intel_71)
  FC = ifc
  #FFLAGS = -O0 -stack_temps
  #FFLAGS = -O2
  FFLAGS = -O3 -stack_temps
  #FFLAGS = -O3 -ipo
  #FFLAGS = -O3 -mp
  #FFLAGS = -zero -CB -CU -CS
  #FFLAGS = -CA -CB -CU -CS
  #FFLAGS = -O0
  ifeq ($(USE_MPI), yes)
    # problems with static linking & MPI
    LDFLAGS =
  else
    LDFLAGS = -static -pthread
    #LDFLAGS = -static  #--this doesnt always work...glibc/pthreads?
    #LDFLAGS = -static-libcxa  #--only intel library static
  endif
  LIBPATH = -L${HOME}/fftw/fftw2/lib
  #LIBS = $(LIBPATH) -lintel71_srfftw -lintel71_sfftw
  LIBS = $(LIBPATH) -lintel71_drfftw -lintel71_dfftw
  #LIBS = $(LIBPATH) -lintel_rfftw -lintel_fftw -lm
  #LIBS = $(LIBPATH) -lintel_drfftw -lintel_dfftw -lm
endif
ifeq ($(ARCH),aix)
  FC = xlf95_r
  #FFLAGS = -qstrict -qsuffix=f=f90 -qsmp -O3 -qreport=smplist
  FFLAGS = -qstrict -qsuffix=f=f90 -O3
  #FFLAGS = -qstrict -qsuffix=f=f90 -O3 -qsmp=omp
  #FFLAGS = -qstrict -qsuffix=f=f90 -O0
  # find out details of how things are stored
  FFLAGS += -qsource -qattr=full -qxref=full
  ifeq ($(q64),yes)
    FFLAGS += -q64 -qarch
    LDFLAGS = 
    LIBPATH = 
    LIBS = 
  else
    LDFLAGS = -bmaxdata:0x80000000 -bmaxstack:0x10000000
    # NOTE: you'll need to modify this!
    #       or specify on command line
    LIBPATH = -L${HOME}/fftw/fftw2/lib
    #LIBS = $(LIBPATH) -lsrfftw -lsfftw -lm
    LIBS = $(LIBPATH) -ldrfftw -ldfftw -lm
  endif
endif
ifeq ($(ARCH),linux_g95)
  FC = g95
  FFLAGS = -O0
  LDFLAGS =
  LIBPATH = -L${HOME}/fftw/fftw2/lib
  LIBS = $(LIBPATH) -lsrfftw -lsfftw
endif
ifeq ($(ARCH),osx_g95)
  FC = g95
  FFLAGS = -O0
  LDFLAGS =
  LIBPATH = -L${HOME}/Work/fftw/fftw2/lib
  #LIBS = $(LIBPATH) -lsrfftw -lsfftw
  LIBS = $(LIBPATH) -ldrfftw -ldfftw
endif

SRCS =  qpost.f90 types.f90

#COMPSTR = '$(FPP) $$< > t.$$<; $$(FC) -c -o $$@ $$(FFLAGS) t.$$<; rm -f t.$$<'
COMPSTR = '$$(FC) -c -o $$@ $$(FFLAGS) $$<'

include .depend

.depend: $(SRCS)
	makedepf90 -r $(COMPSTR) -o $(EXE) $(SRCS) > .depend

#experimental
qpost_lite: qpost_lite.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

qpost_slice: qpost_slice.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

#--almost the same as qpost_f, oops
post_f: post_f.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

qpost_f: qpost_f.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

phi_merge: phi_merge.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

#--the _c versions of these routines are not working on all arch's yet
slice2ppm: types.f90 writeppm_i.f90 colormap_i.f90 slice2ppm.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) types.f90 writeppm_i.f90 \
	  colormap_i.f90 slice2ppm.f90

mean-CD: mean-CD.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<
	
qpost_planes: qpost_planes.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

# This doesn't remove .mod files--should be OK as long a dependency list 
# for the .o files is correct.
# FOBJ is defined in .depend
.PHONY : clean
clean :
	rm -rf $(FOBJ) .depend* *.mod
