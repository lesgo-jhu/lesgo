# This is the main makefile.
# Use with GNU make.
# Relies on makedepf90 to do dependency lists:
# You can obtain this free program from:
#   http://www.helsinki.fi/~eedelman/makedepf90.html
# New $ARCH entries may be required, depending on your compiler.
# Also requires fpx3 fortran preprocessor available at
#   http://wwwuser.gwdg.de/~jbehren/fpx3.html

SHELL = /bin/sh

EXE = LES_code

ARCH = linux_intel_71
#ARCH = linux_g95
ARCH = linux_intel_81
#ARCH = aix
#ARCH = osx_g95

# watch the whitespace here
USE_MPI = no
USE_TREES = no
USE_LVLSET = no

FPP = fpx3
ifeq ($(USE_MPI), yes)
  FPP += -DMPI
endif

ifeq ($(USE_TREES), yes)
  FPP += -DTREES
endif

ifeq ($(USE_LVLSET), yes)
  FPP += -DLVLSET
endif

# Directory for the .o files
OPATH = obj
# Directory for the .mod files, if your compiler generates them
# May want to just make this 'obj' as well
MPATH = mod

ifeq ($(ARCH),linux_intel_81)
  FPP += -DIFORT
  #FC = ifc
  FC = ifort
  #FFLAGS = -O0 -traceback -g
  FFLAGS = -O0 -check bounds
  #FFLAGS = -O3
  #FFLAGS = -fast
  #FFLAGS = -O3 -ipo
  #FFLAGS = -O3 -mp
  #FFLAGS = -g 
  FFLAGS += -warn all 
  FDEBUG = -g
  FPROF = -p
  LDFLAGS = -nothreads
  #LDFLAGS = -static -nothreads
  LIBPATH = -L/opt/fftw2/lib
  #LIBPATH = -L/home/chester/fftw/fftw2/lib
  #LIBS = $(LIBPATH) -lintel81_srfftw -lintel81_sfftw
  #LIBS = $(LIBPATH) -lintel81_drfftw -lintel81_dfftw
  LIBS = $(LIBPATH) -lrfftw -lfftw -lm
  MODDIR = -I$(MPATH) -module $(MPATH)  # where look for/put .mod files
  FFLAGS += $(MODDIR)
endif
ifeq ($(ARCH),linux_intel_71)
  FPP += -DIFC
  ifeq ($(USE_MPI), yes)
    FC = mpif90
  else
    FC = gfortran
  endif
  #FFLAGS = -O0
  FFLAGS = -O0 -stack_temps
  #FFLAGS = -O2
  #FFLAGS = -O3 -stack_temps
  #FFLAGS = -O3 -ipo
  #FFLAGS = -O3 -mp
  #FFLAGS = -zero -CB -CU -CS
  #FFLAGS = -CA -CB -CU -CS
  #FFLAGS = -O0
  FDEBUG = -g
  FPROF = -p
  ifeq ($(USE_MPI), yes)
    # problems with static linking & MPI
    LDFLAGS =
  else
    LDFLAGS = -static -pthread
    #LDFLAGS = -static  #--this doesnt always work...glibc/pthreads?
    #LDFLAGS = -static-libcxa  #--only intel library static
  endif
  LIBPATH =
  #LIBPATH = -L/home/chester/fftw/fftw2/lib
  #LIBS = $(LIBPATH) -lintel71_srfftw -lintel71_sfftw
  #LIBS = $(LIBPATH) -lintel71_drfftw -lintel71_dfftw
  LIBS = $(LIBPATH) -lrfftw -lfftw -lm
  #LIBS = $(LIBPATH) -lintel_drfftw -lintel_dfftw -lm
  #MODDIR = -I$(MPATH)   # where look for/put .mod files
  MODDIR = -I$(MPATH) -module $(MPATH)  # where look for/put .mod files
  FFLAGS += $(MODDIR)
endif
ifeq ($(ARCH),aix)
  FPP += -DXLF
  ifeq ($(USE_MPI), yes)
    FC = mpxlf95_r
  else
    FC = xlf95_r
  endif
  #FFLAGS = -qstrict -qsuffix=f=f90 -qsmp -O3 -qreport=smplist
  FFLAGS = -qstrict -qsuffix=f=f90 -O3
  #FFLAGS = -qstrict -qsuffix=f=f90 -O3 -qsmp=omp
  #FFLAGS = -qstrict -qsuffix=f=f90 -O0
  # find out details of how things are stored
  FFLAGS += -qsource -qattr=full -qxref=full -qmaxmem=-1
  FDEBUG = -g
  FPROF = -p
  LDFLAGS = -bmaxdata:0x80000000 -bmaxstack:0x10000000 -bnoquiet
  # NOTE: you'll need to modify this!
  #       or specify on command line
#  LIBPATH = -L/home/bluesky/schester/fftw/fftw2/lib
LIBPATH = -L/usr/local/lib32/r4i4 -L/home/bluesky/vijayant/fftw/lib
  #LIBS = $(LIBPATH) -lsrfftw -lsfftw -lm
  LIBS = $(LIBPATH) -lrfftw -lfftw -lm
#  LIBS = $(LIBPATH) -ldrfftw -ldfftw -lm
  MODDIR = -I$(MPATH) -qmoddir=$(MPATH)  # where look for/put .mod files
  FFLAGS += $(MODDIR)
endif
ifeq ($(ARCH),linux_g95)
  FPP += -DG95
  FC = g95
  FFLAGS = -O0
  FDEBUG = -g
  FPROF = -p
  LDFLAGS =
  LIBPATH = -L/home/chester/fftw/fftw2/lib
  LIBS = $(LIBPATH) -lsrfftw -lsfftw
  MODDIR =
  FFLAGS += $(MODDIR)
endif
ifeq ($(ARCH),osx_g95)
  FPP += -DG95
  FC = g95
  FFLAGS = -O0
  FDEBUG = -g
  FPROF = -p
  LDFLAGS =
  LIBPATH = -L/Users/stu/Work/fftw/fftw2/lib
  #LIBS = $(LIBPATH) -lsrfftw -lsfftw
  LIBS = $(LIBPATH) -ldrfftw -ldfftw
  MODDIR =
  FFLAGS += $(MODDIR)
endif

SRCS =  \
	bottombc.f90 \
        convec.f90 \
	ddx.f90 ddxy.f90 ddy.f90 ddz_uv.f90 ddz_w.f90 \
        dealias1.f90 dealias2.f90 debug_mod.f90 \
        divstress_uv.f90 divstress_w.f90 dns_stress.f90 \
        energy.f90 \
        fft.f90 filt_da.f90 forcing.f90 \
        ic.f90 ic_dns.f90 immersedbc.f90 initial.f90 \
	interpolag_Sdep.f90 interpolag_Ssim.f90 io.f90 \
	lagrange_Sdep.f90 lagrange_Ssim.f90 \
	main.f90 messages.f90 \
        padd.f90 param.f90 press_stag_array.f90 \
        ran3.f90 rmsdiv.f90 \
        scaledep_dynamic.f90 scalars_module.f90 scalars_module2.f90 \
        sgs_stag.f90 sgsmodule.f90 sim_param.f90 std_dynamic.f90 \
	string_util.f90 \
        test_filtermodule.f90 topbc.f90 \
        tridag.f90 tridag_array.f90 types.f90 \
        unpadd.f90 \
	wallstress.f90 wallstress_dns.f90
#        interpolag_Ssim_VIJ.f90 lagrange_Ssim_VIJ.f90 \
#	output_slice.f90

TREE_SRCS = linear_simple.f90 trees.f90 trees_base.f90 trees_BC.f90 \
            trees_Cd.f90 trees_CV.f90 trees_output.f90

LVLSET_SRCS = level_set.f90

ifeq ($(USE_MPI), yes)
#  SRCS += mpi_transpose_mod.f90
  SRCS += mpi_transpose_mod.f90 tridag_array_pipelined.f90
endif
ifeq ($(USE_TREES), yes)
  SRCS += $(TREE_SRCS)
endif
ifeq ($(USE_LVLSET), yes)
  SRCS += $(LVLSET_SRCS)
endif

#COMPSTR = '$(FPP) $$< > t.$$<; $$(FC) -c -o $$@ $$(FFLAGS) t.$$<; rm -f t.$$<'
COMPSTR = '$(FPP) $$< > t.$$<; $$(FC) -c -o $$@ $$(FFLAGS) t.$$<'

include .depend

.depend: $(SRCS)
	mkdir -p $(OPATH) $(MPATH);
	makedepf90 -r $(COMPSTR) -b $(OPATH) -o $(EXE) $(SRCS) > .depend
#	makedepf90 -r 'fpx3 $$< > t.$$<;$$(FC) -c -o $$(@) $$(FFLAGS) t.$$<;rm -f t.$$<'\
#                   -b $(OPATH) -o $(EXE) $(SRCS) > .depend

#.depend: $(SRCS)
#	mkdir -p $(OPATH) $(MPATH);
#	makedepf90 -r '$$(FC) -c -o $$(@) $$(FFLAGS) $$<'\
#                   -b $(OPATH) -o $(EXE) $(SRCS) > .depend

debug:
	$(MAKE) $(EXE) "FFLAGS = $(FDEBUG) $(FFLAGS)"

prof:
	$(MAKE) $(EXE) "FFLAGS = $(FPROF) $(FFLAGS)"

# Other support programs are listed below this point
#interp: interp.f90
interp: interp_scalars.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<
#interp: interp_scalars.f90
#	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

# This part is experimental
trees_pre:  trees_pre.f90 $(OPATH)/types.o $(OPATH)/param.o \
            $(OPATH)/sim_param.o $(OPATH)/trees_base.o $(OPATH)/trees_CV.o  \
	    $(OPATH)/trees_Cd.o $(OPATH)/trees_output.o $(OPATH)/trees_BC.o \
	    $(OPATH)/messages.o $(OPATH)/trees.o $(OPATH)/linear_simple.o \
	    $(OPATH)/string_util.o $(OPATH)/immersedbc.o $(OPATH)/sgsmodule.o
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $^

# This doesn't remove .mod files--should be OK as long a dependency list 
# for the .o files is correct.
# FOBJ is defined in .depend
.PHONY : clean
clean :
	echo \0>./total_time.dat
	rm -rf $(FOBJ) .depend* $(MPATH)/*.mod
	rm -rf ./output/*
	mkdir ./output/fields_3d
	mkdir output/code_used
	cp $(SRCS) param.f90 Makefile ./output/code_used/
