# This is the main makefile.
# Use with GNU make.
# Relies on makedepf90 to do dependency lists:
# You can obtain this free program from:
#   http://www.helsinki.fi/~eedelman/makedepf90.html
# New $ARCH entries may be required, depending on your compiler.
# Also requires fpx3 fortran preprocessor available at
#   http://wwwuser.gwdg.de/~jbehren/fpx3.html

SHELL = /bin/bash
EXE = lesgo
FCOMP = ifort
LIBPATH = -L${HOME}/lib -L${HOME}/lib64 -L/opt/fftw-2.1.5/lib -L/usr/local/lib -L/usr/local/lib64
LIBS = $(LIBPATH) -lrfftw -lfftw -lm

#--64-bit mode: may want to do export OBJECT_MODE=64
q64 = no

#--Set global DEBUG flag;
#--Still have to set DEBUG in individual routines
DEBUG=no
#--Set global VERBOSE flag;
VERBOSE=no

# watch the whitespace here
USE_MPI = yes
USE_OPENMP = no

#--not fully supported by all parts of the code
USE_DYNALLOC = no
#--still experimental

USE_LVLSET = yes
USE_CYLINDER_SKEW_LS = yes
USE_RNS = yes

USE_TREES_LS = no

FPP = fpx3

ifeq ($(DEBUG), yes)
  FPP += -DDEBUG
endif

ifeq ($(VERBOSE), yes)
  FPP += -DVERBOSE
endif

ifeq ($(USE_MPI), yes)
  FPP += -DMPI
  LIBS += -lmpichf90 -lfmpich -lmpich
endif

ifeq ($(USE_DYNALLOC),yes)
  FPP += -DDYNALLOC
endif

ifeq ($(USE_TREES_LS), yes)
  FPP += -DTREES_LS
endif

ifeq ($(USE_LVLSET), yes)
  FPP += -DLVLSET
endif

ifeq ($(USE_CYLINDER_SKEW_LS), yes)
  FPP += -DCYLINDER_SKEW_LS
endif

ifeq ($(USE_RNS), yes)
  FPP += -DRNS
endif

# Directory for the .o files
OPATH = obj
# Directory for the .mod files, if your compiler generates them
# May want to just make this 'obj' as well
MPATH = mod

ifeq ($(FCOMP),ifort)
  
  FPP += -DIFORT
  
  ifeq ($(USE_MPI), yes)
    FC = mpif90
  else
    FC = ifort
  endif

  FFLAGS = -O0 -check bounds -g -debug all -traceback
  #FFLAGS = -fast
  #FFLAGS = -O3 -ipo
  #FFLAGS = -O3 -ip -ipo -ftz
  #FFLAGS = -axSSE4.2 -xS -ftz -ip -ipo -O3 
  FFLAGS += -warn all -mcmodel=medium
  #FDEBUG = -g -debug all
  FPROF = -p
  LDFLAGS = -threads -shared-intel
  MODDIR = -I$(MPATH) -module $(MPATH)
  FFLAGS += $(MODDIR)
  CYLINDER_SKEW_PRE_LS_FFLAGS = $(FFLAGS) -r8
endif

ifeq ($(FCOMP),gfortran)
  FPP += -DGFORTRAN
  ifeq ($(USE_MPI), yes)
    FC = mpif90
  else
    FC = gfortran
  endif
  FFLAGS = -O0 -fbounds-check
#  FFLAGS = -O2 -ffree-form -ffixed-line-length-none
  FFLAGS += -Wall
  FDEBUG = -g
  FPROF = -p
  LDFLAGS = -static 
  MODDIR = -I$(MPATH) -J$(MPATH)
  FFLAGS += $(MODDIR)  
  CYLINDER_SKEW_PRE_LS_FFLAGS += -fdefault-real-8 -fdefault-double-8
endif

ifeq ($(FCOMP),xlf)
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
  #find out details of how things are stored
  #FFLAGS += -qsource -qattr=full -qxref=full
  #ifeq ($(USE_OPENMP), yes)
    #FFLAGS += -qsmp=omp
  #endif
  #FDEBUG = -g
  #FPROF = -p
  ifeq ($(q64),yes)
    FFLAGS += -q64 -qarch
    LDFLAGS =
    LIBPATH =
    LIBS =
  else
    #LDFLAGS = -bmaxdata:0x80000000 -bmaxstack:0x10000000
    ## NOTE: you'll need to modify this!
    ## or specify on command line
    #LIBPATH = -L${HOME}/fftw/fftw2/lib
    #LIBS = $(LIBPATH) -lsrfftw -lsfftw -lm
    LIBS = $(LIBPATH) -lrfftw -lfftw -lm
  endif
  MODDIR = -I$(MPATH) -qmoddir=$(MPATH)  # where look for/put .mod files
  FFLAGS += $(MODDIR)
  CYLINDER_SKEW_PRE_LS_FFLAGS = $(FFLAGS) 
  #-qautodbl=dbl4 -qrealsize=8
endif

ifeq ($(FCOMP),g95)
  FPP += -DG95
  FC = g95
  FFLAGS = -O0
  FDEBUG = -g
  FPROF = -p
  LDFLAGS =
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
        fft.f90 filt_da.f90 forcing.f90 functions.f90 grid.f90 \
        ic.f90 ic_dns.f90 immersedbc.f90 initial.f90 \
	interpolag_Sdep.f90 interpolag_Ssim.f90 io.f90 \
        lagrange_Sdep.f90 lagrange_Ssim.f90 \
	main.f90 messages.f90 \
        padd.f90 param.f90 \
	press_stag_array.f90 \
        ran3.f90 rmsdiv.f90 \
        scaledep_dynamic.f90 scalars_module.f90 scalars_module2.f90 \
        sgs_stag.f90 sgsmodule.f90 sim_param.f90 stat_defs.f90 stats.f90 stats_init.f90 \
	std_dynamic.f90 string_util.f90 \
        test_filtermodule.f90 topbc.f90 \
        tridag.f90 tridag_array.f90 types.f90 \
        unpadd.f90 \
	wallstress.f90 wallstress_dns.f90

#--these also depend on linear_simple.f90
TREES_LS_SRCS = string_util.f90 \
		trees_ls.f90 trees_base_ls.f90 trees_fmodel_ls.f90 \
                trees_io_ls.f90 trees_post_mod_ls.f90 trees_setup_ls.f90 \
		trees_global_fmask_ls.f90

LVLSET_SRCS = level_set_base.f90 level_set.f90 linear_simple.f90

CYLINDER_SKEW_LS_SRCS = cylinder_skew_base_ls.f90 cylinder_skew_ls.f90

ifeq ($(USE_MPI), yes)
  SRCS += mpi_transpose_mod.f90 tridag_array_pipelined.f90
endif

ifeq ($(USE_TREES_LS), yes)
  SRCS += $(TREES_LS_SRCS)
endif

ifeq ($(USE_LVLSET), yes)
  SRCS += $(LVLSET_SRCS)
endif

ifeq ($(USE_CYLINDER_SKEW_LS), yes)
  SRCS += $(CYLINDER_SKEW_LS_SRCS)
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

tsum_post: utils/tsum_post.f90 $(OPATH)/types.o $(OPATH)/param.o $(OPATH)/stat_defs.o $(OPATH)/grid.o $(OPATH)/cylinder_skew_base_ls.o
	$(FPP) utils/mpi_defs.f90 > t.mpi_defs.f90; $(FC) -c -o $(OPATH)/mpi_defs.o $(FFLAGS) t.mpi_defs.f90
	$(FPP) $< > t.tsum_post.f90; $(FC) -o $@ $(FFLAGS) $(LIBPATH) t.tsum_post.f90 \
	$(OPATH)/param.o $(OPATH)/stat_defs.o $(OPATH)/mpi_defs.o $(OPATH)/grid.o $(OPATH)/cylinder_skew_base_ls.o

cylinder_skew_pre_ls: utils/cylinder_skew_pre_ls.f90 $(OPATH)/param.o $(OPATH)/cylinder_skew_base_ls.o
	$(FPP) utils/mpi_defs.f90 > t.mpi_defs.f90; $(FC) -c -o $(OPATH)/mpi_defs.o $(FFLAGS) t.mpi_defs.f90
	$(FPP) $< > t.cylinder_skew_pre_ls.f90; $(FC) -o $@ \
	$(CYLINDER_SKEW_PRE_LS_FFLAGS) $(LIBPATH) \
	-lgeometry t.cylinder_skew_pre_ls.f90 $(OPATH)/mpi_defs.o $(OPATH)/cylinder_skew_base_ls.o

cylinder_skew_post_ls: utils/cylinder_skew_post_ls.f90 $(OPATH)/types.o \
	$(OPATH)/param.o $(OPATH)/cylinder_skew_base_ls.o 
	$(FPP) $< > t.cylinder_skew_post_ls.f90; $(FC) -o $@ \
	$(CYLINDER_SKEW_PRE_LS_FFLAGS) $(LIBPATH) t.cylinder_skew_post_ls.f90

# Other support programs are listed below this point
interp: utils/interp.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

convert_endian:	utils/convert_endian.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

# This part is experimental
trees_pre:  trees_pre.f90 $(OPATH)/types.o \
            $(OPATH)/param.o $(OPATH)/messages.o \
            $(OPATH)/string_util.o \
	    $(OPATH)/trees_base.o $(OPATH)/trees_setup.o  \
	    $(OPATH)/trees_output.o
	$(FC) -o $@ $(FFLAGS) $^ $(LDFLAGS)

trees_pre_ls:  trees_pre_ls.f90 $(OPATH)/types.o $(OPATH)/param.o $(OPATH)/messages.o \
	$(OPATH)/string_util.o $(OPATH)/trees_base_ls.o $(OPATH)/trees_setup_ls.o \
	$(OPATH)/trees_io_ls.o $(OPATH)/trees_global_fmask_ls.o
	$(FC) -o $@ $(FFLAGS) $^ $(LDFLAGS)

trees_apri_ls:  trees_apri_ls.f90 $(OPATH)/types.o \
	$(OPATH)/param.o $(OPATH)/messages.o \
	$(OPATH)/string_util.o \
	$(OPATH)/trees_base_ls.o $(OPATH)/trees_setup_ls.o \
	$(OPATH)/trees_io_ls.o $(OPATH)/linear_simple.o \
	$(OPATH)/trees_aprioriCD_ls.o
	$(FC) -o $@ $(FFLAGS) $^ $(LDFLAGS)

trees_full_apri_ls:  trees_full_apri_ls.f90 $(OPATH)/types.o \
	        $(OPATH)/param.o $(OPATH)/sim_param_vel.o \
                $(OPATH)/messages.o $(OPATH)/string_util.o \
                $(OPATH)/trees_base_ls.o $(OPATH)/trees_setup_ls.o \
                $(OPATH)/trees_io_ls.o $(OPATH)/linear_simple.o \
                $(OPATH)/trees_fmodel_ls.o $(OPATH)/trees_ls.o \
		$(OPATH)/level_set_base.o $(OPATH)/trees_global_fmask_ls.o
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $^

trees_post_ls:  trees_post_ls.f90 $(OPATH)/types.o \
                $(OPATH)/param.o $(OPATH)/messages.o \
                $(OPATH)/string_util.o \
                $(OPATH)/trees_base_ls.o $(OPATH)/trees_setup_ls.o \
                $(OPATH)/trees_io_ls.o $(OPATH)/linear_simple.o \
                $(OPATH)/trees_post_mod_ls.o
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $^

merge_phi:  merge_phi.f90 $(OPATH)/types.o \
            $(OPATH)/param.o $(OPATH)/messages.o \
            $(OPATH)/string_util.o \
            $(OPATH)/trees_base_ls.o
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $^

# This doesn't remove .mod files--should be OK as long a dependency list 
# for the .o files is correct.
# FOBJ is defined in .depend
.PHONY : clean
clean :
	echo \0>./total_time.dat
	rm -rf $(OPATH)/* $(FOBJ) .depend* $(MPATH)/*.mod
	rm -f t.*
