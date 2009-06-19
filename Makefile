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
ARCH = linux_intel
FCOMP = ifort
LIBPATH = -L/opt/fftw-2.1.5/lib -L/opt/mpich2-1.1-ifort/lib/
LIBS = $(LIBPATH) -lrfftw -lfftw -lm -lmpichf90 -lmpichf90 -lfmpich -lmpich

#--64-bit mode: may want to do export OBJECT_MODE=64
q64 = no

# watch the whitespace here
USE_MPI = no
USE_OPENMP = no
    #--not fully supported by all parts of the code
USE_DYNALLOC = no
    #--still experimental

USE_TREES_LS = no
USE_LVLSET = yes

FPP = fpx3
ifeq ($(USE_MPI), yes)
  FPP += -DMPI
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

# Directory for the .o files
OPATH = obj
# Directory for the .mod files, if your compiler generates them
# May want to just make this 'obj' as well
MPATH = mod

ifeq ($(FCOMP),ifort)
  FPP += -DIFORT
ifeq ($(USE_MPI), yes)
  FC = /opt/mpich2-1.1-ifort/bin/mpif90
else
  FC = ifort
endif

#  FFLAGS = -O0 -traceback -g -r8
  FFLAGS = -O0 -r8 -check all -g -traceback -debug all
#  FFLAGS = -fast
#  FFLAGS = -O3 -ipo
#  FFLAGS = -O3 -r8
#  FFLAGS = -O2 
#  FFLAGS = -axSSE4.2 -xS -ftz -ip -ipo -O3 
  FFLAGS += -warn all 
  #FDEBUG = -g -debug all
  FPROF = -p
  LDFLAGS = -threads
  ifeq ($(USE_MPI), yes)
    MODDIR = -I/opt/mpich2-1.1-ifort/include -I$(MPATH) -module /opt/mpich2-1.1-ifort/include -module $(MPATH)  
  else
    MODDIR = -I$(MPATH) -module $(MPATH)
  endif
  FFLAGS += $(MODDIR)
endif

ifeq ($(FCOMP),gfortran)
  FPP += -DGFORTRAN
  FC = gfortran
  FFLAGS = -O2
  FFLAGS += -Wall
  FDEBUG = -g
  FPROF = -p
  LDFLAGS = -static -nothreads
  MODDIR = -I$(MPATH) -J$(MPATH)  
  FFLAGS += $(MODDIR)  
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

ifeq ($(USE_MPI), yes)
  SRCS += mpi_transpose_mod.f90 tridag_array_pipelined.f90
endif

ifeq ($(USE_TREES_LS), yes)
  SRCS += $(TREES_LS_SRCS)
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

cylinder_skew: cylinder_skew.f90
	$(FC) -o $@ $(FFLAGS) -lgeometry $<

# Other support programs are listed below this point
interp: interp.f90
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
	rm -rf $(FOBJ) .depend* $(MPATH)/*.mod
	rm -f t.*
