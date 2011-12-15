# This is the main makefile.
# Use with GNU make.
# Relies on makedepf90 to do dependency lists:
# You can obtain this free program from:
#   http://www.helsinki.fi/~eedelman/makedepf90.html
# New $ARCH entries may be required, depending on your compiler.
# Also requires fpx3 fortran preprocessor available at
#   http://wwwuser.gwdg.de/~jbehren/fpx3.html

include Makefile.in

EXE := lesgo

SRCS =  cfl_util.f90 \
	clocks.f90 \
	convec.f90 \
        derivatives.f90 \
        dealias1.f90 \
	dealias2.f90 \
	debug_mod.f90 \
        divstress_uv.f90 \
	divstress_w.f90 \
	dns_stress.f90 \
        emul_complex.f90 \
        energy.f90 \
        fft.f90 \
	forcing.f90 \
	functions.f90 \
	grid.f90 \
        ic.f90 \
	ic_dns.f90 \
	initial.f90 \
	initialize.f90 \
	input_util.f90 \
	interpolag_Sdep.f90 \
	interpolag_Ssim.f90 \
	io.f90 \
        lagrange_Sdep.f90 \
	lagrange_Ssim.f90 \
	main.f90 \
	messages.f90 \
        padd.f90 \
	param.f90 \
	param_output.f90\
	press_stag_array.f90 \
        ran3.f90 rmsdiv.f90 \
        scaledep_dynamic.f90 \
	sgs_hist.f90 \
	sgs_param.f90 \
        sgs_stag_util.f90 \
	sim_param.f90 \
	stat_defs.f90 \
	std_dynamic.f90 \
	string_util.f90 \
        test_filtermodule.f90 \
        tridag.f90 \
	tridag_array.f90 \
	types.f90 \
        unpadd.f90 \
	wallstress.f90 \
	wallstress_dns.f90

LVLSET_SRCS = level_set_base.f90 level_set.f90 linear_simple.f90

CYL_SKEW_LS_SRCS = cyl_skew_base_ls.f90 cyl_skew_ls.f90

RNS_LS_SRCS = rns_base_ls.f90 rns_ls.f90 rns_cyl_skew_ls.f90

TURBINES_SRCS = turbines.f90

CPS_SRCS = concurrent_precursor.f90

ifeq ($(USE_MPI), yes)
  SRCS += mpi_transpose_mod.f90 tridag_array_pipelined.f90 mpi_defs.f90
  EXE := $(EXE)-mpi
endif

ifeq ($(USE_LVLSET), yes)
  SRCS += $(LVLSET_SRCS)
  EXE := $(EXE)-ls
endif

ifeq ($(USE_RNS_LS), yes)
  SRCS += $(RNS_LS_SRCS)
  EXE := $(EXE)-rns
endif

ifeq ($(USE_CYL_SKEW_LS), yes)
  SRCS += $(CYL_SKEW_LS_SRCS)
  EXE := $(EXE)-cs
endif

ifeq ($(USE_TURBINES), yes)
  SRCS += $(TURBINES_SRCS)
endif

ifeq ($(USE_CPS), yes)
  SRCS += $(CPS_SRCS)
  EXE := $(EXE)-cps
endif

ifeq ($(OUTPUT_EXTRA), yes)
  EXE := $(EXE)-exout
endif

ifeq ($(USE_DYN_TN), yes)
  EXE := $(EXE)-dyntn
endif

#COMPSTR = '$(FPP) $$< > t.$$<; $$(FC) -c -o $$@ $$(FFLAGS) t.$$<; rm -f t.$$<'
COMPSTR = '$(FPP) $$< > t.$$<; $$(FC) -c -o $$@ $$(FFLAGS) t.$$<'


include .depend

.depend: $(SRCS)
	mkdir -p $(OPATH) $(MPATH);
	makedepf90 -r $(COMPSTR) -b $(OPATH) -o $(EXE) $(SRCS) > .depend

debug:
	$(MAKE) $(EXE) "FFLAGS = $(FDEBUG) $(FFLAGS)"

prof:
	$(MAKE) $(EXE) "FFLAGS = $(FPROF) $(FFLAGS)"

# Other support programs are listed below this point
convert_endian:	utils/convert_endian.f90
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $<

# This doesn't remove .mod files--should be OK as long a dependency list 
# for the .o files is correct.
# FOBJ is defined in .depend
.PHONY : clean
clean :
	echo \0,0.,0.,0.,0.>./total_time.dat
	rm -rf $(OPATH)/* $(FOBJ) .depend* $(MPATH)/*.mod
	rm -f t.*
