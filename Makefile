##
##  Copyright (C) 2009-2013  Johns Hopkins University
##
##  This file is part of lesgo.
##
##  lesgo is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  lesgo is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
##

# This is the main makefile.
# Use with GNU make.
# Relies on makedepf90 to do dependency lists:
# You can obtain this free program from:
#   http://www.helsinki.fi/~eedelman/makedepf90.html
# New $ARCH entries may be required, depending on your compiler.
# Also requires fpx3 fortran preprocessor available at
#   http://wwwuser.gwdg.de/~jbehren/fpx3.html

MK_INCL_PATH=.
include Makefile.in

EXE := lesgo

SRCS =  cfl_util.f90 \
	clocks.f90 \
	convec.f90 \
        derivatives.f90 \
	debug_mod.f90 \
        divstress_uv.f90 \
	divstress_w.f90 \
	dns_stress.f90 \
        emul_complex.f90 \
        fft.f90 \
	finalize.f90 \
	forcing.f90 \
	fringe_util.f90 \
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
	open_file.f90 \
        padd.f90 \
	param.f90 \
	param_output.f90\
	press_stag_array.f90 \
        ran3.f90 rmsdiv.f90 \
        scaledep_dynamic.f90 \
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

RNS_LS_SRCS = rns_base_ls.f90 rns_ls.f90

RNS_CYL_SKEW_LS_SRCS = rns_cyl_skew_ls.f90

TURBINES_SRCS = turbines.f90 turbines_base.f90

#################################################### Tony ATM
ATM_SRCS = atm_base.f90 atm_input_util.f90                           \
           actuator_turbine_model.f90 atm_lesgo_interface.f90
#################################################### Tony ATM

CPS_SRCS = concurrent_precursor.f90

ifeq ($(USE_MPI), yes)
  SRCS += mpi_transpose_mod.f90 tridag_array_pipelined.f90 mpi_defs.f90
  EXE := $(EXE)-mpi
endif

ifeq ($(USE_CPS), yes)
  SRCS += $(CPS_SRCS)
  EXE := $(EXE)-cps
endif

ifeq ($(USE_STREAKS), yes)
  EXE := $(EXE)-streaks
endif

ifeq ($(USE_LVLSET), yes)
  SRCS += $(LVLSET_SRCS)
  EXE := $(EXE)-ls
endif

ifeq ($(USE_RNS_LS), yes)
  SRCS += $(RNS_LS_SRCS)
  ifeq ($(USE_CYL_SKEW_LS), yes)
    SRCS += $(RNS_CYL_SKEW_LS_SRCS)
  endif
  EXE := $(EXE)-rns
endif

ifeq ($(USE_CYL_SKEW_LS), yes)
  SRCS += $(CYL_SKEW_LS_SRCS)
  EXE := $(EXE)-cs
endif

ifeq ($(USE_TURBINES), yes)
  SRCS += $(TURBINES_SRCS)
  EXE := $(EXE)-turbines
endif

#################################################### Tony ATM
ifeq ($(USE_ATM), yes)
  SRCS += $(ATM_SRCS)
  EXE := $(EXE)-ATM
endif
#################################################### Tony ATM

ifeq ($(OUTPUT_EXTRA), yes)
  EXE := $(EXE)-exout
endif

ifeq ($(USE_DYN_TN), yes)
  EXE := $(EXE)-dyntn
endif

ifeq ($(USE_BINARY), yes)
  EXE := $(EXE)-binary
endif

ifeq ($(USE_SAFETYMODE), no)
  EXE := $(EXE)-safety_off
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
	rm -rf $(OPATH)/* $(FOBJ) .depend* $(MPATH)/*.mod
	rm -f t.*
