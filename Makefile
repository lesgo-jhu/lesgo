# This is the main makefile.
# Use with GNU make.
# Relies on makedepf90 to do dependency lists:
# You can obtain this free program from:
#   http://www.helsinki.fi/~eedelman/makedepf90.html
# New $ARCH entries may be required, depending on your compiler.
# Also requires fpx3 fortran preprocessor available at
#   http://wwwuser.gwdg.de/~jbehren/fpx3.html

include Makefile.in

EXE = lesgo

SRCS =  \
	bottombc.f90 \
        cfl.f90 convec.f90 \
	ddx.f90 ddxy.f90 ddy.f90 ddz_uv.f90 ddz_w.f90 \
        dealias1.f90 dealias2.f90 debug_mod.f90 \
        divstress_uv.f90 divstress_w.f90 dns_stress.f90 \
        energy.f90 \
        fft.f90 filt_da.f90 forcing.f90 functions.f90 grid.f90 \
        ic.f90 ic_dns.f90 immersedbc.f90 initial.f90 \
	interpolag_Sdep.f90 interpolag_Ssim.f90 io.f90 \
        lagrange_Sdep.f90 lagrange_Ssim.f90 \
	main.f90 messages.f90 \
        padd.f90 param.f90 param_output.f90\
	press_stag_array.f90 \
        ran3.f90 rmsdiv.f90 \
        scaledep_dynamic.f90 scalars_module.f90 scalars_module2.f90 \
        sgs_stag.f90 sgsmodule.f90 sim_param.f90 stat_defs.f90 \
	strmod.f90 \
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

CYL_SKEW_LS_SRCS = cyl_skew_base_ls.f90 cyl_skew_ls.f90

RNS_LS_SRCS = rns_base_ls.f90 rns_ls.f90 rns_cyl_skew_ls.f90

TURBINES_SRCS = turbines.f90

TSUM_POST_DEPS = utils/tsum_post.f90 $(OPATH)/types.o $(OPATH)/param.o $(OPATH)/stat_defs.o $(OPATH)/grid.o
 
TSUM_POST_COMP2 = $(FPP) $< > t.tsum_post.f90; $(FC) -o $@ $(FFLAGS) $(LIBPATH) t.tsum_post.f90 \
	$(OPATH)/param.o $(OPATH)/stat_defs.o $(OPATH)/grid.o

ifeq ($(USE_MPI), yes)
  SRCS += mpi_transpose_mod.f90 tridag_array_pipelined.f90 mpi_defs.f90
  TSUM_POST_COMP1 = $(FPP) mpi_defs.f90 > t.mpi_defs.f90; $(FC) -c -o $(OPATH)/mpi_defs.o $(FFLAGS) t.mpi_defs.f90; 
  TSUM_POST_COMP2 += $(OPATH)/mpi_defs.o
endif

ifeq ($(USE_TREES_LS), yes)
  SRCS += $(TREES_LS_SRCS)
endif

ifeq ($(USE_LVLSET), yes)
  SRCS += $(LVLSET_SRCS)
endif

ifeq ($(USE_CYL_SKEW_LS), yes)
  SRCS += $(CYL_SKEW_LS_SRCS)
  TSUM_POST_DEPS += $(OPATH)/cyl_skew_base_ls.o
  TSUM_POST_COMP2 += $(OPATH)/cyl_skew_base_ls.o
endif

ifeq ($(USE_RNS_LS), yes)
	SRCS += $(RNS_LS_SRCS)
endif

ifeq ($(USE_TURBINES), yes)
  SRCS += $(TURBINES_SRCS)
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
tsum_post: $(TSUM_POST_DEPS)
	$(TSUM_POST_COMP1)
	$(TSUM_POST_COMP2)

cyl_skew_post_ls: utils/cyl_skew_post_ls.f90 $(OPATH)/types.o \
	$(OPATH)/param.o $(OPATH)/cyl_skew_base_ls.o 
	$(FPP) $< > t.cyl_skew_post_ls.f90; $(FC) -o $@ \
	$(CYLINDER_SKEW_PRE_LS_FFLAGS) $(LIBPATH) t.cyl_skew_post_ls.f90

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
	echo \0,0.,0.>./total_time.dat
	rm -rf $(OPATH)/* $(FOBJ) .depend* $(MPATH)/*.mod
	rm -f t.*
