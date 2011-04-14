#
# This file contains checks for setting the FPP flags
#
ifeq ($(USE_DBLPREC), yes)
  FPP += -DDBLPREC
endif

ifeq ($(DEBUG), yes)
  FPP += -DDEBUG
endif

ifeq ($(VERBOSE), yes)
  FPP += -DVERBOSE
endif

ifeq ($(DEVEL), yes)
  FPP += -DDEVEL
endif

ifeq ($(USE_MPI), yes)
  FPP += -DMPI
  LIBS += -lmpichf90 -lfmpich -lmpich
endif

ifeq ($(USE_CFL_DT), yes)
  FPP += -DCFL_DT
endif

ifeq ($(USE_DYNALLOC),yes)
  FPP += -DDYNALLOC
endif

ifeq ($(WRITE_ENDIAN),LITTLE)
  FPP += -DWRITE_LITTLE_ENDIAN
endif

ifeq ($(WRITE_ENDIAN),BIG)
  FPP += -DWRITE_BIG_ENDIAN
endif

ifeq ($(READ_ENDIAN),LITTLE)
  FPP += -DREAD_LITTLE_ENDIAN
endif

ifeq ($(READ_ENDIAN),BIG)
  FPP += -DREAD_BIG_ENDIAN
endif


ifeq ($(USE_TREES_LS), yes)
  FPP += -DTREES_LS
endif

ifeq ($(USE_PC_SCHEME), DEFAULT)
  FPP += -DPC_SCHEME_1
endif

ifeq ($(USE_PC_SCHEME), PC_SCHEME_0)
  FPP += -DPC_SCHEME_0
endif

ifeq ($(USE_PC_SCHEME), PC_SCHEME_1)
  FPP += -DPC_SCHEME_1
endif 

ifeq ($(USE_PC_SCHEME), PC_SCHEME_2)
  FPP += -DPC_SCHEME_2
endif  

ifeq ($(USE_LVLSET), yes)
  FPP += -DLVLSET
endif

ifeq ($(USE_CYL_SKEW_LS), yes)
  FPP += -DCYL_SKEW_LS
endif

ifeq ($(USE_RNS_LS), yes)
  FPP += -DRNS_LS
endif

ifeq ($(USE_TURBINES), yes)
  FPP += -DTURBINES
endif

ifeq ($(USE_DYN_TN), yes)
  FPP += -DDYN_TN
endif
