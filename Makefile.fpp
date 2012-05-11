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

ifeq ($(OUTPUT_EXTRA), yes)
  FPP += -DOUTPUT_EXTRA
endif

ifeq ($(USE_MPI), yes)
  FPP += -DMPI
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

ifeq ($(USE_CPS), yes)
  FPP += -DCPS
endif

ifeq ($(USE_BINARY), yes)
  FPP += -DBINARY
endif
