##
##  Copyright (C) 2011-2013  Johns Hopkins University
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

ifeq ($(USE_SAFETYMODE), yes)
  FPP += -DSAFETYMODE
endif

ifeq ($(USE_FFTW3), yes)
  FPP += -DFFTW3
endif

ifeq ($(USE_CGNS), yes)
  FPP += -DCGNS
endif
