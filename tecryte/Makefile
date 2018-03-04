# Environment 
# compiler
#FC = gfortran
#FC = ifort
#FC = pgf90
FC = g95
#FC = xlf

LIBNAME=libtecryte
HEADER = tecryte.h

PREFIX=/home/eramudu
LIBDIR=$(PREFIX)/lib
INCLDIR=$(PREFIX)/include

ifeq ($(FC),gfortran)
OUT = $(LIBNAME).a
FFLAGS = -O3 -Wall
LIBS =
LDFLAGS =
endif

ifeq ($(FC),ifort)
OUT = $(LIBNAME)-intel.a
FFLAGS = -O3 -warn all
LIBS =
LDFLAGS =
endif

ifeq ($(FC),pgf90)
OUT = $(LIBNAME)-pgi.a
FFLAGS = -O3 -Minform=warn
LIBS = 
LDFLAGS =
endif

ifeq ($(FC),g95)
OUT = $(LIBNAME)-g95.a
FFLAGS = -O3 -Wall 
LIBS =
LDFLAGS =
endif

ifeq ($(FC),xlf)
FC=xlf95_r
OUT = $(LIBNAME)-xlf.a
FFLAGS = -O3
LIBS = 
LDFLAGS =
endif

RM = rm -f
CP = cp
AR = ar rcs
FTPTOF90 = ./ftptof90.sh

SRC = tecryte.f90 \
open_file.f90 \
write_real_data.f90 \
write_real_data_1d.f90 \
write_real_data_2d.f90 \
write_real_data_3d.f90 \
write_tecplot_header_xyline.f90 \
write_tecplot_header_nd.f90 

OBJ = $(SRC:.f90=.o)
MOD = $(SRC:.f90=.mod)

%.o : %.mod

.SUFFIXES: .f90

default: $(OUT)

.f90.o:
	$(FC) -fno-second-underscore $(INCLUDES) $(FFLAGS) -c $< -o $@

$(OUT): $(OBJ)
	$(AR) $(OUT) $(OBJ)

source:
	$(FTPTOF90)

clean:
	$(RM) $(OBJ) $(MOD) $(OUT)

install: $(OUT)
	$(CP) $(OUT) $(LIBDIR)
	$(CP) $(HEADER) $(INCLDIR)
uninstall:
	$(RM) $(LIBDIR)/$(OUT) 
	$(RM) $(INCLDIR)/$(HEADER)
