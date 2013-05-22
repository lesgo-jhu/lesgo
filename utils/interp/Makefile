##
##  Copyright 2012 Johns Hopkins University
##
##  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
##  use this file except in compliance with the License. You may obtain a copy of
##  the License at:
##
##    http://www.apache.org/licenses/LICENSE-2.0
##
##  Unless required by applicable law or agreed to in writing, software 
##  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
##  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
##  License for the specific language governing permissions and limitations under
##  the License.
##

# This is the main makefile.
# Use with GNU make.
# Relies on makedepf90 to do dependency lists:
# You can obtain this free program from:
#   http://www.helsinki.fi/~eedelman/makedepf90.html
# New $ARCH entries may be required, depending on your compiler.
# Also requires fpx3 fortran preprocessor available at
#   http://wwwuser.gwdg.de/~jbehren/fpx3.html

MK_INCL_PATH=../
include $(MK_INCL_PATH)/Makefile.in

EXE = interp

SRCS =  interp.f90

COMPSTR = '$(FPP) $$< > t.$$<; $$(FC) -c -o $$@ $$(FFLAGS) t.$$<'

include .depend_interp

.depend_interp: $(SRCS)
	mkdir -p $(OPATH) $(MPATH);
	makedepf90 -r $(COMPSTR) -b $(OPATH) -o $(EXE) $(SRCS) > .depend_interp

.PHONY : clean
clean:
	rm -rf $(OPATH)/* $(FOBJ) .depend* $(MPATH)/*.mod
	rm -f t.*
