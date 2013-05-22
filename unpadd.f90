!!
!!  Copyright 2009,2010,2011,2012 Johns Hopkins University
!!
!!  Licensed under the Apache License, Version 2.0 (the "License"); you may not 
!!  use this file except in compliance with the License. You may obtain a copy of
!!  the License at:
!!
!!    http://www.apache.org/licenses/LICENSE-2.0
!!
!!  Unless required by applicable law or agreed to in writing, software 
!!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!!  License for the specific language governing permissions and limitations under
!!  the License.
!!

!*******************************************************************************
subroutine unpadd(cc,cc_big)
!*******************************************************************************
use types,only:rprec
use param,only:ld,nx,ny,ny2,nz,ld_big,lh
implicit none

!  cc and cc_big are interleaved as complex arrays
real(rprec), dimension( ld, ny ) :: cc
real(rprec), dimension( ld_big, ny2 ) :: cc_big

integer :: ny_h, j_s, j_big_s

ny_h = ny/2

cc(:nx,:ny_h) = cc_big(:nx,:ny_h)

! oddballs
cc(ld-1:ld,:) = 0._rprec
cc(:,ny_h+1) = 0._rprec

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2
cc(:nx,j_s:ny) = cc_big(:nx,j_big_s:ny2)

end subroutine unpadd
