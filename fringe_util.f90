!!
!!  Copyright 2011,2012 Johns Hopkins University
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

!***********************************************************************
module fringe_util
!***********************************************************************
implicit none

save
private

public fringe_init, &
       fringe_weighting

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine fringe_init( istart, iplateau, iend )
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Sets the beginning, ending and plateau index within the fringe
! region. Provides a common routine to do this.
!
use types, only : rprec
use param, only : nx, fringe_region_end, fringe_region_len
implicit none

integer, intent(out) :: istart, iplateau, iend

iend = floor (fringe_region_end * nx + 1.0_rprec)
iplateau = floor (( fringe_region_end - fringe_region_len / 4 ) * nx + 1.0_rprec)
istart = floor ((fringe_region_end - fringe_region_len) * nx + 1.0_rprec)

return
end subroutine fringe_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function fringe_weighting( i, istart, iend ) result( w ) 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Provides the weigting for the sampled (desired) velocity field in the
! fringe region. The result from this function (w) is used to compute
! the fringe region velocity as
!
!   u_fringe = w * u_sample + (1-w)*u
!
! where u_fringe is the fringe velocity, u_sample is the sampled
! velocity and u is the unmodulated velocity in the fringe region right
! after the projection step.
!
use types, only : rprec
use param, only : pi
implicit none

integer, intent(in) :: i, istart, iend
real(rprec) :: w

! Linear profile
!beta = real ( i - istart, rprec ) / real ( iend - istart, rprec )
! Sine profile
!beta = 0.5_rprec * ( 1._rprec - cos (pi * real (i - istart, rprec)  &
!                                       / (iend - istart)) )
! Sine profile with plateau
if ( i > iend ) then 
   w = 1.0_rprec
else
   w = 0.5_rprec * ( 1.0_rprec - cos (pi * real (i - istart, rprec)  &
        / (iend - istart)) )
endif

return
end function fringe_weighting

end module fringe_util
