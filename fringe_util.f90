!!
!!  Copyright (C) 2011-2013  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
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
