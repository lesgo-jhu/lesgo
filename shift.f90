!!
!!  Copyright (C) 2010-2020 Johns Hopkins University
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

!*******************************************************************************
module shift
!*******************************************************************************
! Spanwise shifting of variables in the domain

use types, only : rprec

implicit none

save
private

! Use shifting to eliminate streaks in time average
! The domain will be shifted by shift_n spanwise grid point every shift_base steps
logical, public :: use_shift = .false.
integer, public :: shift_base = 1000
integer, public :: shift_n = 1

! Dummy variables
real(rprec), dimension(:,:,:), allocatable :: dummy

public :: shift_domain, shift_init

contains

!*******************************************************************************
subroutine shift_init()
!*******************************************************************************
use param, only : ld, ny, nz, lbz

! Only actually allocate if the shifting is used
if (use_shift) allocate( dummy(ld, ny, lbz:nz) )

end subroutine shift_init

!*******************************************************************************
subroutine shift_domain()
!*******************************************************************************
use param, only : jt_total, coord
use sim_param, only : u, v, w, RHSx, RHSy, RHSz
use sgs_param, only : F_LM, F_MM, F_QN, F_NN

! Shift the domain in the y (spanwise) direction
if (use_shift) then
    if ( (shift_base < 0 .and. jt_total == 1) .or.                             &
        (modulo (jt_total, shift_base) == 0 .and. shift_base >0)) then

        if (coord == 0) then
            write(*,*) 'Shifting domain ', shift_n, ' grid point in spanwise direction'
        end if
        call shift_variable(u)
        call shift_variable(v)
        call shift_variable(w)

        call shift_variable(RHSx)
        call shift_variable(RHSy)
        call shift_variable(RHSz)

        call shift_variable(F_LM)
        call shift_variable(F_MM)
        call shift_variable(F_QN)
        call shift_variable(F_NN)
    end if
end if

end subroutine shift_domain

!*******************************************************************************
subroutine shift_variable(var)
!*******************************************************************************
use param, only : ny
real(rprec), dimension(:,:,:), intent(inout) :: var

if (shift_n >= ny/2) then
    dummy(:,:,:) = var(:,:,:)
    var(:,shift_n+1:ny,:) = dummy(:,1:ny-shift_n,:)
    var(:,1:shift_n,:) = dummy(:,ny-shift_n+1:ny,:)
else
    dummy(:,:,:) = var(:,:,:)
    var(:,shift_n+1:ny,:) = dummy(:,1:ny-shift_n,:)
    var(:,1:shift_n,:) = dummy(:,ny-shift_n+1:ny,:)
endif

end subroutine shift_variable

end module shift
