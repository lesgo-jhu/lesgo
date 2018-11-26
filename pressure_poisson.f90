!!
!!  Copyright (C) 2009-2018  Johns Hopkins University
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
module pressure_poisson_m
!*******************************************************************************
use tridiagonal_m
use mpi
use types, only : rprec, cprec
use param, only : grid, BOGUS

implicit none

private
public :: pressure_poisson, pressure_poisson_init

type(tridiagonal_t), dimension(:,:), allocatable :: td
integer :: Nind, istart, iend

contains

!*******************************************************************************
subroutine pressure_poisson
!*******************************************************************************



end subroutine pressure_poisson


!*******************************************************************************
subroutine pressure_poisson_init()
!*******************************************************************************
integer :: i, j, Nplus
real(rprec), dimension(:), allocatable :: a, b, c

#ifdef PPMPI
! Split wavenumbers across the x-direction
! Baseline number of wavenumbers per processors
Nind = grid%Nkx / grid%nproc
! Number of processors that will need to have Nind+1 wavenumbers to sum to total
Nplus = grid%Nkx - Nind*grid%nproc
! Update Nind if we're on the right coord
if (grid%coord < Nplus) Nind = Nind + 1
! Find start and stopping indices
if (grid%coord < Nplus) then
    istart = grid%coord*Nind + 1
    iend = (grid%coord+1)*Nind
else
    istart = (Nind+1)*Nplus + (grid%coord-Nplus)*Nind + 1
    iend = (Nind+1)*Nplus + (grid%coord-Nplus+1)*Nind
end if

#else
Nind = grid%Nkx
istart = 1
iend = this%grid%Nkx
#endif

! Allocate
allocate( td(Nind, grid%ny) )
allocate( a(0:grid%nz) )
allocate( b(0:grid%nz) )
allocate( c(0:grid%nz) )

! Create TDMAs
a(:) = 1._rprec/grid%dz**2
c(:) = 1._rprec/grid%dz**2
#ifdef PPSAFETYMODE
a(0) = BOGUS
c(grid%nz) = BOGUS
#endif
do i = 1, Nind
    do j = 1, grid%ny
        b(:) = -2._rprec/grid%dz**2 - grid%kx(istart+i-1)**2 - grid%ky(j)**2
        td(i,j) = tridiagonal_t(a, b, c)
    end do
end do

! Deallocate
deallocate(a, b, c)


end subroutine pressure_poisson_init




end module pressure_poisson_m
