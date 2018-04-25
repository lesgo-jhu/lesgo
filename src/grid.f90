!!
!!  Copyright (C) 2009-2017  Johns Hopkins University
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
module grid_m
!*******************************************************************************
use param, only : rprec
implicit none
save
private

public grid

type grid_t
    logical :: built
    real(rprec), pointer, dimension(:) :: x, y, z, zw
    integer, pointer, dimension(:) :: autowrap_i, autowrap_j
contains
    procedure, public :: build
end type grid_t

! The uv grid
type(grid_t) :: grid

contains

!*******************************************************************************
subroutine build(this)
!*******************************************************************************
!  This subroutine creates the uv grid for the domain.

use param, only : nx, ny, nz, jzmin, jzmax, dx, dy, dz, lbz
#ifdef PPMPI
use param, only : nproc, coord
#endif

class(grid_t) :: this
integer :: i,j,k
real(rprec), pointer, dimension(:) :: x, y, z, zw
integer, pointer, dimension(:) :: autowrap_i, autowrap_j

nullify(x, y, z, zw)
nullify(autowrap_i,autowrap_j)

!  x and y go to nx+1, ny+1 respectively for adding
!  the buffered points for periodicity
allocate(grid % x(nx+1),grid % y(ny+1))
allocate(grid % z(lbz:nz), grid % zw(lbz:nz))
allocate(grid % autowrap_i(0:nx+1), grid % autowrap_j(0:ny+1))

! Initialize built
grid % built = .false.

! Set pointers
x => this % x
y => this % y
z => this % z
zw => this % zw

autowrap_i => this % autowrap_i
autowrap_j => this % autowrap_j

do k = lbz, nz
#ifdef PPMPI
    z(k) = (coord*(nz-1) + k - 0.5_rprec) * dz
#else
    z(k) = (k - 0.5_rprec) * dz
#endif
enddo

do j = 1, ny+1
    y(j) = (j-1)*dy
enddo

do i = 1, nx+1
    x(i) = (i - 1)*dx
enddo

zw = z - dz/2._rprec

! Set index autowrapping arrays
autowrap_i(0) = nx
autowrap_j(0) = ny
autowrap_i(nx+1) = 1
autowrap_j(ny+1) = 1
do i=1,nx; autowrap_i(i) = i; enddo
do j=1,ny; autowrap_j(j) = j; enddo

this % built = .true.

nullify(x,y,z,zw)
nullify(autowrap_i,autowrap_j)

! Set jzmin and jzmax - the levels that this processor "owns"
#ifdef PPMPI
if (coord == 0) then
    jzmin = 0
    jzmax = nz-1
elseif (coord == nproc-1) then
    jzmin = 1
    jzmax = nz
else
    jzmin = 1
    jzmax = nz-1
endif
#else
jzmin = 1
jzmax = nz
#endif

end subroutine build

end module grid_m
