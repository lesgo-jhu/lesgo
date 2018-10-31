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
use types, only : rprec
use messages
implicit none

private
public grid_t

type grid_t
    integer :: Nx, Ny, Nz, Nz_tot, jzmin, jzmax, ld, Nkx
    real(rprec) :: dx, dy, dz, L_x, L_y, L_z
    logical :: uniform_spacing
    real(rprec), dimension(:), allocatable :: x, y, z, zw
    integer, dimension(:), allocatable :: autowrap_i, autowrap_j
contains

end type grid_t

interface grid_t
    module procedure constructor
end interface grid_t

contains

!*******************************************************************************
function constructor(Nx, Ny, Nz_tot, L_x, L_y, L_z, nproc, coord, uniform_spacing)&
    result(this)
!*******************************************************************************
!  Constructor for grid_t

type(grid_t) :: this
integer, intent(in) :: Nx, Ny, Nz_tot, nproc, coord
real(rprec), intent(in) :: L_x, L_y, L_z
logical, intent(in) :: uniform_spacing
integer :: i,j,k
integer :: ival_read
real(rprec) :: val_read
! Thresh hold for evaluating differences in floating point values.
real(rprec), parameter :: thresh = 1.0e-6_rprec

! set input arguments
this%Nx = Nx
this%Ny = Ny
this%Nz_tot = Nz_tot
this%L_x = L_x
this%L_y = L_y
this%L_z = L_z
this%uniform_spacing = uniform_spacing

! Set the processor owned vertical grid spacing
this%Nz = floor(real(this%Nz_tot, rprec)/nproc) + 1

! Recompute nz_tot to be compliant with computed nz
ival_read = this%Nz_tot
this%Nz_tot = ( this%Nz - 1 ) * nproc + 1
if (coord == 0 .AND. ival_read /= this%Nz_tot )                             &
   write(*,*) 'Reseting Nz (total) to: ', this%Nz_tot

! Grid size for FFTs
this%Nkx = this%Nx / 2 + 1
this%ld = 2 * this%Nkx

! Grid spacing (x direction)
this%dx = this%L_x / this%Nx

! Check if we are to enforce uniform grid spacing
if (this%uniform_spacing) then
    ! Adjust L_y
    val_read = this%L_y
    this%L_y = this%Ny * this%dx
    if (coord == 0 .AND. abs( val_read - this%L_y ) >= thresh)              &
        call mesg( 'grid_t:constructor', 'Reseting Ly to: ', this%L_y )

    ! Adjust L_z
    val_read = this%L_z
    this%L_z = (nz_tot - 1 ) * this%dx
    if (coord == 0 .AND. abs( val_read - this%L_z ) >= thresh)              &
        call mesg( 'grid_t:constructor', 'Reseting Lz to: ', this%L_z )
endif

! Grid spacing (y and z directions)
this%dy = this%L_y / this%Ny
this%dz = this%L_z / ( this%Nz_tot - 1 )

!  x and y go to nx+1, ny+1 respectively for adding
!  the buffered points for periodicity
allocate(this%x(this%Nx+1),this%y(this%Ny+1))
allocate(this%z(0:this%Nz), this%zw(0:this%Nz))
allocate(this%autowrap_i(0:this%Nx+1), this%autowrap_j(0:this%Ny+1))

do k = 0, this%Nz
#ifdef PPMPI
    this%z(k) = (coord*(this%Nz-1) + k - 0.5_rprec) * this%dz
#else
    this%z(k) = (k - 0.5_rprec) * this%dz
#endif
enddo

do j = 1, this%Ny+1
    this%y(j) = (j-1)*this%dy
enddo

do i = 1, this%Nx+1
    this%x(i) = (i-1)*this%dx
enddo

this%zw = this%z - this%dz/2._rprec

! Set index autowrapping arrays
this%autowrap_i(0) = this%Nx
this%autowrap_j(0) = this%Ny
this%autowrap_i(this%Nx+1) = 1
this%autowrap_j(this%Ny+1) = 1
do i=1,this%Nx; this%autowrap_i(i) = i; enddo
do j=1,this%Ny; this%autowrap_j(j) = j; enddo

! Set jzmin and jzmax - the levels that this processor "owns"
#ifdef PPMPI
if (coord == 0) then
    this%jzmin = 0
    this%jzmax = this%Nz-1
elseif (coord == nproc-1) then
    this%jzmin = 1
    this%jzmax = this%Nz
else
    this%jzmin = 1
    this%jzmax = this%Nz-1
endif
#else
this%jzmin = 1
this%jzmax = this%Nz
#endif

end function constructor

end module grid_m
