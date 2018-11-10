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
use types, only : rprec, cprec
use messages
#ifdef PPMPI
use mpi
#endif
implicit none

private
public grid_t

type grid_t
    integer :: Nx, Ny, Nz, Nz_tot, jzmin, jzmax, ld, Nkx
    real(rprec) :: dx, dy, dz, L_x, L_y, L_z
    logical :: uniform_spacing
    real(rprec), dimension(:), allocatable :: x, y, z, zw
    real(rprec), allocatable, dimension(:) :: kx, ky
    integer, dimension(:), allocatable :: autowrap_i, autowrap_j
    integer :: nproc=1, rank=0, coord=0
    integer :: comm, up, down
    integer, allocatable, dimension(:) :: rank_of_coord, coord_of_rank
    integer, dimension(:,:), allocatable :: nyquist
    complex(cprec), dimension(:,:), allocatable :: G_test
#ifdef PPCGNS
    integer :: cgnsParallelComm
#endif
contains

end type grid_t

interface grid_t
    module procedure constructor
end interface grid_t

contains

!*******************************************************************************
function constructor(Nx, Ny, Nz_tot, L_x, L_y, L_z, uniform_spacing, globalComm)&
    result(this)
!*******************************************************************************
!  Constructor for grid_t

type(grid_t) :: this
integer, intent(in) :: Nx, Ny, Nz_tot, globalComm
real(rprec), intent(in) :: L_x, L_y, L_z
logical, intent(in) :: uniform_spacing
integer :: i, j, k, ierr
! Thresh hold for evaluating differences in floating point values.
integer :: coords(1)
real(rprec) :: pi, delta_test
! the implicit filter (1=grid size)
integer, parameter :: filter_size=1
! alpha is ratio of test filter to grid filter widths
real(rprec) :: alpha_test = 2.0_rprec * filter_size
real(rprec), dimension(:,:), allocatable :: k2

! set input arguments
this%Nx = Nx
this%Ny = Ny
this%Nz_tot = Nz_tot
this%L_x = L_x
this%L_y = L_y
this%L_z = L_z
this%uniform_spacing = uniform_spacing

! pi
pi = acos(-1._rprec)

#ifdef PPMPI
! Set the local communicator
! set up a 1d cartesian topology
call mpi_comm_size (globalComm, this%nproc, ierr)
call mpi_cart_create(globalComm, 1, (/ this%nproc /), (/ .false. /),                &
    .false., this%comm, ierr)

! slight problem here for ghost layers:
! u-node info needs to be shifted up to proc w/ rank "up",
! w-node info needs to be shifted down to proc w/ rank "down"
call mpi_cart_shift(this%comm, 0, 1, this%down, this%up, ierr)
call mpi_comm_size (this%comm, this%nproc, ierr)
call mpi_comm_rank(this%comm, this%rank, ierr)
call mpi_cart_coords(this%comm, this%rank, 1, coords, ierr)

! use coord (NOT rank) to determine global position
this%coord = coords(1)

! rank->coord and coord->rank conversions
allocate(this%rank_of_coord(0:this%nproc-1), this%coord_of_rank(0:this%nproc-1))
do i = 0, this%nproc-1
    call mpi_cart_rank(this%comm, (/ i /), this%rank_of_coord(i), ierr)
    call mpi_cart_coords(this%comm, i, 1, coords, ierr)
    this%coord_of_rank(i) = coords(1)
end do

#ifdef PPCGNS
! Set the CGNS parallel Communicator
this%cgnsParallelComm = globalComm

! Set the parallel communicator
call cgp_mpi_comm_f(this%cgnsParallelComm, ierr)
#endif

#endif

! Set the processor owned vertical grid spacing
this%Nz = floor(real(this%Nz_tot, rprec)/this%nproc) + 1

! Recompute nz_tot to be compliant with computed nz
this%Nz_tot = ( this%Nz - 1 ) * this%nproc + 1
if (this%coord == 0) write(*,*) 'Nz (total) = ', this%Nz_tot

! Grid spacing (x direction)
this%dx = this%L_x / this%Nx

! Check if we are to enforce uniform grid spacing
if (this%uniform_spacing) then
    ! Adjust L_y
    this%L_y = this%Ny * this%dx
    ! Adjust L_z
    this%L_z = (nz_tot - 1 ) * this%dx
end if

if (this%coord == 0) then
    write(*,*) "Lx = ", this%L_x
    write(*,*) "Ly = ", this%L_y
    write(*,*) "Lz = ", this%L_z
end if

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
    this%z(k) = (this%coord*(this%Nz-1) + k - 0.5_rprec) * this%dz
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

! Wavenumber grid
this%Nkx = this%Nx / 2 + 1
this%ld = 2 * this%Nkx
allocate(this%kx(this%Nkx))
allocate(this%ky(this%Ny))

do i = 1, this%Nkx
    this%kx(i) = real(i-1, kind=rprec)
end do
do i = 1, this%Ny
    this%ky(i) = real(modulo(i - 1 + this%Ny/2,this%Ny) - this%Ny/2, kind=rprec)
end do

! Aspect ratio change
this%kx = 2._rprec*pi/this%L_x*this%kx
this%ky = 2._rprec*pi/this%L_y*this%ky

! Find nyqusit frequencies
allocate(this%nyquist(this%Ny+this%Nkx,3))
do i = 1, this%Nkx
    this%nyquist(i+this%Ny,1) = i
    this%nyquist(i+this%Ny,2) = mod(this%Ny/2+1, this%Ny)
end do

do i = 1, this%Ny
    this%nyquist(i,1) = this%Nkx
    this%nyquist(i,2) = i
end do

! Set index autowrapping arrays
this%autowrap_i(0) = this%Nx
this%autowrap_j(0) = this%Ny
this%autowrap_i(this%Nx+1) = 1
this%autowrap_j(this%Ny+1) = 1
do i=1,this%Nx; this%autowrap_i(i) = i; enddo
do j=1,this%Ny; this%autowrap_j(j) = j; enddo

! Set jzmin and jzmax - the levels that this processor "owns"
#ifdef PPMPI
if (this%coord == 0) then
    this%jzmin = 0
    this%jzmax = this%Nz-1
elseif (this%coord == this%nproc-1) then
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

! Create test filter

! Allocate the arrays
allocate( this%G_test(this%Nkx,this%Ny), k2(this%Nkx,this%Ny) )
do i = 1, this%Nkx
    k2(i,:) = this%kx(i)**2 + this%ky**2
end do

! Include the normalization for the forward FFT
this%G_test = 1._rprec/(this%Nx*this%Ny)

! Filter characteristic width
! "2d-delta", not full 3d one
delta_test = alpha_test * sqrt(this%dx*this%dy)

! Calculate the kernel
! spectral cutoff filter
! if(ifilter==1) then
where (k2 >= (pi/(delta_test))**2) this%G_test = 0._rprec

! ! Gaussian filter
! else if(ifilter==2) then
!     G_test=exp(-(delta_test)**2*k2/(4._rprec*6._rprec))*G_test
!
! ! Top-hat (Box) filter
! else if(ifilter==3) then
!     G_test = (sin(kx*delta_test/2._rprec)*sin(ky*delta_test/2._rprec)+1E-8)/   &
!         (kx*delta_test/2._rprec*ky*delta_test/2._rprec+1E-8)*G_test
! endif

! ! since our k2 has zero at Nyquist, we have to do this by hand
! G_test(lh,:) = 0._rprec
! G_test(:,ny/2+1) = 0._rprec

end function constructor

end module grid_m
