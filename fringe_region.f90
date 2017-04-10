!!
!!  Copyright (C) 2011-2017  Johns Hopkins University
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
module fringe_region
!*******************************************************************************
!
!
!
!
! Determine the weights 
!
!   u_fringe = alpha * u + beta * u_sample
!   alpha = 1-beta
!
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
implicit none

private
public fringe_region_t

type fringe_region_t
    ! check if initialized
    
    ! number of points in the fringe region
    integer :: n
    ! fringe weights
    real(rprec), allocatable, dimension(:) :: alpha, beta
    ! wrapped index
    integer, allocatable, dimension(:) :: iwrap
    ! velocity sample block
    real(rprec), allocatable, dimension(:,:,:) :: u, v, w
contains
    procedure, private :: initialize
    procedure, public :: sample_vel
    procedure, public :: apply_vel
    final :: destructor
end type fringe_region_t

! constructor for fringe_region_t
interface fringe_region_t
    module procedure constructor_len
    module procedure constructor_num
end interface fringe_region_t

! Fringe region
type(fringe_region_t), public :: fringe
! recycling region, if using shifted periodic boundary conditions
type(fringe_region_t), public :: recycl

contains

!*******************************************************************************
function constructor_len(fr_end, fr_len) result(this)
!*******************************************************************************
! Constructs the fringe region type.
use param, only : nx
implicit none
real(rprec), intent(in) :: fr_end, fr_len
type(fringe_region_t) :: this
integer :: istart, iend

! Set indices for start and end of region
iend = floor (fr_end * nx + 1.0_rprec)
istart = floor ((fr_end - fr_len) * nx + 1.0_rprec)

! Number of points in the fringe region
this%n = iend - istart

call this%initialize(istart)

end function constructor_len

!*******************************************************************************
function constructor_num(fr_end, fr_num) result(this)
!*******************************************************************************
! Constructs the fringe region type.
use param, only : nx, L_x
implicit none
real(rprec), intent(in) :: fr_end
integer, intent(in) :: fr_num
type(fringe_region_t) :: this
integer :: istart, iend

! Number of points in the fringe region
this%n = fr_num

! Set indices for start and end of region
iend = floor (fr_end * nx + 1.0_rprec)
istart = iend - this%n

call this%initialize(istart)

end function constructor_num

!*******************************************************************************
subroutine initialize(this, istart)
!*******************************************************************************
use param, only : nx, ny, nz, pi
implicit none
class(fringe_region_t), intent(inout) :: this
integer, intent(in) :: istart
integer :: i, np

! number of points in the non-plateau region
np = 3*this%n/4

! Calculate weighting functions. beta is a cosine profile with plateau
if (allocated(this%alpha)) deallocate(this%alpha)
if (allocated(this%beta)) deallocate(this%beta)
allocate(this%alpha(this%n))
allocate(this%beta(this%n))
this%beta = 1._rprec
do i = 1, np
    this%beta(i) = 0.5_rprec * (1.0_rprec - cos(pi*real(i, rprec)/np))
end do
this%alpha = 1.0_rprec - this%beta

! Wrapped index
if (allocated(this%iwrap)) deallocate(this%iwrap)
allocate(this%iwrap(this%n))
do i = 1, this%n
   this%iwrap(i) = modulo(istart+i-1, nx ) + 1
end do

! Allocate the sample block
if (allocated(this%u)) deallocate(this%u)
if (allocated(this%v)) deallocate(this%v)
if (allocated(this%w)) deallocate(this%w)
allocate(this%u(this%n, ny, nz))
allocate(this%v(this%n, ny, nz))
allocate(this%w(this%n, ny, nz))

end subroutine initialize

!*******************************************************************************
subroutine destructor(this)
!*******************************************************************************
! Deallocate any allocated arrays
implicit none
type(fringe_region_t), intent(inout) :: this

if (allocated(this%alpha)) deallocate(this%alpha)
if (allocated(this%beta)) deallocate(this%beta)
if (allocated(this%iwrap)) deallocate(this%iwrap)
if (allocated(this%u)) deallocate(this%u)
if (allocated(this%v)) deallocate(this%v)
if (allocated(this%w)) deallocate(this%w)

end subroutine destructor

!*******************************************************************************
subroutine sample_vel(this)
!*******************************************************************************
! Takes the velocities from the flow field and copies to the buffer arrays 
! in fringe_region_t
use sim_param, only : ny, nz, u, v, w
implicit none
class(fringe_region_t) :: this

this%u(:,:,:) = u(this%iwrap, 1:ny, 1:nz)
this%v(:,:,:) = v(this%iwrap, 1:ny, 1:nz)
this%w(:,:,:) = w(this%iwrap, 1:ny, 1:nz)

end subroutine sample_vel

!*******************************************************************************
subroutine apply_vel(this)
!*******************************************************************************
! Takes the velocities in the buffer arrays of fringe_region_t and applies the 
! velocities using the weighting functions
use sim_param, only : ny, nz, u, v, w
implicit none
class(fringe_region_t) :: this
integer :: j, k

do k=1,nz
    do j=1,ny
        u(this%iwrap,j,k) = this%alpha * u(this%iwrap, j, k)                   &
            + this%beta * this%u(:, j, k)
        v(this%iwrap,j,k) = this%alpha * v(this%iwrap, j, k)                   &
            + this%beta * this%u(:, j, k)
        w(this%iwrap,j,k) = this%alpha * w(this%iwrap, j, k)                   &
            + this%beta * this%u(:, j, k)
    enddo
enddo

end subroutine apply_vel

end module fringe_region
