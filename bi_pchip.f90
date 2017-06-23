!!
!!  Copyright (C) 2017 Johns Hopkins University
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
module bi_pchip
!*******************************************************************************
! The bi_pchip_t class performs cubic spline interpolation for a 2D
! function v(x,y). The evaluated data is passed on construction. Interpolation
! is evaluated for sample points (xq,yq) passed as values or arrays. The
! first partial derivatives can also be returned as optional arguments.
!
! The resulting interpolation has continuous first and second derivatives.
! Linear extrapolation is used for evaluation points (xq,yq) outside of the
! interval of (x,y)'s passed at construction. The second partial derivatives
! vanish at the boundaries
use types, only : rprec
use messages
use pchip
implicit none

private
public bi_pchip_t

type :: bi_pchip_t
    type(pchip_t), dimension(:), allocatable :: xspline, yspline
    real(rprec), dimension(:), allocatable :: x, y
    integer :: N, M
    character(9) :: type_name = "bi_pchip_t"
contains
    procedure, private :: interp_scalar
    procedure, private :: interp_array
    generic, public :: interp => interp_scalar, interp_array
end type bi_pchip_t

interface bi_pchip_t
    module procedure constructor
end interface bi_pchip_t

contains

!*******************************************************************************
function constructor(x, y, v) result(this)
!*******************************************************************************
implicit none
type(bi_pchip_t) :: this
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
character(14) :: proc_name = "constructor"
integer :: i

! Check size of inputs
this%N = size(x)
this%M = size(y)
if ( size(v, 1) /= this%N .or. size(v, 2) /= this%M) then
    call error(this%type_name // '.' // proc_name,'v must be an NxM array')
end if

! Allocate class members
allocate( this%xspline(this%M) )
allocate( this%yspline(this%N) )
allocate( this%x(this%N) )
allocate( this%y(this%M) )
this%x = x
this%y = y

! Create a pchip for each x location along the y coordinate
do i = 1, this%N
    this%yspline(i) = pchip_t(y, v(i,:))
end do

! Create a pchip for each y location along the x coordinate
do i = 1, this%M
    this%xspline(i) = pchip_t(x, v(:,i))
end do

end function constructor

!*******************************************************************************
subroutine interp_scalar(this, xq, yq, vq, vqpx, vqpy)
!*******************************************************************************
implicit none
class(bi_pchip_t) :: this
real(rprec), intent(in) :: xq, yq
real(rprec), intent(out) :: vq
real(rprec), intent(out), optional :: vqpx, vqpy
real(rprec), dimension(:), allocatable :: vqx, vqy
type(pchip_t) :: xqspline, yqspline
integer :: i

! First evaluate the value along the x-direction
! Also evaluate partial derivative in x if needed
allocate( vqy(this%N) )
do i = 1, this%N
    call this%yspline(i)%interp(yq, vqy(i))
end do
xqspline = pchip_t(this%x, vqy)
if ( present(vqpx) ) then
    call xqspline%interp(xq, vq, vqpx)
else
    call xqspline%interp(xq, vq)
end if
deallocate(vqy)

! Evaluate partial derivative in y if needed
! If this interpolation is done, the value of vq will also be recalculated, but
! will not change
if ( present(vqpy) ) then
    allocate( vqx(this%M) )
    do i = 1, this%M
        call this%xspline(i)%interp(xq, vqx(i))
    end do
    yqspline = pchip_t(this%y, vqx)
    call yqspline%interp(yq, vq, vqpy)
    deallocate(vqx)
end if

end subroutine interp_scalar

!*******************************************************************************
subroutine interp_array(this, xq, yq, vq, vqpx, vqpy)
!*******************************************************************************
implicit none
class(bi_pchip_t) :: this
real(rprec), dimension(:), intent(in) :: xq, yq
real(rprec), dimension(:), intent(out) :: vq
real(rprec), dimension(:), intent(out), optional :: vqpx, vqpy
integer :: i, Nq

! Check input dimensions
Nq = size(xq)
if (size(yq) /= Nq) then
    call error('bi_pchip_t/interp_array','xq and yq must be the same size')
end if

! Determine which call is needed, and then iterate through the arrays
if ( present(vqpx) .and. present(vqpy) ) then
    do i = 1, Nq
        call this%interp(xq(i), yq(i), vq(i), vqpx(i), vqpy(i))
    end do
else if ( present(vqpx) ) then
    do i = 1, Nq
        call this%interp(xq(i), yq(i), vq(i), vqpx(i))
    end do
else if ( present(vqpy) ) then
    do i = 1, Nq
        call this%interp(xq(i), yq(i), vq(i), vqpy=vqpy(i))
    end do
else
    do i = 1, Nq
        call this%interp(xq(i), yq(i), vq(i))
    end do
end if

end subroutine interp_array

end module bi_pchip
