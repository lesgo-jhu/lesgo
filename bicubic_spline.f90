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
module bicubic_spline
!*******************************************************************************
! The bicubic_spline_t class performs cubic spline interpolation for a 2D
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
use cubic_spline
implicit none

private
public bicubic_spline_t

type :: bicubic_spline_t
    type(cubic_spline_t), dimension(:), allocatable :: xspline, yspline
    real(rprec), dimension(:), allocatable :: x, y
    character(:), allocatable :: x_low_bc, x_high_bc, y_low_bc, y_high_bc
    type(cubic_spline_t) :: x_low_fspline, x_high_fspline
    type(cubic_spline_t) :: y_low_fspline, y_high_fspline
    integer :: N, M
    character(16) :: type_name = "bicubic_spline_t"
contains
    procedure, private :: interp_scalar
    procedure, private :: interp_array
    generic, public :: interp => interp_scalar, interp_array
end type bicubic_spline_t

interface bicubic_spline_t
    module procedure constructor
end interface bicubic_spline_t

contains

!*******************************************************************************
function constructor(x, y, v, x_low_bc, x_high_bc, y_low_bc, y_high_bc,        &
        x_low_f, x_high_f, y_low_f, y_high_f) result(this)
!*******************************************************************************
implicit none
type(bicubic_spline_t) :: this
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
character(*), intent(in), optional :: x_low_bc, x_high_bc, y_low_bc, y_high_bc
real(rprec), dimension(:), intent(in), optional  :: x_low_f, x_high_f
real(rprec), dimension(:), intent(in), optional  :: y_low_f, y_high_f
real(rprec), dimension(:), allocatable :: vbc
character(14) :: proc_name = "constructor"
integer :: i
real(rprec) :: low_f, high_f

! Check size of inputs
this%N = size(x)
this%M = size(y)
if ( size(v, 1) /= this%N .or. size(v, 2) /= this%M) then
    call error(this%type_name // '.' // proc_name,'v must be an NxM array')
end if

! Check size of optional boundary condition inputs
if ( present(x_low_f) ) then
    if ( size(x_low_f) /= this%M ) then
        call error(this%type_name // '.' // proc_name,                         &
            'x_low_f must be an array of length M')
    end if
end if
if ( present(x_high_f) ) then
    if ( size(x_high_f) /= this%M ) then
        call error(this%type_name // '.' // proc_name,                         &
            'x_high_f must be an array of length M')
    end if
end if
if ( present(y_low_f) ) then
    if ( size(y_low_f) /= this%N ) then
        call error(this%type_name // '.' // proc_name,                         &
            'y_low_f must be an array of length N')
    end if
end if
if ( present(y_high_f) ) then
    if ( size(y_high_f) /= this%N ) then
        call error(this%type_name // '.' // proc_name,                         &
            'y_high_f must be an array of length N')
    end if
end if

! Allocate class members
allocate( this%xspline(this%M) )
allocate( this%yspline(this%N) )
allocate( this%x(this%N) )
allocate( this%y(this%M) )
this%x = x
this%y = y

! Specify boundary condition types
this%x_low_bc = 'natural'
this%x_high_bc = 'natural'
this%y_low_bc = 'natural'
this%y_high_bc = 'natural'
if ( present(x_low_bc) ) this%x_low_bc = x_low_bc
if ( present(x_high_bc) ) this%x_high_bc = x_high_bc
if ( present(y_low_bc) ) this%y_low_bc = y_low_bc
if ( present(y_high_bc) ) this%y_high_bc = y_high_bc

! Create a cubic spline for each x location along the y coordinate
low_f = 0._rprec
high_f = 0._rprec
do i = 1, this%N
    if ( present(y_low_f) ) low_f = y_low_f(i)
    if ( present(y_high_f) ) high_f = y_high_f(i)
    this%yspline(i) = cubic_spline_t(y, v(i,:), this%y_low_bc, this%y_high_bc, &
        low_f, high_f)
end do

! Create a cubic spline for each y location along the x coordinate
low_f = 0._rprec
high_f = 0._rprec
do i = 1, this%M
    if ( present(x_low_f) ) low_f = x_low_f(i)
    if ( present(x_high_f) ) high_f = x_high_f(i)
    this%xspline(i) = cubic_spline_t(x, v(:,i), this%x_low_bc, this%x_high_bc, &
        low_f, high_f)
end do

! Create cubic splines for x boundary conditions (evaluated along y axis)
allocate( vbc(this%M) )
if ( present(x_low_f) ) then
    vbc = x_low_f
else
    vbc = 0._rprec
end if
this%x_low_fspline = cubic_spline_t(y, vbc, 'clamped', 'clamped')
if ( present(x_high_f) ) then
    vbc = x_high_f
else
    vbc = 0._rprec
end if
this%x_high_fspline = cubic_spline_t(y, vbc, 'clamped', 'clamped')
deallocate(vbc)

! Create cubic splines for y boundary conditions (evaluated along x axis)
allocate( vbc(this%N) )
if ( present(y_low_f) ) then
    vbc = y_low_f
else
    vbc = 0._rprec
end if
this%y_low_fspline = cubic_spline_t(x, vbc, 'clamped', 'clamped')
if ( present(y_high_f) ) then
    vbc = y_high_f
else
    vbc = 0._rprec
end if
this%y_high_fspline = cubic_spline_t(x, vbc, 'clamped', 'clamped')
deallocate(vbc)

end function constructor

!*******************************************************************************
subroutine interp_scalar(this, xq, yq, vq, vqpx, vqpy)
!*******************************************************************************
implicit none
class(bicubic_spline_t) :: this
real(rprec), intent(in) :: xq, yq
real(rprec), intent(out) :: vq
real(rprec), intent(out), optional :: vqpx, vqpy
real(rprec), dimension(:), allocatable :: vqx, vqy
type(cubic_spline_t) :: xqspline, yqspline
real(rprec) :: low_f, high_f
integer :: i

! First evaluate the value along the x-direction
! Also evaluate partial derivative in x if needed
allocate( vqy(this%N) )
do i = 1, this%N
    call this%yspline(i)%interp(yq, vqy(i))
end do
call this%x_low_fspline%interp(yq, low_f)
call this%x_high_fspline%interp(yq, high_f)
xqspline = cubic_spline_t(this%x, vqy, this%x_low_bc, this%x_high_bc,          &
    low_f, high_f)
if ( present(vqpx) ) then
    call xqspline%interp(xq+0.1_rprec, vq, vqpx)
    write(*,*) low_f, high_f, xq+0.1_rprec, vq, vqpx
    write(*,*) xqspline%v
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
    call this%y_low_fspline%interp(xq, low_f)
    call this%y_high_fspline%interp(xq, high_f)
    yqspline = cubic_spline_t(this%y, vqx, this%y_low_bc, this%y_high_bc,      &
        low_f, high_f)
    write(*,*) low_f, high_f
    call yqspline%interp(yq, vq, vqpy)
    deallocate(vqx)
end if

end subroutine interp_scalar

!*******************************************************************************
subroutine interp_array(this, xq, yq, vq, vqpx, vqpy)
!*******************************************************************************
implicit none
class(bicubic_spline_t) :: this
real(rprec), dimension(:), intent(in) :: xq, yq
real(rprec), dimension(:), intent(out) :: vq
real(rprec), dimension(:), intent(out), optional :: vqpx, vqpy
integer :: i, Nq

! Check input dimensions
Nq = size(xq)
if (size(yq) /= Nq) then
    call error('bicubic_spline_t/interp_array',                                &
        'xq and yq must be the same size')
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

end module bicubic_spline
