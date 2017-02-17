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
module pchip
!*******************************************************************************
use types, only : rprec
use messages
implicit none

private
public pchip_t

type :: pchip_t
    real(rprec), dimension(:), allocatable :: x, v, vp
    integer :: N
    character(7) :: type_name = "pchip_t"
contains
    procedure, private :: interp_scalar
    procedure, private :: interp_array
    generic, public :: interp => interp_scalar, interp_array
end type pchip_t

interface pchip_t
    module procedure constructor
end interface pchip_t

contains

!*******************************************************************************
function constructor(x, v) result(this)
!*******************************************************************************
implicit none

type(pchip_t) :: this
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), dimension(:), allocatable :: delta
logical, dimension(:), allocatable :: zeroed
real(rprec) :: a, b, tau
character(14) :: proc_name = "constructor"
integer :: i

! Set the size of the matrix
this%N = size(x)

! Check that all input arguments are the same size
if ( size(v) /= this%N ) then
    call error(this%type_name // '.' // proc_name,                             &
        'x and v must be the same size')
end if

! Check that x is sorted
do i = 2, this%N
    if ( x(i) < x(i-1) ) then
        call error(this%type_name // '.' // proc_name, 'x must be increasing')
    end if
end do

! Allocate and assign variables
allocate( this%x(this%N) )
allocate( this%v(this%N) )
allocate( this%vp(this%N) )
this%x = x
this%v = v
this%vp = 0._rprec

! Calculate secant lines between points
allocate( Delta(this%N - 1) )
do i  = 1, this%N-1
    Delta(i) = (this%v(i+1) - this%v(i)) / (this%x(i+1) - this%x(i))
end do

! Set initial guess for the derivatives
!   end points use one-sided differences
!   if the secant lines are different signs, set to zero
!   interior points just average secant lines on both sides
this%vp(1) = Delta(1)
this%vp(this%N) = Delta(this%N-1)
do i = 2, this%N-1
    if ( (Delta(i-1) * Delta(i)) <= 0._rprec ) then
        this%vp(i) = 0._rprec
    else
        this%vp(i) = 0.5_rprec * ( Delta(i-1) + Delta(i) )
    end if
end do

! An array to keep track of zeroed derivative segments
allocate( zeroed(this%N-1) )
zeroed = .false.

! Check if a segment has zero slope and set the derivatives at the end to zero
do i = 1, this%N-1
    if (Delta(i) == 0._rprec) then
        zeroed(i) = .true.
        this%vp(i) = 0._rprec
        this%vp(i+1) = 0._rprec
    end if
end do

! Check that the derivates at the end of each segment have the same sign as the
! secant line. Otherwise, set the derivatives to zero
do i = 1, this%N-1
    if ( .not.zeroed(i) ) then
        a = this%vp(i) / Delta(i) 
        b = this%vp(i+1) / Delta(i) 
        if ( a < 0._rprec .or. b < 0._rprec ) this%vp(i) = 0._rprec
    end if
end do

! Ensure monotonicity by keeping a and b < 3
do i = 1, this%N-1
    if ( .not.zeroed(i) ) then
        a = this%vp(i) / Delta(i) 
        b = this%vp(i+1) / Delta(i)
        tau = 3._rprec / sqrt(a**2 + b**2)
        if (tau < 1._rprec) then
            this%vp(i) = tau * a * Delta(i)
            this%vp(i+1) = tau * b * Delta(i)
        end if
    end if
end do

! Cleanup
deallocate(Delta)
deallocate(zeroed)

end function constructor

!*******************************************************************************
subroutine interp_scalar(this, xq, vq, vqp)
!*******************************************************************************
! Perform interpolation for a single point. Uses binary_search to find the
! interval on which the sample point lies. This is a guaranteed log2(N) search
! method.
use functions, only : binary_search
implicit none

class(pchip_t) :: this
real(rprec), intent(in) :: xq
real(rprec), intent(out) :: vq
real(rprec), intent(out), optional :: vqp
integer :: i
real(rprec) :: t, h
real(rprec) :: a, b, c, d

i = binary_search(this%x, xq)
if (i == 0) then
    vq = this%v(1) + this%vp(1) * (xq - this%x(1))
    if ( present(vqp) ) vqp = this%vp(1)
else if (i == this%N) then
    vq = this%v(this%N) + this%vp(this%N) * (xq - this%x(this%N))
    if ( present(vqp) ) vqp = this%vp(this%N)
else
    h = this%x(i+1)-this%x(i)
    t = ( xq - this%x(i) ) / h
    A = (1._rprec + 2._rprec*t) * (1._rprec - t)**2
    B = t * (1._rprec - t)**2
    C = t**2 * (3._rprec - 2._rprec*t)
    D = t**2 * (t - 1._rprec)
    vq = A*this%v(i) + B*h*this%vp(i) + C*this%v(i+1) + D*h*this%vp(i+1)
    if ( present(vqp) ) then
        A = 6._rprec * t * (t - 1._rprec)
        B = 1._rprec + t * (3._rprec*t - 4._rprec)
        C = 6._rprec * t * (1._rprec - t)
        D = t * (3._rprec*t - 2._rprec)
        vqp = A*this%v(i)/h + B*this%vp(i) + C*this%v(i+1)/h + D*this%vp(i+1)
    end if
end if

end subroutine interp_scalar

!*******************************************************************************
subroutine interp_array(this, xq, vq, vqp)
!*******************************************************************************
! Perform interpolation for an array of points. This simply calls interp_scalar
! for each of the sample points.
implicit none

class(pchip_t) :: this
real(rprec), dimension(:), intent(in) :: xq
real(rprec), dimension(:), intent(out) :: vq
real(rprec), dimension(:), intent(out), optional :: vqp
integer :: i, N

N = size(xq)

if ( present(vqp) ) then
    do i = 1, N
        call this%interp(xq(i), vq(i), vqp(i))
    end do
else
    do i = 1, N
        call this%interp(xq(i), vq(i))
    end do
end if

end subroutine interp_array

end module pchip