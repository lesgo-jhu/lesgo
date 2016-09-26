!!
!!  Copyright (C) 2016  Johns Hopkins University
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
module cubic_spline
!*******************************************************************************
! The cubic_spline_t class performs cubic spline interpolation for a 1D function
! v(x). The evaluated data is passed on construction. Interpolation is evaluated
! for real sample points xq passed as values or arrays. The first derivatives 
! can also be returned as optional arguments.
!
! The resulting interpolation has continuous first and second derivatives. 
! Linear extrapolation is used for evaluation points xq outside of the interval 
! of x's passed at construction. The second derivative vanishes as x(1) and x(N)
use types, only : rprec
use messages
implicit none

private
public cubic_spline_t

type :: cubic_spline_t
    real(rprec), dimension(:), allocatable :: x, v, vpp
    integer :: N
contains
    procedure, private :: interp_scalar, interp_array
    generic, public :: interp => interp_scalar, interp_array
end type cubic_spline_t

interface cubic_spline_t
    module procedure constructor
end interface cubic_spline_t

contains

!*******************************************************************************
function constructor(i_x, i_v) result(this)
!*******************************************************************************
! Constructor for cubic_spline_t. Takes points v(x) that are used for the 
! interpolation. This function also evaluates the second derivative as these x's
! using the tridiagonal matrix algorithm. The second derivative at x(1) and x(N)
! are set to 0.
use tridiagonal
implicit none
type(cubic_spline_t) :: this
real(rprec), dimension(:), intent(in) :: i_x, i_v
real(rprec), dimension(:), allocatable :: a, b, c, d
type(tridiagonal_t) :: M
integer :: i, j

! Set the size of the matrix
this%N = size(i_x)

! Check that all input arguments are the same size
if ( size(i_v) /= this%N ) then
    call error('cubic_spline_t/constructor','x and v must be the same size')
end if

! Check that i_x is sorted
do i = 2, this%N
    if ( i_x(i) < i_x(i-1) ) then
        call error('cubic_spline_t/constructor', 'x must be increasing')
    end if
end do

! Allocate and assign variables
allocate( this%x(this%N) )
allocate( this%v(this%N) )
allocate( this%vpp(this%N) )
this%x = i_x
this%v = i_v
this%vpp = 0._rprec

! Allocate tridiagonal matrix arrays. 
! Dimension is 2 less because vpp(0)=0 and vpp(N)=0.
allocate( a(this%N-2) )
allocate( b(this%N-2) )
allocate( c(this%N-2) )
allocate( d(this%N-2) )

! Calculate the tridiagonal matrix arrays. 
do i = 1, this%N-2
    j = i+1
    a(i) = (this%x(j) - this%x(j-1)) / 6._rprec
    b(i) = (this%x(j+1) - this%x(j-1)) / 3._rprec
    c(i) = (this%x(j+1) - this%x(j)) / 6._rprec
    d(i) = (this%v(j+1) - this%v(j)) / (this%x(j+1) - this%x(j))               &
         - (this%v(j) - this%v(j-1)) / (this%x(j) - this%x(j-1)) 
end do

! Solve the system and place answer into vpp
M = tridiagonal_t(a, b, c)
call M%solve(d)
do i = 2, this%N-1
    this%vpp(i) = d(i-1)
end do

! Cleanup
deallocate(a)
deallocate(b)
deallocate(c)
deallocate(d)

end function constructor

!*******************************************************************************
subroutine interp_scalar(this, xq, vq, vqp)
!*******************************************************************************
! Perform interpolation for a single point. Uses binary_search to find the 
! interval on which the sample point lies. This is a guarnateed log2(N) search
! method.
use functions, only : binary_search
implicit none
class(cubic_spline_t) :: this
real(rprec), intent(in) :: xq
real(rprec), intent(out) :: vq
real(rprec), intent(out), optional :: vqp
real(rprec) :: vqpi
integer :: i
real(rprec) :: A, B, C, D

i = binary_search(this%x, xq)
if (i == 0) then
    vqpi = (this%v(2) - this%v(1)) / (this%x(2) - this%x(1))                   &
            - (this%x(2) - this%x(1)) * this%vpp(1)/ 3.                        &
            - (this%x(2) - this%x(1)) * this%vpp(2)/ 6.
    vq = this%v(1) + vqpi * (xq - this%x(1))
    if ( present(vqp) ) vqp = vqpi
else if (i == this%N) then
    vqpi = (this%v(this%N) - this%v(this%N-1))                                 &
        / (this%x(this%N) - this%x(this%N-1))                                  &
        + (this%x(this%N) - this%x(this%N-1)) * this%vpp(this%N-1)/ 6.         &
        + (this%x(this%N) - this%x(this%N-1)) * this%vpp(this%N)/ 3.
    vq = this%v(this%N) + vqpi * (xq - this%x(this%N))
    if ( present(vqp) ) vqp = vqpi
else
    A = (this%x(i+1) - xq) / (this%x(i+1) - this%x(i))
    B = (1-A)
    C = (A**3 - A) * (this%x(i+1) - this%x(i))**2 / 6._rprec
    D = (B**3 - B) * (this%x(i+1) - this%x(i))**2 / 6._rprec
    vq = A*this%v(i) + B*this%v(i+1) + C*this%vpp(i) + D*this%vpp(i+1)
    if ( present(vqp) ) then
        vqp = (this%v(i+1) - this%v(i)) / (this%x(i+1) - this%x(i))            &
            - (3.*A*A - 1.) * (this%x(i+1) - this%x(i)) * this%vpp(i)/ 6.      &
            + (3.*B*B - 1.) * (this%x(i+1) - this%x(i)) * this%vpp(i+1)/ 6.
    end if
end if

end subroutine interp_scalar

!*******************************************************************************
subroutine interp_array(this, xq, vq, vqp)
!*******************************************************************************
! Perform interpolation for an array of points. This simply calls interp_scalar
! for each of the sample points.
implicit none
class(cubic_spline_t) :: this
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

end module cubic_spline