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
module cubic_spline
!*******************************************************************************
! The cubic_spline_t class performs cubic spline interpolation for a 1D function
! v(x). The evaluated data is passed on construction. Interpolation is evaluated
! for real sample points xq passed as values or arrays. The first derivatives
! can also be returned as optional arguments.
!
! The resulting interpolation has continuous first and second derivatives.
! Linear extrapolation is used for evaluation points xq outside of the interval
! of x's passed at construction.
!
! Three boundary conditions are available:
!   Natural:    Second derivatives are at x(1) or x(N)
!   Clamped:    First derivatives are specified at x(1) or x(N)
!   Not-a-knot: A single spline is used on the interval between x(1) and x(3)
!               or x(N-2) and x(N). This can alternatively be viewed as
!               requiring continuity of the third derivative at x(2) or x(N-1)
!
use types, only : rprec
use messages
use string_util, only : uppercase
implicit none

private
public cubic_spline_t

type :: cubic_spline_t
    real(rprec), dimension(:), allocatable :: x, v, vpp
    integer :: N
    character(:), allocatable :: low_bc, high_bc
    real(rprec) :: low_f = 0._rprec, high_f = 0._rprec
    character(14) :: type_name = "cubic_spline_t"
contains
    procedure, private :: interp_scalar
    procedure, private :: interp_array
    generic, public :: interp => interp_scalar, interp_array
end type cubic_spline_t

interface cubic_spline_t
    module procedure constructor
end interface cubic_spline_t

contains

!*******************************************************************************
function constructor(x, v, low_bc, high_bc, low_f, high_f) result(this)
!*******************************************************************************
! Constructor for cubic_spline_t. Takes points v(x) that are used for the
! interpolation. This function also evaluates the second derivative as these x's
! using the tridiagonal matrix algorithm.
!
! By default, the constructor using natural boundary conditions with vanishing
! second derivatives. Other boundary conditions can be specified by supplying
! the optional arguments low_bc or high_bc. The values of the specified
! derivatives can be passed sing the optional arguments low_f and high_f.
!
use tridiagonal
implicit none

type(cubic_spline_t) :: this
real(rprec), dimension(:), intent(in) :: x, v
character(*), intent(in), optional :: low_bc, high_bc
real(rprec), intent(in), optional :: low_f, high_f
real(rprec), dimension(:), allocatable :: a, b, c, d
real(rprec) :: aa, bb
integer :: Nm, offset
type(tridiagonal_t) :: M
integer :: i, j
character(14) :: proc_name = "constructor"

! set optional arguments
this%low_bc = "natural"
this%high_bc = "natural"
if (present(low_bc)) this%low_bc = low_bc
if (present(high_bc)) this%high_bc = high_bc
if (present(low_f)) this%low_f = low_f
if (present(high_f)) this%high_f = high_f

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
allocate( this%vpp(this%N) )
this%x = x
this%v = v
this%vpp = 0._rprec

! Natural and clamped boundary conditions
Nm = this%N
offset = 0

! Not-a-knot BCs have one less dimension per boundary condition
if (uppercase(this%low_bc) == 'NOT-A-KNOT') then
    Nm = Nm - 1
    offset = 1
end if
if  (uppercase(this%high_bc) == 'NOT-A-KNOT') then
    Nm = Nm - 1
end if

! Allocate tridiagonal matrix arrays.
allocate( a(Nm) )
allocate( b(Nm) )
allocate( c(Nm) )
allocate( d(Nm) )

! Calculate interior of the tridiagonal matrix arrays
do i = 2, Nm-1
    j = i + offset
    a(i) = (this%x(j) - this%x(j-1)) / 6._rprec
    b(i) = (this%x(j+1) - this%x(j-1)) / 3._rprec
    c(i) = (this%x(j+1) - this%x(j)) / 6._rprec
    d(i) = (this%v(j+1) - this%v(j)) / (this%x(j+1) - this%x(j))               &
         - (this%v(j) - this%v(j-1)) / (this%x(j) - this%x(j-1))
end do

! Apply low boundary conditions
select case (uppercase(this%low_bc))
    case('NATURAL')
        b(1) = 1._rprec
        c(1) = 0._rprec
        d(1) = this%low_f
    case('CLAMPED')
        b(1) = (this%x(2) - this%x(1)) / 3._rprec
        c(1) = (this%x(2) - this%x(1)) / 6._rprec
        d(1) = (this%v(2) - this%v(1)) / (this%x(2) - this%x(1)) - this%low_f
    case('NOT-A-KNOT')
        a(2) = (this%x(3) - this%x(1)) / 6._rprec
        b(2) = (this%x(4) - this%x(1)) / 3._rprec
        c(2) = (this%x(4) - this%x(3)) / 6._rprec
        d(2) = (this%v(4) - this%v(3)) / (this%x(4) - this%x(3))               &
             - (this%v(3) - this%v(1)) / (this%x(3) - this%x(1))
        aa = (this%x(3) - this%x(2)) / (this%x(3) - this%x(1))
        bb = (1 - aa)
        b(1) = (aa**3 - aa) * (this%x(3) - this%x(1))**2 / 6._rprec
        c(1) = (bb**3 - bb) * (this%x(3) - this%x(1))**2 / 6._rprec
        d(1) = this%v(2) - aa*this%v(1) - bb*this%v(3)
    case default
        call error(this%type_name // '.' // proc_name,                         &
            'Invalid low boundary condition type ' // this%low_bc)
end select

! Apply high boundary conditions
select case (uppercase(this%high_bc))
    case('NATURAL')
        b(Nm) = 1._rprec
        a(Nm) = 0._rprec
        d(Nm) = this%high_f
    case('CLAMPED')
        b(Nm) = (this%x(this%N) - this%x(this%N-1)) / 3._rprec
        a(Nm) = (this%x(this%N) - this%x(this%N-1)) / 6._rprec
        d(Nm) = this%high_f - (this%v(this%N) - this%v(this%N-1)) /                 &
            (this%x(this%N) - this%x(this%N-1))
    case('NOT-A-KNOT')
        a(Nm-1) = (this%x(this%N-2) - this%x(this%N-3)) / 6._rprec
        b(Nm-1) = (this%x(this%N) - this%x(this%N-3)) / 3._rprec
        c(Nm-1) = (this%x(this%N) - this%x(this%N-2)) / 6._rprec
        d(Nm-1) = (this%v(this%N) - this%v(this%N-2))                          &
                    / (this%x(this%N) - this%x(this%N-2))                      &
                - (this%v(this%N-2) - this%v(this%N-3))                        &
                    / (this%x(this%N-2) - this%x(this%N-3))
        aa = (this%x(this%N) - this%x(this%N-1))                               &
            / (this%x(this%N) - this%x(this%N-2))
        bb = (1 - aa)
        b(Nm) = (aa**3 - aa) * (this%x(this%N) - this%x(this%N-2))**2 / 6._rprec
        a(Nm) = (bb**3 - bb) * (this%x(this%N) - this%x(this%N-2))**2 / 6._rprec
        d(Nm) = this%v(this%N-1) - aa*this%v(this%N-2) - bb*this%v(this%N)
    case default
        call error(this%type_name // '.' // proc_name,                         &
            'Invalid high boundary condition type ' // this%high_bc)
end select

! Solve the system and place answer into vpp
M = tridiagonal_t(a, b, c)
call M%solve(d)
do i = 1, Nm
    j = i + offset
    this%vpp(j) = d(i)
end do

! fix lower boundary condition if not-a-knot
if (uppercase(this%low_bc) == 'NOT-A-KNOT') then
    this%vpp(1) = d(1)
    aa = (this%x(3) - this%x(2)) / (this%x(3) - this%x(1))
    bb = (1 - aa)
    this%vpp(2) = aa * this%vpp(1) + bb*this%vpp(3)
end if

! fix upper boundary condition if not-a-knot
if (uppercase(this%high_bc) == 'NOT-A-KNOT') then
    this%vpp(this%N) = d(Nm)
    aa = (this%x(this%N) - this%x(this%N-1))                                   &
        / (this%x(this%N) - this%x(this%N-2))
    bb = (1 - aa)
    this%vpp(this%N-1) = aa * this%vpp(this%N-2) + bb*this%vpp(this%N)
end if


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
! interval on which the sample point lies. This is a guaranteed log2(N) search
! method.
!
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
!
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
