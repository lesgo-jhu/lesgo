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
module tridiagonal_m
!*******************************************************************************
! The tridiagonal_t class solves a tridiagonal matrix equation Ax = d, where A
! is an NxN matrix and x and d are vectors of length N. The matrix M is composed
! of three vectors of length N:
!
!       a:  subdiagonal     (The first element is not used)
!       b:  diagonal
!       c:  superdiagonal   (The last element is not used)
!
use types, only : rprec, cprec
implicit none

private
public tridiagonal_t

type :: tridiagonal_t
    real(rprec), dimension(:), allocatable :: a, b, c, p, m
    integer :: N
contains
    procedure, public :: solve => solve_real, solve_cmplx
end type tridiagonal_t

interface tridiagonal_t
    module procedure :: constructor
end interface tridiagonal_t

contains

!*******************************************************************************
function constructor(i_a, i_b, i_c) result(this)
!*******************************************************************************
! Constructor for tridiagonal_t that takes the three diagonal vectors are
! input arguments.
!
real(rprec), dimension(:), intent(in) :: i_a, i_b, i_c
type(tridiagonal_t) :: this
integer :: i

! Set the size of the matrix
this%N = size(i_a)

! Check that all input arguments are the same size
if ( size(i_b) /= this%N .or. size(i_c) /= this%N ) then
    write(*,*) "ERROR: tridiagonal_t%constructor: " //                         &
        "a, b, and c must be the same size"
    stop 9
end if

! Allocate arrays
allocate( this%a(this%N) )
allocate( this%b(this%N) )
allocate( this%c(this%N) )
allocate( this%p(this%N) )
allocate( this%m(this%N) )

! Place input arguments into object
this%a = i_a
this%b = i_b
this%c = i_c

! Check diagonal for zeros
do i = 1, this%N
    if (this%b(i) == 0._rprec) then
        write(*,*) "ERROR: tridiagonal_t%constructor: " //                    &
            "found zero element along diagonal on row ",  i
        stop 9
    end if
end do

! Perform forward matrix evaluations
this%p(1) = this%b(1)
this%m(1) = 0._rprec
do i = 2, this%N
    this%m(i) = this%a(i) / this%p(i-1)
    this%p(i) = this%b(i) - this%m(i) * this%c(i-1)
    if (this%p(i) == 0._rprec) then
        write(*,*) "ERROR: tridiagonal_t%constructor: " //                     &
            "found zero pivot on row ", i
        stop 9
    end if
end do

end function constructor

!*******************************************************************************
subroutine solve_real(this, d)
!*******************************************************************************
! Solves the matrix equation when given the right hand side d. The answer is
! returned in the same input array.
!
class(tridiagonal_t) :: this
real(rprec), dimension(:), intent(inout) :: d
integer :: i

! Check input argument size
if ( size(d) /= this%N ) then
    write(*,*) "ERROR: tridiagonal_t%csolve: d must the same size as the matrix"
    stop 9
end if

! Forward sweep
do i = 2, this%N
    d(i) = d(i) - this%m(i)*d(i-1)
end do

! Backward sweep
d(this%N) = d(this%N) / this%p(this%N)
do i = this%N-1, 1, -1
    d(i) = (d(i) - this%c(i) * d(i+1)) / this%p(i);
end do

end subroutine solve_real

!*******************************************************************************
subroutine solve_cmplx(this, d)
!*******************************************************************************
! Solves the matrix equation when given the right hand side d. The answer is
! returned in the same input array.
!
class(tridiagonal_t) :: this
real(rprec), dimension(:), intent(inout) :: d
integer :: i

! Check input argument size
if ( size(d) /= this%N ) then
    write(*,*) "ERROR: tridiagonal_t%csolve: d must the same size as the matrix"
    stop 9
end if

! Forward sweep
do i = 2, this%N
    d(i) = d(i) - this%m(i)*d(i-1)
end do

! Backward sweep
d(this%N) = d(this%N) / this%p(this%N)
do i = this%N-1, 1, -1
    d(i) = (d(i) - this%c(i) * d(i+1)) / this%p(i);
end do

end subroutine solve_cmplx

end module tridiagonal_m
