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
module minimize_m
!*******************************************************************************
use types, only : rprec
use messages
implicit none

private
public minimize_t

type :: minimize_t
    ! pointer to function to minimize. Not used if extending this base class
    procedure(minimize_function), pointer, nopass :: fun => NULL()
contains
    ! This procedure is used by the minimizer to find the minimum
    procedure, public :: eval
end type minimize_t

interface minimize_t
    module procedure :: constructor
end interface minimize_t

! Interface for the scalar-valued function to be minimize_t. This interface
! must evaluate the function at x and return the function value f and gradient g
! at this location.
abstract interface
    subroutine minimize_function(x, f, g)
        use types, only : rprec
        implicit none
        real(rprec), dimension(:), intent(in) :: x      ! Point to evaluate
        real(rprec), intent(inout) :: f                 ! Function value (scalar)
        real(rprec), dimension(:), intent(inout) :: g   ! Function gradient
    end subroutine minimize_function
end interface

contains

!*******************************************************************************
function constructor(i_fun) result(this)
!*******************************************************************************
! Constructor for minimize_t class. This constructor takes a procedure pointer
! matching the minimize_function interface. This constructor can be overloaded
! for a class that extends this base class.
implicit none
procedure(minimize_function) :: i_fun
type(minimize_t) :: this

this%fun => i_fun

end function  constructor

!*******************************************************************************
subroutine eval(this, x, f, g)
!*******************************************************************************
! Evaluates the function pointer. Overload this procedure if extending this base
! class
implicit none
class(minimize_t), intent(inout) :: this
real(rprec), dimension(:), intent(in) :: x
real(rprec), intent(inout) :: f
real(rprec), dimension(:), intent(inout) :: g

call this%fun(x, f, g)

end subroutine eval

end module minimize_m
