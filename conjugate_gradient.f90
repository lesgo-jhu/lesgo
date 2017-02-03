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
module conjugate_gradient
!*******************************************************************************
use types, only : rprec
use line_search
use minimize
use messages
implicit none

private
public conjugate_gradient_t

type :: conjugate_gradient_t
    class(minimize_t), pointer :: mini => NULL()
    real(rprec) :: eps = 1E-10      ! A small number
    integer :: maxIter = 10000      ! maximum number of CG iterations
    real(rprec) :: tol = 1E-6       ! convergence level
    real(rprec) :: f = 0.0          ! current function evaluation
    real(rprec) :: fp = 0.0         ! previous function evaluation
    integer :: fnev = 0             ! number of function evaluations
    real(rprec) :: gamma = 0.0      ! conjugate direction step
    type(line_search_t) :: ls       ! line search class
    real(rprec), dimension(:), allocatable :: gd    ! conjugate direction
    real(rprec), dimension(:), allocatable :: x, xp ! current and previous location
    real(rprec), dimension(:), allocatable :: g, gp ! current and previous gradients
    real(rprec), dimension(:), allocatable :: h, hp ! current and previous search direction
    real(rprec) :: lb = -10000000                   ! -Infinity
    real(rprec) :: ub = 10000000                    ! Infinity
contains
   procedure, public :: minimize
   procedure, private :: evaluate_gamma
end type conjugate_gradient_t

interface conjugate_gradient_t
    module procedure :: constructor
end interface conjugate_gradient_t

contains

!*******************************************************************************
function constructor(i_mini, i_maxIter, i_lb, i_ub, i_tol) result(this)
!*******************************************************************************
implicit none
type(conjugate_gradient_t) :: this
class(minimize_t), target :: i_mini
integer, intent(in), optional :: i_maxIter
real(rprec), intent(in), optional :: i_lb, i_ub, i_tol

! Assign input arguments
if ( present(i_maxIter) )   this%maxIter = i_maxIter
if ( present(i_tol) )       this%tol     = i_tol
if ( present(i_lb) )        this%lb      = i_lb
if ( present(i_ub) )        this%ub      = i_ub
this%mini => i_mini
this%ls = line_search_t(i_mini)

end function constructor

!*******************************************************************************
subroutine evaluate_gamma(this)
!*******************************************************************************
implicit none
class(conjugate_gradient_t), intent(inout) :: this

this%gd = this%g - this%gp
this%gamma = sum(this%g * this%gd) / sum(this%gp * this%gp)

end subroutine evaluate_gamma

!*******************************************************************************
subroutine minimize(this, i_x, o_x)
!*******************************************************************************
implicit none
class(conjugate_gradient_t), intent(inout) :: this
real(rprec), dimension(:), intent(in) :: i_x
real(rprec), dimension(:), intent(out), optional :: o_x
real(rprec) :: d, delta_f, stp
integer :: i, j
integer :: dummy = 0

! Allocate arrays
allocate(this%xp(size(i_x)))
allocate(this%gp(size(i_x)))
allocate(this%hp(size(i_x)))
allocate(this%x(size(i_x)))
allocate(this%g(size(i_x)))
allocate(this%h(size(i_x)))
allocate(this%gd(size(i_x)))

! Evaluate gradient and function at starting position
this%xp = i_x
call this%mini%eval(this%xp, this%fp, this%gp)
this%fnev = this%fnev + 1
this%hp = -this%gp

! Assign other arrays
this%f  = this%fp
this%g  = this%gp
this%h  = this%hp
this%x  = this%xp
this%gd = this%g

d = sum(this%g * this%h)
delta_f = -0.5 * d

do i = 1, this%maxIter
    ! Compute derivative with respect to stp at origin
    d = sum(this%g * this%h)

    ! If derivative is positive, reset at the gradient
    if (d > 0) then
        this%h = -this%g
        d = sum(this%g * this%h)
        delta_f = -0.5 * d
    end if

    ! Search in the direction hp
    stp = min(1.0, -2.0 * delta_f / d)
    call this%ls%search(this%x, this%f, this%g, this%h, stp, dummy)
    this%fnev = this%fnev + dummy

    ! Safeguard step against bounds
    do j = 1, size(this%x)
        this%x(j) = min( max(this%x(j), this%lb), this%ub)
    end do
    call this%mini%eval(this%x, this%f, this%g)

    ! Check for convergence
    if ( 2.0 * abs(this%fp - this%f) <= this%tol * ( abs(this%fp) +            &
    abs(this%f) + this%eps ) .or. sum(this%g*this%g) == 0 ) then
        ! Set output if present
        if ( present(o_x) ) o_x = this%x

        ! Evaluate minimization at current point
        call this%mini%eval(this%x, this%f, this%g)

        ! Print result
        write(*,*) 'Conjugate gradient terminated after ', i,                  &
            'iterations. Minimum f = ', this%f

        return
    else if (this%f > this%fp) then
        ! Reset to previous point
        this%x = this%xp;
        this%f = this%fp;
        this%g = this%gp;

        ! Set output if present
        if ( present(o_x) ) o_x = this%x

        ! Evaluate minimization at current point
        call this%mini%eval(this%x, this%f, this%g)

        ! Print result
        write(*,*) 'Conjugate gradient terminated after ', i,                  &
            'iterations. Minimum f = ', this%f
        return
    end if

    ! Compute next search direction
    call this%evaluate_gamma
    this%h = -this%g + this%gamma * this%hp

    ! Swap arrays
    delta_f = max(this%fp - this%f, 10 * epsilon(this%f))

    this%gp = this%g
    this%hp = this%h
    this%xp = this%x
    this%fp = this%f
end do

! Evaluate minimization at current point
call this%mini%eval(this%x, this%f, this%g)

! Set output if present
if ( present(o_x) ) o_x = this%x

! Print result
write(*,*) 'Conjugate gradient terminated after ', this%maxIter,               &
    'iterations. Minimum f = ', this%f

end subroutine minimize

end module conjugate_gradient
