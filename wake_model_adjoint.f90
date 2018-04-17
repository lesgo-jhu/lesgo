!!
!!  Copyright (C) 2018  Johns Hopkins University
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
module wake_model_adjoint
!*******************************************************************************
use types, only : rprec
use util,  only : logistic, softplus, gaussian
use wake_model_base
use messages
implicit none

private
public wake_model_adjoint_t

type, extends(wake_model_base_t) :: wake_model_adjoint_t
    ! adjoint velocity deficit (turbine, space)
    real(rprec), dimension(:,:), allocatable :: dustar
contains
    procedure, public  :: initialize_val
    procedure, public  :: makeDimensionless
    procedure, public  :: makeDimensional
    procedure, public  :: retract
    procedure, private :: rhs

end type wake_model_adjoint_t

interface wake_model_adjoint_t
    module procedure :: constructor_val
end interface wake_model_adjoint_t

contains

!*******************************************************************************
function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx) result(this)
!*******************************************************************************
! Constructor for wake model with values given
implicit none
type(wake_model_adjoint_t)                :: this
real(rprec), intent(in)               :: i_U_infty, i_Delta, i_Dia
real(rprec), dimension(:), intent(in) :: i_s, i_k
integer, intent(in)                   :: i_Nx

call this%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx)

end function constructor_val

!*******************************************************************************
subroutine initialize_val(this, i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx)
!*******************************************************************************
use wake_model_base
implicit none
class(wake_model_adjoint_t), intent(inout) :: this
real(rprec), intent(in)                :: i_U_infty, i_Delta, i_Dia
real(rprec), dimension(:), intent(in)  :: i_s, i_k
integer, intent(in)                    :: i_Nx

! Call base class initializer
call this%wake_model_base_t%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx)

allocate( this%dustar(this%N, this%Nx) )

this%dustar(:,:) = 0.d0

end subroutine initialize_val

!*******************************************************************************
subroutine makeDimensionless(this)
!*******************************************************************************
implicit none
class(wake_model_adjoint_t), intent(inout)  :: this

if (.not.this%isDimensionless) then
    call this%wake_model_base_t%makeDimensionless
    this%dustar = this%dustar / this%VELOCITY
end if

end subroutine makeDimensionless

!*******************************************************************************
subroutine makeDimensional(this)
!*******************************************************************************
implicit none
class(wake_model_adjoint_t), intent(inout)  :: this

if (this%isDimensionless) then
    call this%wake_model_base_t%makeDimensional
    this%dustar = this%dustar * this%VELOCITY
end if

end subroutine makeDimensional

!*******************************************************************************
subroutine retract(this, fstar, dt, g)
!*******************************************************************************
implicit none
class(wake_model_adjoint_t), intent(inout)   :: this
real(rprec), dimension(:,:), intent(in)  :: fstar
real(rprec), intent(in)                  :: dt
real(rprec), dimension(:), intent(inout) :: g
integer                                  :: i

do i = 1, this%N
    this%dustar(i,:) = this%dustar(i,:) - this%rhs(this%dustar(i,:), fstar(i,:), i) * dt
    g(i) = sum(this%G(i,:) * this%dustar(i,:) / this%d(i,:) / this%d(i,:)) * this%dx
end do

end subroutine retract

!*******************************************************************************
function rhs(this, du, f, i) result(ddudt)
!*******************************************************************************
! Evaluates RHS of wake equation using 3rd-order biased downwind differencing
use util, only : ddx_downwind1
implicit none
class(wake_model_adjoint_t), intent(in)    :: this
real(rprec), dimension(:), intent(in)  :: f, du
real(rprec), dimension(:), allocatable :: ddudt, ddudx
integer, intent(in)                    :: i

allocate(ddudt(this%Nx))
allocate(ddudx(this%Nx))

ddudx = ddx_downwind1(du, this%dx)
ddudt = -this%U_infty * ddudx + this%w(i,:) * du - f

end function rhs

end module wake_model_adjoint
