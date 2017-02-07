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
module turbines_mpc
!*******************************************************************************
use types, only : rprec
use minimize
use wake_model
use wake_model_adjoint
use functions, only : linear_interp

type, extends(minimize_t) :: turbines_mpc_t
    type(wake_model_t) :: w                 ! wake model
    type(wake_model_t) :: iw                ! wake model initial condition
    type(wake_model_adjoint_t) :: wstar     ! adjoint wake model
    type(wake_model_adjoint_t) :: iwstar    ! adjoint wake model initial condition
    integer :: N, Nt
    logical :: isDimensionless = .false.
    real(rprec) :: cfl, dt
    real(rprec), dimension(:), allocatable :: t, Pref, Pfarm
    real(rprec), dimension(:,:), allocatable :: beta, gen_torque ! (turbine, time)
    real(rprec), dimension(:,:), allocatable :: grad_beta, grad_gen_torque      ! gradients
    real(rprec), dimension(:,:), allocatable :: fdgrad_beta, fdgrad_gen_torque ! finite difference gradients
    real(rprec) :: cost = 0._rprec
contains
    procedure, public :: initialize
    procedure, public :: makeDimensionless
    procedure, public :: makeDimensional
    procedure, public :: eval
    procedure, public :: get_control_vector
    procedure, public :: run
    procedure, public :: finite_difference_gradient
end type turbines_mpc_t

interface turbines_mpc_t
    module procedure :: constructor
end interface turbines_mpc_t

contains

!*******************************************************************************
function constructor(i_wm, i_t0, i_T, i_cfl, i_time, i_Pref) result(this)
!*******************************************************************************
implicit none
type(turbines_mpc_t) :: this
class(wake_model_t), intent(in) :: i_wm
real(rprec), dimension(:), intent(in) :: i_time, i_Pref
real(rprec), intent(in) :: i_t0, i_T, i_cfl

call this%initialize(i_wm, i_t0, i_T, i_cfl, i_time, i_Pref) 
end function constructor

!*******************************************************************************
subroutine initialize(this, i_wm, i_t0, i_T, i_cfl, i_time, i_Pref)
!*******************************************************************************
use functions, only : linear_interp
implicit none
class(turbines_mpc_t) :: this
type(wake_model_t), intent(in) :: i_wm
real(rprec), dimension(:), intent(in) :: i_time, i_Pref
real(rprec), intent(in) :: i_t0, i_T, i_cfl
integer :: i

! number of turbine rows
this%N = i_wm%N

! Wake model
this%iw = i_wm
call this%iw%makeDimensional
this%w = this%iw

! Adjoint wake model
this%iwstar = wake_model_adjoint_t(this%w%s, this%w%U_infty, this%w%Delta,     &
    this%w%k, this%w%Dia, this%w%rho, this%w%inertia, this%w%Nx,               &
    this%w%Ctp_spline, this%w%Cpp_spline)
this%wstar = this%iwstar

! Create time for the time horizon
this%cfl = i_cfl
this%dt = this%cfl * this%w%dx / this%w%U_infty
this%Nt = ceiling(i_T / this%dt)
allocate( this%t(this%Nt) )
do i = 1, this%Nt
    this%t(i) = i_t0 + this%dt * (i - 1)
end do

! allocate and set initial conditions for control variables
allocate( this%beta(this%N, this%Nt) )
allocate( this%gen_torque(this%N, this%Nt) )
allocate( this%fdgrad_beta(this%N, this%Nt) )
allocate( this%fdgrad_gen_torque(this%N, this%Nt) )
allocate( this%grad_beta(this%N, this%Nt) )
allocate( this%grad_gen_torque(this%N, this%Nt) )
this%beta(:,1) = this%iw%beta
this%gen_torque(:,1) = this%iw%gen_torque

! Interpolate the power signals and set initial condition for Pfarm
allocate( this%Pref(this%Nt) )
allocate( this%Pfarm(this%Nt) )
this%Pref = linear_interp(i_time, i_Pref, this%t)
this%Pfarm(1) = sum(this%w%Phat)
! 
! ! Allocate other variables
! allocate( this%Ctp(this%N, this%Nt)    )
! allocate( this%phi(this%N, this%Nt)    )
! allocate( this%grad(this%N, this%Nt)   )
! allocate( this%fdgrad(this%N, this%Nt) )

end subroutine initialize

!*******************************************************************************
subroutine makeDimensionless(this)
!*******************************************************************************
implicit none
class(turbines_mpc_t), intent(inout)  :: this

if (.not.this%isDimensionless) then
    this%isDimensionless = .true.
    this%t = this%t / this%w%TIME
    this%dt = this%dt / this%w%TIME
    this%Pref = this%Pref / this%w%POWER
    this%Pfarm = this%Pfarm / this%w%POWER
    this%cost = this%cost / this%w%POWER**2 / this%w%TIME
    this%gen_torque = this%gen_torque / this%w%TORQUE
    this%fdgrad_gen_torque = this%fdgrad_gen_torque / this%w%POWER**2 /        &
        this%w%TIME * this%w%TORQUE
    call this%w%makeDimensionless
    call this%wstar%makeDimensionless
    call this%iw%makeDimensionless
    call this%iwstar%makeDimensionless
end if
end subroutine makeDimensionless

!*******************************************************************************
subroutine makeDimensional(this)
!*******************************************************************************
implicit none
class(turbines_mpc_t), intent(inout)  :: this

if (this%isDimensionless) then
    this%isDimensionless = .false.
    this%t = this%t * this%w%TIME
    this%dt = this%dt * this%w%TIME
    this%Pref = this%Pref * this%w%POWER
    this%Pfarm = this%Pfarm * this%w%POWER
    this%cost = this%cost * this%w%POWER**2 * this%w%TIME
    this%gen_torque = this%gen_torque * this%w%TORQUE
    this%fdgrad_gen_torque = this%fdgrad_gen_torque * this%w%POWER**2 *        &
        this%w%TIME / this%w%TORQUE
    call this%w%makeDimensional
    call this%wstar%makeDimensional
    call this%iw%makeDimensional
    call this%iwstar%makeDimensional
end if
end subroutine makeDimensional

!*******************************************************************************
subroutine run(this)
!*******************************************************************************
implicit none
class(turbines_mpc_t), intent(inout) :: this
integer :: i, k
real(rprec), dimension(:,:,:), allocatable :: fstar
real(rprec), dimension(:,:), allocatable :: Adu, Aw, Bj, Bdu, Bw, Gdu, Gw

! allocate adjoint forcing terms
allocate(fstar(this%Nt,this%N,this%w%Nx))
allocate(Adu(this%Nt,this%N))
allocate(Aw(this%Nt,this%N))
allocate(Bj(this%Nt,this%N))
allocate(Bdu(this%Nt,this%N))
allocate(Bw(this%Nt,this%N))
allocate(Gdu(this%Nt,this%N))
allocate(Gw(this%Nt,this%N))
fstar =  0._rprec
Adu = 0._rprec
Aw = 0._rprec
Bj = 0._rprec
Bdu = 0._rprec
Bw = 0._rprec
Gdu = 0._rprec
Gw = 0._rprec

! reset costs and gradients
this%cost = 0._rprec
this%grad_beta = 0._rprec
this%grad_gen_torque = 0._rprec

! Run forward in time
this%w = this%iw
call this%w%adjoint_values(this%Pref(1), fstar(1,:,:), Adu(1,:), Aw(1,:),  &
    Bj(1,:), Bdu(1,:), Bw(1,:), Gdu(1,:), Gw(1,:))
do k = 2, this%Nt
    call this%w%advance(this%beta(:,k), this%gen_torque(:,k), this%dt)
    this%Pfarm(k) = sum(this%w%Phat)
    call this%w%adjoint_values(this%Pref(k), fstar(k,:,:), Adu(k,:), Aw(k,:),  &
        Bj(k,:), Bdu(k,:), Bw(k,:), Gdu(k,:), Gw(k,:))
    ! calculate contribution of cost function to gradient
    this%grad_gen_torque(:,k) = 2._rprec * (this%Pfarm(k) - this%Pref(k))      &
        * this%w%omega * this%dt
end do
this%cost = sum( (this%Pref-this%Pfarm)**2 ) * this%dt

! Run backwards in time
this%wstar = this%iwstar
do k = this%Nt-1, 1, -1
    call this%wstar%retract(fstar(k+1,:,:), Adu(k+1,:), Aw(k+1,:), Bj(k+1,:),  &
        Bdu(k+1,:), Bw(k+1,:), this%dt)
    this%grad_gen_torque(:,k) = this%grad_gen_torque(:,k)                      &
        + this%wstar%omega_star / this%wstar%inertia * this%dt
    this%grad_beta(:,k) = Gw(k,:) * this%wstar%omega_star * this%dt
    do i = 1, this%N
        this%grad_beta(i,k) = this%grad_beta(i,k) + Gdu(k,i)  * this%dt        &
            * sum(this%wstar%du_star(i,:) * this%wstar%f(i,:)) * this%wstar%dx
    end do
end do

! deallocate adjoint forcing terms
deallocate(fstar)
deallocate(Adu)
deallocate(Aw)
deallocate(Bj)
deallocate(Bdu)
deallocate(Bw)
deallocate(Gdu)
deallocate(Gw)

end subroutine run

!*******************************************************************************
subroutine eval(this, x, f, g)
!*******************************************************************************
implicit none
class(turbines_mpc_t), intent(inout) :: this
real(rprec), dimension(:), intent(in) :: x
real(rprec), intent(inout) :: f
real(rprec), dimension(:), intent(inout) :: g
integer :: k, istart, istop, iskip

! Place x in control variables
iskip = (this%Nt-1) * this%N
do k = 1, this%Nt-1
    istart = (k-1)*this%N+1
    istop = this%N*k
    this%beta(:,k) = x(istart:istop)
    this%gen_torque(:,k) = x(istart+iskip:istop+iskip)
end do

! Run model
call this%run

! Return cost function
f = this%cost

! Return gradient as vector
g = 0._rprec
do k = 1, this%Nt-1
istart = (k-1)*this%N+1
    istart = (k-1)*this%N+1
    istop = this%N*k
    g(istart:istop) = this%grad_beta(:,k)
    g(istart+iskip:istop+iskip) = this%grad_gen_torque(:,k)
end do
end subroutine eval

!*******************************************************************************
function get_control_vector(this) result(x)
!*******************************************************************************
implicit none
class(turbines_mpc_t) :: this
real(rprec), dimension(:), allocatable :: x
integer :: k, istart, istop, iskip

allocate( x(2*(size(this%beta) - this%N)) )
x = -1000

iskip = (this%Nt-1) * this%N
do k = 1, this%Nt-1
    istart = (k - 1) * this%N + 1
    istop = this%N * k
    write(*,*) istart, istop    
    x(istart:istop) = this%beta(:,k)
    x(istart+iskip:istop+iskip) = this%gen_torque(:,k)
end do
    
end function get_control_vector

!*******************************************************************************
subroutine finite_difference_gradient(this)
!*******************************************************************************
implicit none
class(turbines_mpc_t), intent(inout) :: this
type(turbines_mpc_t) :: mf
integer :: n, k
real(rprec) :: dphi

! finite difference
dphi = sqrt( epsilon( this%beta(1,1) ) )

! calculate gradient for beta
this%fdgrad_beta = 0._rprec
do n = 1, this%N
    do k = 1, this%Nt
        mf = this
        mf%beta(n,k) = mf%beta(n,k) + dphi
        call mf%run
        this%fdgrad_beta(n,k) = (mf%cost - this%cost) / dphi
    end do
end do

! calculate gradient for gen_torque
this%fdgrad_gen_torque = 0._rprec
do n = 1, this%N
    do k = 1, this%Nt
        mf = this
        mf%gen_torque(n,k) = mf%gen_torque(n,k) + dphi
        call mf%run
        this%fdgrad_gen_torque(n,k) = (mf%cost - this%cost) / dphi
    end do
end do

end subroutine finite_difference_gradient

end module turbines_mpc