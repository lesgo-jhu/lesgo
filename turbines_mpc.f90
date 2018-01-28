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
use minimize_m
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
    ! Control variables (turbine, time)
    real(rprec), dimension(:,:), allocatable :: beta, torque_gain
    ! gradients
    real(rprec), dimension(:,:), allocatable :: grad_beta, grad_torque_gain
    ! finite difference gradients
    real(rprec), dimension(:,:), allocatable :: fdgrad_beta, fdgrad_torque_gain
    real(rprec) :: cost = 0._rprec
    ! scaling constants for gradient
    real(rprec) :: Ca = 1._rprec
    real(rprec) :: Cb = 1._rprec
    ! Constraints on controls
    real(rprec) :: beta_min = -45._rprec
    real(rprec) :: beta_max = 45._rprec
    real(rprec) :: torque_gain_min = 0._rprec
    real(rprec) :: torque_gain_max = 1000000000000._rprec
    ! Constraints on rotational speed (penalized in cost function)
    real(rprec) :: omega_min = 0._rprec
    real(rprec) :: omega_max = 1000000000000._rprec
    real(rprec) :: eta1 = 1000._rprec
    ! Penalty away from optimal TSR and pitch angle
    real(rprec) :: eta2 = 0.000001_rprec
    real(rprec) :: beta_star = 0._rprec
    real(rprec) :: eta3 = 0.000001_rprec
    real(rprec) :: lambda_prime_star = 11._rprec
contains
    procedure, public :: initialize
    procedure, public :: makeDimensionless
    procedure, public :: makeDimensional
    procedure, public :: eval
    procedure, public :: get_control_vector
    procedure, public :: get_lower_bound
    procedure, public :: get_upper_bound
    procedure, public :: run
    procedure, public :: finite_difference_gradient
    procedure, public :: rescale_gradient
end type turbines_mpc_t

interface turbines_mpc_t
    module procedure :: constructor
end interface turbines_mpc_t

contains

!*******************************************************************************
function constructor(i_wm, i_t0, i_T, i_cfl, i_time, i_Pref, i_beta_star,      &
    i_lambda_prime_star, i_omega_min, i_omega_max) result(this)
!*******************************************************************************
implicit none
type(turbines_mpc_t) :: this
class(wake_model_t), intent(in) :: i_wm
real(rprec), dimension(:), intent(in) :: i_time, i_Pref
real(rprec), intent(in) :: i_t0, i_T, i_cfl, i_beta_star, i_lambda_prime_star
real(rprec), intent(in), optional :: i_omega_min, i_omega_max

if (present(i_omega_max)) then
    call this%initialize(i_wm, i_t0, i_T, i_cfl, i_time, i_Pref, i_beta_star,  &
        i_lambda_prime_star, i_omega_min, i_omega_max)
elseif (present(i_omega_min)) then
    call this%initialize(i_wm, i_t0, i_T, i_cfl, i_time, i_Pref, i_beta_star,  &
        i_lambda_prime_star, i_omega_min)
else
    call this%initialize(i_wm, i_t0, i_T, i_cfl, i_time, i_Pref, i_beta_star,  &
        i_lambda_prime_star)
end if

end function constructor

!*******************************************************************************
subroutine initialize(this, i_wm, i_t0, i_T, i_cfl, i_time, i_Pref,            &
    i_beta_star, i_lambda_prime_star, i_omega_min, i_omega_max)
!*******************************************************************************
use functions, only : linear_interp
implicit none
class(turbines_mpc_t) :: this
type(wake_model_t), intent(in) :: i_wm
real(rprec), dimension(:), intent(in) :: i_time, i_Pref
real(rprec), intent(in) :: i_t0, i_T, i_cfl, i_beta_star, i_lambda_prime_star
real(rprec), intent(in), optional :: i_omega_min, i_omega_max
integer :: i

! number of turbine rows
this%N = i_wm%N

! Wake model
this%iw = i_wm
call this%iw%makeDimensional
this%w = this%iw

! Adjoint wake model
this%iwstar = wake_model_adjoint_t(this%w%sx, this%w%sy, this%w%U_infty,                        &
    this%w%Delta, this%w%k, this%w%Dia, this%w%rho, this%w%inertia, this%w%Nx,                  &
    this%w%Ny, this%w%Ctp_spline, this%w%Cpp_spline, this%w%torque_gain(1))
this%wstar = this%iwstar

! Create time for the time horizon
this%cfl = i_cfl
this%dt = this%cfl * this%w%dx / this%w%U_infty
this%Nt = ceiling(i_T / this%dt)
allocate( this%t(this%Nt) )
do i = 1, this%Nt
    this%t(i) = i_t0 + this%dt * (i - 1)
end do

! Penalty terms away from optimal
this%beta_star = i_beta_star
this%lambda_prime_star = i_lambda_prime_star

! Set bounds of rotational speed
if (present(i_omega_min)) this%omega_min = i_omega_min
if (present(i_omega_max)) this%omega_max = i_omega_max

! allocate and set initial conditions for control variables
allocate( this%beta(this%N, this%Nt) )
allocate( this%torque_gain(this%N, this%Nt) )
allocate( this%grad_beta(this%N, this%Nt) )
allocate( this%grad_torque_gain(this%N, this%Nt) )
allocate( this%fdgrad_beta(this%N, this%Nt) )
allocate( this%fdgrad_torque_gain(this%N, this%Nt) )
this%beta(:,1) = this%iw%beta
this%torque_gain(:,1) = this%iw%torque_gain

! Interpolate the power signals and set initial condition for Pfarm
allocate( this%Pref(this%Nt) )
allocate( this%Pfarm(this%Nt) )
this%Pref = linear_interp(i_time, i_Pref, this%t)
this%Pfarm = 0._rprec
this%Pfarm(1) = sum(this%w%Phat)

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
    this%torque_gain = this%torque_gain / this%w%TORQUE / this%w%TIME**2
    this%torque_gain_min = this%torque_gain_min / this%w%TORQUE / this%w%TIME**2
    this%torque_gain_max = this%torque_gain_max / this%w%TORQUE / this%w%TIME**2
    this%omega_min = this%omega_min * this%w%TIME
    this%omega_max = this%omega_max * this%w%TIME
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
    this%torque_gain = this%torque_gain * this%w%TORQUE * this%w%TIME**2
    this%torque_gain_min = this%torque_gain_min * this%w%TORQUE * this%w%TIME**2
    this%torque_gain_max = this%torque_gain_max * this%w%TORQUE * this%w%TIME**2
    this%omega_min = this%omega_min / this%w%TIME
    this%omega_max = this%omega_max / this%w%TIME
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
real(rprec), dimension(:,:,:), allocatable, save :: du
real(rprec), dimension(:,:), allocatable, save :: Udu, Uw, Wdu, Wu, Ww
real(rprec), dimension(:,:), allocatable, save :: Bdu, Bu, Bw, Adu, Au, Aw
real(rprec), dimension(:,:), allocatable, save :: Uj, Wj
real(rprec), dimension(:), allocatable, save :: dCt_dbeta, dCt_dlambda
real(rprec), dimension(:), allocatable, save :: dCp_dbeta, dCp_dlambda
real(rprec), dimension(:), allocatable, save :: reg
real(rprec) :: dummy

! allocate adjoint forcing terms
if (.not.allocated(du)) then
    allocate(du(this%Nt, this%N, this%w%Nx))
    allocate(Udu(this%Nt,this%N))
    allocate(Uw(this%Nt,this%N))
    allocate(Wdu(this%Nt,this%N))
    allocate(Wu(this%Nt,this%N))
    allocate(Ww(this%Nt,this%N))
    allocate(Bdu(this%Nt,this%N))
    allocate(Bu(this%Nt,this%N))
    allocate(Bw(this%Nt,this%N))
    allocate(Adu(this%Nt,this%N))
    allocate(Au(this%Nt,this%N))
    allocate(Aw(this%Nt,this%N))
    allocate(Uj(this%Nt,this%N))
    allocate(Wj(this%Nt,this%N))
    allocate(dCt_dbeta(this%N))
    allocate(dCt_dlambda(this%N))
    allocate(dCp_dbeta(this%N))
    allocate(dCp_dlambda(this%N))
    allocate(reg(this%N))
end if
du = 0._rprec
Udu = 0._rprec
Uw = 0._rprec
Wdu = 0._rprec
Wu = 0._rprec
Ww = 0._rprec
Bdu = 0._rprec
Bu = 0._rprec
Bw = 0._rprec
Adu = 0._rprec
Au = 0._rprec
Aw = 0._rprec
Uj = 0._rprec
Wj = 0._rprec

! reset costs and gradients
this%cost = 0._rprec
this%grad_beta = 0._rprec
this%grad_torque_gain = 0._rprec

! Run forward in time
this%w = this%iw
do k = 2, this%Nt
    ! advance and save values that are independent of the cost function
    call this%w%adjoint_advance(this%dt, this%beta(:,k), this%torque_gain(:,k),&
        Udu(k,:), Uw(k,:), Wdu(k,:), Wu(k,:), Ww(k,:),                         &
        Adu(k,:), Au(k,:), Aw(k,:), Bdu(k,:), Bu(k,:), Bw(k,:),                &
        dCt_dbeta, dCt_dlambda, dCp_dbeta, dCp_dlambda)
    du(k,:,:) = this%w%du
    ! calculate contribution to cost function
    this%Pfarm(k) = sum(this%w%Phat)
    this%cost = this%cost + this%dt * (this%Pfarm(k) - this%Pref(k))**2
    ! Add regularizations to cost function
    reg = 0._rprec
    where (this%w%omega(:) < this%omega_min)
        reg(:) = reg(:) + this%w%omega(:) - this%omega_min
    end where
    where (this%w%omega(:) > this%omega_max)
        reg(:) = reg(:) + this%w%omega(:) - this%omega_max
    end where
    this%cost = this%cost + this%dt * this%eta1 * sum(reg**2)!                  &
!         + this%dt * this%eta2 * sum( (this%w%beta - this%beta_star)**2 )
    ! calculate adjoint values that depend on cost function
    Uj(k,:) = 0._rprec
    Wj(k,:) = -6._rprec * (this%Pfarm(k) - this%Pref(k))                       &
        * this%w%torque_gain * this%w%omega**2 - 2._rprec * this%eta1 * reg
    ! calculate part of gradients that depends only on the cost function and   &
    ! states
    this%grad_torque_gain(:,k) = 2._rprec * (this%Pfarm(k) - this%Pref(k))     &
        * this%w%omega**3 * this%dt
!     this%grad_beta(:,k) = 2._rprec * this%eta2 * (this%w%beta - this%beta_star)&
!         * this%dt
end do

! Run backwards in time
this%wstar = this%iwstar
do k = this%Nt-1, 1, -1
    ! retract adjoint wake model. Indices should definitely should be k+1
    call this%wstar%retract(du(k+1,:,:), Udu(k+1,:), Uw(k+1,:), Uj(k+1,:),     &
        Wdu(k+1,:), Ww(k+1,:), Wj(k+1,:), this%dt)
    do i = 1, this%N
        dummy = sum(this%wstar%du_star(i,:) * this%w%f(i,:)) * this%w%dx
        this%grad_beta(i,k) = this%grad_beta(i,k)                              &
            + Bw(k,i) * this%wstar%omega_star(i) * this%dt                     &
            + Bdu(k,i) * dummy * this%dt                                       &
            + Bu(k,i) * this%wstar%uhat_star(i) * this%dt
        this%grad_torque_gain(i,k) = this%grad_torque_gain(i,k)                &
            + Aw(k,i) * this%wstar%omega_star(i) * this%dt                     &
            + Adu(k,i) * dummy * this%dt                                       &
            + Au(k,i) * this%wstar%uhat_star(i) * this%dt
    end do
end do

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
    istart = (k-1) * this%N + 1
    istop = this%N * k
    this%beta(:,k+1) = x(istart:istop) * this%Cb
    this%torque_gain(:,k+1) = x(istart+iskip:istop+iskip) * this%Ca
end do

! Run model
call this%run

! Return cost function
f = this%cost

! Return gradient as vector
g = 0._rprec
do k = 1, this%Nt-1
    istart = (k-1) * this%N + 1
    istop = this%N * k
    g(istart:istop) = this%grad_beta(:,k+1) * this%Cb
    g(istart+iskip:istop+iskip) = this%grad_torque_gain(:,k+1) * this%Ca
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
    istart = (k-1) * this%N + 1
    istop = this%N * k
    x(istart:istop) = this%beta(:,k+1) / this%Cb
    x(istart+iskip:istop+iskip) = this%torque_gain(:,k+1) / this%Ca
end do

end function get_control_vector

!*******************************************************************************
function get_lower_bound(this) result(x)
!*******************************************************************************
implicit none
class(turbines_mpc_t) :: this
real(rprec), dimension(:), allocatable :: x
integer :: k, istart, istop, iskip

allocate( x(2*(size(this%beta) - this%N)) )
! combined with get_upper_bound default, would cause an error
x = 1._rprec

iskip = (this%Nt-1) * this%N
do k = 1, this%Nt-1
    istart = (k-1) * this%N + 1
    istop = this%N * k
    x(istart:istop) = this%beta_min / this%Cb
    x(istart+iskip:istop+iskip) = this%torque_gain_min / this%Ca
end do

end function get_lower_bound

!*******************************************************************************
function get_upper_bound(this) result(x)
!*******************************************************************************
implicit none
class(turbines_mpc_t) :: this
real(rprec), dimension(:), allocatable :: x
integer :: k, istart, istop, iskip

allocate( x(2*(size(this%beta) - this%N)) )
! combined with get_lower_bound default, would cause an errors
x = -1._rprec

iskip = (this%Nt-1) * this%N
do k = 1, this%Nt-1
    istart = (k-1) * this%N + 1
    istop = this%N * k
    x(istart:istop) = this%beta_max / this%Cb
    x(istart+iskip:istop+iskip) = this%torque_gain_max / this%Ca
end do

end function get_upper_bound


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

! calculate gradient for torque gain
this%fdgrad_torque_gain = 0._rprec
do n = 1, this%N
    do k = 1, this%Nt
        mf = this
        mf%torque_gain(n,k) = mf%torque_gain(n,k) + dphi
        call mf%run
        this%fdgrad_torque_gain(n,k) = (mf%cost - this%cost) / dphi
    end do
end do

end subroutine finite_difference_gradient

!*******************************************************************************
subroutine rescale_gradient(this, Ca, Cb)
!*******************************************************************************
implicit none
class(turbines_mpc_t), intent(inout) :: this
real(rprec), intent(in), optional :: Ca, Cb

call this%run()
if ( present(Ca) ) then
    this%Ca = Ca
else
    this%Ca = 1._rprec / maxval(abs(this%grad_torque_gain))
end if
if ( present(Cb) ) then
    this%Cb = Cb
else
    this%Cb = 1._rprec / maxval(abs(this%grad_beta))
end if

end subroutine rescale_gradient

end module turbines_mpc
