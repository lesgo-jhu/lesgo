!!
!!  Copyright (C) 2009-2016  Johns Hopkins University
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

module rh_control

use types, only : rprec
use param, only : pi
use wake_model
use wake_model_adjoint
use minimize_m
use messages
implicit none

private
public rh_control_t

type, extends(minimize_t) :: rh_control_t
    type(wake_model_t)         :: w       ! wake (turbine)
    type(wake_model_t)         :: iw      ! wake initial condition (turbine)
    type(wake_model_adjoint_t)  :: wstar   ! adjoint wake (turbine)
    type(wake_model_adjoint_t)  :: iwstar  ! adjoint wake initial condition  (turbine)
    real(rprec), dimension(:,:), allocatable     :: Ctp     ! local thrust coefficient (turbine, time)
    real(rprec), dimension(:,:), allocatable     :: phi     ! control variable coefficient (turbine, time)
    real(rprec), dimension(:), allocatable       :: t       ! times
    real(rprec), dimension(:), allocatable       :: Pref    ! reference signal
    real(rprec), dimension(:), allocatable       :: Pfarm   ! farm power
    real(rprec)                                  :: cost = 0._rprec    ! cost
    ! Gradient is indexed with time, but only 2:end counts since Ctp is fixed at 1
    real(rprec), dimension(:,:), allocatable     :: grad    ! gradient (turbine, time)
    real(rprec), dimension(:,:), allocatable     :: fdgrad  ! finite-difference gradient (turbine, time)
    real(rprec) :: Ctp0 = 1.33_rprec
    real(rprec) :: POWER = 0._rprec
    real(rprec) :: cfl, dt              ! Constant time step
    integer     :: N, Nt
    logical     :: isDimensionless = .false.
    real(rprec) :: tau = 120._rprec
contains
    procedure, public  :: eval
!    procedure, public  :: get_Ctp_vector
    procedure, public  :: get_phi_vector
    procedure, public  :: get_lower_bound
    procedure, public  :: get_upper_bound
    procedure, public  :: initialize
    procedure, public  :: makeDimensionless
    procedure, public  :: makeDimensional
    procedure, public  :: finiteDifferenceGradient
    procedure, public  :: run => run_input
    procedure, private :: run_noinput
end type rh_control_t

interface rh_control_t
    module procedure :: constructor
end interface rh_control_t

contains

function constructor(i_wm, i_t0, i_T, i_cfl, i_time, i_Pref, i_tau) result(this)
    implicit none
    type(rh_control_t)                         :: this
    class(wake_model_t), intent(in)                :: i_wm
    real(rprec), dimension(:), intent(in)       :: i_time, i_Pref
    real(rprec), intent(in)                     :: i_t0, i_T, i_cfl, i_tau

    call this%initialize(i_wm, i_t0, i_T, i_cfl, i_time, i_Pref, i_tau)
end function constructor

subroutine initialize(this, i_wm, i_t0, i_T, i_cfl, i_time, i_Pref, i_tau)
    use functions, only : linear_interp
    implicit none
    class(rh_control_t)                        :: this
    type(wake_model_t), intent(in)                 :: i_wm
    real(rprec), dimension(:), intent(in)       :: i_time, i_Pref
    real(rprec), intent(in)                     :: i_t0, i_T, i_cfl, i_tau
    integer                                     :: i

    ! Set reference values
    this%tau = i_tau

    ! Set initial condition for wake model
    this%iw = i_wm
    call this%iw%makeDimensional
    this%w = this%iw
    this%N = this%w%N

    ! Create adjoint wake models
    this%iwstar = wake_model_adjoint_t(this%w%s, this%w%U_infty, this%w%Delta, this%w%k, this%w%Dia, this%w%Nx)
    this%wstar = this%iwstar

    ! Dimension for power
    this%POWER = this%w%VELOCITY**3

    ! Create time for the time horizon
    this%cfl      = i_cfl
    this%dt       = this%cfl * this%w%dx / this%w%U_infty
    this%Nt       = ceiling(i_T / this%dt)
    allocate( this%t(this%Nt) )
    do i = 1, this%Nt
        this%t(i) = i_t0 + this%dt * (i - 1)
    end do

    ! Interpolate the reference signal
    allocate( this%Pref(this%Nt) )
    allocate( this%Pfarm(this%Nt) )
    this%Pref = linear_interp(i_time, i_Pref, this%t)

    ! Allocate other variables
    allocate( this%Ctp(this%N, this%Nt)    )
    allocate( this%phi(this%N, this%Nt)    )
    allocate( this%grad(this%N, this%Nt)   )
    allocate( this%fdgrad(this%N, this%Nt) )

    ! Put state of wake model Ctp into first array
    do i = 1, this%Nt
        this%Ctp(:, i) = this%iw%Ctp
    end do

end subroutine initialize

subroutine run_input(this, i_t, i_phi)
    use functions, only : linear_interp
    implicit none
    class(rh_control_t), intent(inout)        :: this
    real(rprec), dimension(:), intent(in)      :: i_t
    real(rprec), dimension(:,:), intent(in)    :: i_phi
    integer                                    :: n

    ! Make dimensionless
    call this%makeDimensionless

    ! Interpolate input onto object
    do n = 1, this%N
        this%phi(n,:) = linear_interp(i_t, i_phi(n,:), this%t * this%w%TIME)
    end do

    call this%run_noinput
end subroutine run_input

subroutine run_noinput(this)
    use functions, only : linear_interp
    implicit none
    class(rh_control_t), intent(inout)        :: this
    integer                                    :: i, n, k
    real(rprec), dimension(:), allocatable     :: du_super
    real(rprec), dimension(:), allocatable     :: uhatstar, ustar ! adjoint forcing terms
    real(rprec), dimension(:,:,:), allocatable :: fstar           ! adjoint forcing (turbine, time, space)
    real(rprec), dimension(:,:), allocatable   :: phistar         ! adjoint forcing (turbine, time)
    real(rprec), dimension(:), allocatable     :: g

    ! Allocate some variables
    allocate( du_super(this%w%Nx) )
    allocate( uhatstar(this%N) )
    allocate( ustar(this%w%Nx) )
    allocate( fstar(this%N, this%Nt, this%w%Nx) )
    allocate( phistar(this%N, this%Nt) )

    ! Make dimensionless
    call this%makeDimensionless

    ! Reset cost and gradient
    this%cost = 0._rprec
    this%grad = 0._rprec
    this%Pfarm = 0._rprec
    fstar = 0._rprec

    ! Assign initial conditions
    this%w = this%iw
    this%wstar = this%iwstar
    this%Ctp(:,1) = this%w%Ctp

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Forward Simulation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Iterate in time
    do k = 2, this%Nt
        ! Calculate next thrust coefficient
        do n = 1, this%N
            this%Ctp(n,k) = this%Ctp(n,k-1) + this%dt * ( this%phi(n,k-1) - this%Ctp(n,k-1) ) / this%tau
        end do

        ! Calculate next step
        call this%w%advance(this%Ctp(:,k), this%dt)

        ! Calculate farm power
        this%Pfarm(k) = sum(this%w%Phat)

        ! Calculate contribution to cost function
        this%cost = this%cost + (this%Pfarm(k) - this%Pref(k))**2 * this%dt

        ! Calculate the contribution dJ/dCt'_n to phistar at time k
        do n = 1, this%N
            phistar(n, k) = 2.d0 * (this%Pfarm(k) - this%Pref(k)) * this%w%uhat(n)**3 * this%dt
        end do

        ! Calculate adjoint forcing
        uhatstar = -6.d0 * (this%Pfarm(k) - this%Pref(k)) * this%Ctp(:,k) * this%w%uhat**2 * this%dt
        ustar = 0
        do n = 1, this%N
            ustar = ustar + this%w%G(n,:) * uhatstar(n)
        end do
        du_super = this%w%U_infty - this%w%u
        do n = 1, this%N
            fstar(n, k, :) = - ustar * this%w%du(n,:) / du_super
            do i = 1, this%w%Nx
                if ( du_super(i) <= 1E-10 )   fstar(n, k, i) = 0._rprec
            end do
        end do

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Backward (Adjoint) Simulation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate( g(this%N) )
    do k = this%Nt-1, 1, -1
        ! Integrate Ctstar (which is held in grad)
        do n = 1, this%N
            this%grad(n, k) = this%grad(n, k+1) - this%dt * ( this%grad(n, k+1) / this%tau + phistar(n, k+1) )
        end do

        ! Integrate adjoint wake model
        g = 0._rprec
        call this%wstar%retract(fstar(:,k+1,:), this%dt, g)

        ! Add additional term to phistar
        do n = 1, this%N
            phistar(n, k) = phistar(n, k) - g(n) * 8.d0 * this%w%U_infty**2 &
            / ( ( 4.d0 + this%Ctp(n,k))**2 )
        end do
    end do

    ! Multiply Ctstar to get gradient
    this%grad = -this%grad / this%tau

end subroutine run_noinput

subroutine makeDimensionless(this)
    implicit none
    class(rh_control_t), intent(inout)  :: this

    if (.not.this%isDimensionless) then
        this%isDimensionless = .true.
        this%Pref  = this%Pref / this%POWER
        this%Pfarm = this%Pfarm / this%POWER
        this%t     = this%t    / this%w%TIME
        this%dt    = this%dt   / this%w%TIME
        this%tau   = this%tau / this%w%TIME
        call this%w%makeDimensionless
        call this%wstar%makeDimensionless
        call this%iw%makeDimensionless
        call this%iwstar%makeDimensionless
    end if
end subroutine makeDimensionless

subroutine makeDimensional(this)
    implicit none
    class(rh_control_t), intent(inout)  :: this

    if (this%isDimensionless) then
        this%isDimensionless = .false.
        this%Pref  = this%Pref * this%POWER
        this%Pfarm = this%Pfarm * this%POWER
        this%t     = this%t    * this%w%TIME
        this%dt    = this%dt   * this%w%TIME
        this%tau   = this%tau * this%w%TIME
        call this%w%makeDimensional
        call this%wstar%makeDimensional
        call this%iw%makeDimensional
        call this%iwstar%makeDimensional
    end if
end subroutine makeDimensional

subroutine eval(this, x, f, g)
    implicit none
    class(rh_control_t), intent(inout)                   :: this
    real(rprec), dimension(:), intent(in)                 :: x
    real(rprec), intent(inout)                            :: f
    real(rprec), dimension(:),  intent(inout)             :: g
    integer                         :: k

    ! Place x in this%phi
    do k = 1, this%Nt-1
        this%phi(:,k) = x((k - 1) * this%N + 1 : this%N * k)
    end do

    ! Run model
    call this%run_noinput

    ! Return cost function
    f = this%cost

    ! Return gradient as vector
    g = 0.d0
    do k = 1, this%Nt-1
        g((k - 1) * this%N + 1 : this%N * k) = this%grad(:,k + 1)
    end do
end subroutine eval

subroutine finiteDifferenceGradient(this)
    implicit none
    class(rh_control_t), intent(inout)     :: this
    type(rh_control_t)                     :: mf
    integer :: n, k
    real(rprec) :: dphi

    this%fdgrad = 0.d0
    dphi = sqrt( epsilon( this%Ctp(n,k) ) )
    do n = 1, this%N
        do k = 1, this%Nt
            mf = this
            mf%phi(n,k) = mf%phi(n,k) + dphi
            call mf%run_noinput
            this%fdgrad(n,k) = (mf%cost - this%cost) / dphi
        end do
    end do
end subroutine finiteDifferenceGradient

function get_phi_vector(this) result(phi_vec)
    implicit none
    class(rh_control_t) :: this
    real(rprec), dimension(:), allocatable :: phi_vec
    integer :: k

    allocate(phi_vec(size(this%phi) - this%N))

    do k = 1, this%Nt-1
        phi_vec((k - 1) * this%N + 1 : this%N * k) = this%phi(:,k)
    end do

end function get_phi_vector

function get_lower_bound(this) result(lb)
    implicit none
    class(rh_control_t) :: this
    real(rprec), dimension(:), allocatable :: lb

    allocate(lb(size(this%phi) - this%N))
    lb = 0._rprec

end function get_lower_bound

function get_upper_bound(this) result(ub)
    implicit none
    class(rh_control_t) :: this
    real(rprec), dimension(:), allocatable :: ub

    allocate(ub(size(this%phi) - this%N))
    ub = 2._rprec

end function get_upper_bound



! function get_Ctp_vector(this) result(Ctp_vec)
!     implicit none
!     class(rh_control_t) :: this
!     real(rprec), dimension(:), allocatable :: Ctp_vec
!     integer :: k
!
!     allocate(Ctp_vec(size(this%Ctp) - this%N))
!
!     do k = 1, this%Nt-1
!         Ctp_vec((k - 1) * this%N + 1 : this%N * k) = this%Ctp(:,k)
!     end do
!
! end function get_Ctp_vector

end module rh_control
